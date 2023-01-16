#!/usr/bin/env python
## -*- python -*-
"""PGDB Workflow

Usage:
	pgdb_build_wf.py --mp_out <mp_dir> --sif <sif_file> --tmp_dir <tmp_dir> [--taxprune <taxprune> ] [--tag <tag>]

Options:
	-h --help	Show this screen.
	--version	Show version.
	--mp_out=DIR	MP3 output directory.
	--sif=FILE	Path to the Ptools SIF file.
	--tmp_dir=DIR	TMP working dir for Ptools to save intermediates.
	--taxprune	Use taxonomic pruning when building PGDBs [True or False; default: False]
	--tag=STR	Tag for metagenome PGDB [default: community].
"""


import sys
import pandas as pd
import os
import numpy as np
import subprocess
import shutil
import glob
from pathlib import Path
from sexpdata import loads, dumps, Symbol
from docopt import docopt
from camelot_frs.camelot_frs import get_kb, get_frame, get_frame_all_children, frame_parent_of_frame_p, frame_object_p
from camelot_frs.pgdb_loader import load_pgdb, make_camelot_file
from camelot_frs.pgdb_api import genes_of_pathway
import html2text



def create_pgdb(pt_inputs, pt_outputs, sif_file,
				tmp_dir,  tprune, tag
				):

	rename_pgdb(pt_inputs, tag)

	# Create output dir if doesn't exist
	Path(pt_outputs).mkdir(parents=True, exist_ok=True)

	# Get $PATH
	my_path = os.environ.copy()['PATH']

	# Run Singularity images for Ptools
	bind_str = ''.join([pt_inputs, ':/pt_inputs,',
						pt_outputs, ':/pt_outputs,',
						tmp_dir, ':/data'
						])
	if tprune == 'True':
		pt_cmd = ['singularity', 'run', '--env', 'APPEND_PATH=' + my_path,
					'-B', bind_str, sif_file,
					'run-pathway-tools-and-copy-pgdb-singularity_taxprune.sh',
					'/pt_inputs', '/pt_outputs'
					]
	elif tprune == 'False':
		pt_cmd = ['singularity', 'run', '--env', 'APPEND_PATH=' + my_path,
					'-B', bind_str, sif_file,
					'run-pathway-tools-and-copy-pgdb-singularity.sh',
					'/pt_inputs', '/pt_outputs'
					]
	pt_out = subprocess.run(' '.join(pt_cmd),
							shell=True
							) #, capture_output=True, text=True).stdout
	
	# Uncompress PGDB to create PWYs table
	pgdb_arc = glob.glob(pt_outputs + '/*.tar.bz2')[0]
	tar_cmd = ['tar', '-xf', pgdb_arc, '-C', pt_outputs]
	tar_out = subprocess.run(tar_cmd) #, capture_output=True, text=True).stdout


def rename_pgdb(pt_inputs, tag):
	o_params = os.path.join(pt_inputs, 'organism-params.dat')
	attr_list = ['ID', 'NAME', 'ABBREV-NAME']
	with open(o_params, 'r') as in_par:
		data = in_par.readlines()
		with open(o_params + '.tmp', 'w') as out_par:
			for line in data:
				line = line.strip('\n')
				split_line = line.split('\t')
				if split_line[0] in attr_list:
					split_line[1] = tag
				new_line = '\t'.join(split_line) + '\n'
				out_par.write(new_line)
	os.rename(o_params + '.tmp', o_params)


def extract_pwy(pt_outputs):
	## version.dat file is not in expected directory, create it
	pt_id = os.path.basename(glob.glob(pt_outputs + '/*.tar.bz2')[0]).split('cyc', 1)[0]
	flatpath = os.path.join(pt_outputs, '1.0/data')
	pwy_outfile = os.path.join(pt_outputs, pt_id + '_pwy.tsv')
	verfile = os.path.join(flatpath.rsplit('/', 2)[0], 'default-version')
	new_verfile = os.path.join(flatpath, 'version.dat')
	shutil.copyfile(verfile, new_verfile)

	## Need to create a custom sample_id since org_id is blank
	
	## Create the .camelot file:
	org_id = make_camelot_file(flatpath, pt_outputs)
	print('SAMPLE_ID:', pt_id)
	print('ORGANISM_ID:', org_id)

	## Load the PGDB:
	load_pgdb(pt_outputs + '/' + org_id + '.camelot')

	curr_kb = get_kb(org_id)

	# Build Pathway Inference Data Dictionary from contents of ./reports/ dirextory
	reportspath = os.path.dirname(flatpath) + '/reports'
	pwy_inf_data = get_pwy_inf(reportspath)

	# Build Pathway/Superpathway Data Dictionary for reactions:
	pwy_evi_data = get_superpath_rxns(reportspath)
	## Generate the report:
	headers = [ "SAMPLE",
				"PWY_NAME",
				"PWY_COMMON_NAME",
				"PWY_SCORE",
				"NUM_REACTIONS",
				"NUM_COVERED_REACTIONS",
				"ORF_COUNT",
				"ORFS" 
			   ]

	with open(pwy_outfile,"w") as report_fp:

		print('\t'.join(headers),
			  file=report_fp)
		for pwy in get_frame_all_children(get_frame(curr_kb, 'Pathways'), frame_types='instance'):
			# if you wanna check for Super-Pathways: frame_parent_of_frame_p(get_frame(curr_kb, 'Super-Pathways'), pwy)
			#try:
			rxns_covered = len(pwy_evi_data[pwy.frame_id]['RXNs'])
	
			rxns_total = 0
			end = False
			rxns = pwy.get_slot_values('REACTION-LIST')
			for rxn in rxns:
				r_slots = rxn.slots
				if 'REACTION-LIST' in r_slots:
					rrxns = rxn.get_slot_values('REACTION-LIST')
					for rrxn in rrxns:
						rr_slots = rrxn.slots
						if 'REACTION-LIST' in rr_slots:
							rrrxns = rrxn.get_slot_values('REACTION-LIST')
							for rrrxn in rrrxns:
								rxns_total += len(rrrxn.get_slot_values('REACTION-LIST'))	
						else:
							rxns_total += 1
				else:
					rxns_total += 1

			if pwy.frame_id in pwy_inf_data:
				pscore = pwy_inf_data[pwy.frame_id]['SCORE']
			else:
				pscore = 'NaN'
			print(pwy, pscore, rxns_total, rxns_covered)
			try:
				pwy_gene_names = [ str(gene.get_slot_values('COMMON-NAME')[0]).lstrip('frame:') for gene in genes_of_pathway(pwy) ]
			except Exception:
				pwy_genes_names = ['EcoCyc','error']
			print('\t'.join([pt_id, #  curr_kb.kb_name,
							 pwy.frame_id,
							 html2text.html2text(pwy.get_slot_values('COMMON-NAME')[0]).split('\n')[0],
							 pscore,
							 str(rxns_total),
							 str(rxns_covered),
							 str(len(pwy_gene_names)),
							 ','.join(pwy_gene_names)]),
				  file = report_fp)
			#except:
			#	print('Warning: ', pwy, ' does not have reaction infomation available...')
			#	print(pwy.get_slot_values('COMMON-NAME')[0])
			#	print(pwy.get_slot_values('REACTION-LIST'))
			#	for rxn in pwy.get_slot_values('REACTION-LIST'):
			#		print(rxn.slots)
			#	print(pwy.slots)
def get_pwy_inf(reports_dir):
	"""
	Accepts the path to the 'reports' directory within
	Ptools flatfile output.

	Returns dictionary of all values found in
	'pwy-inference-report_YYYY-MM-DD.txt' file.
	"""
	pwy_inf_rec_list = []
	pwy_inf_file = glob.glob(os.path.join(reports_dir, 'pwy-inference-report_*.txt'))[0]

	with open(pwy_inf_file, 'r') as pwy_inf_in:

		data = pwy_inf_in.read()
		trim_dat = data.split('::: Pathway Inference Report')
		if len(trim_dat) == 3:
			keep_dat = trim_dat[2]
		else:
			keep_dat = trim_dat[1]
		keep_dat = keep_dat.split('List of pathways pruned')[0]

		pwy_inf_rec = ''
		start = False
		for line in keep_dat.split('\n'):
			if line[:2] == ' (': # start of record
				if pwy_inf_rec != '': # add if there is something to add
					pwy_inf_rec_list.append(pwy_inf_rec)
					pwy_inf_rec = line
				else: # start a new record
					pwy_inf_rec = line
				start = True
			elif start == True:
				pwy_inf_rec = pwy_inf_rec + line
		pwy_inf_rec_list.append(pwy_inf_rec) # add last record
	pwy_inf_dict = {}
	for p_rec in pwy_inf_rec_list:
		parsed_sexpr = [r.value() if isinstance(r, Symbol) else str(r) for r in loads(p_rec)]
		pwy_id = parsed_sexpr[0]
		pwy_conf = parsed_sexpr[2]
		pwy_score = parsed_sexpr[5]
		pwy_inf_dict[pwy_id] = {'SCORE': pwy_score, 'CONFIDENCE': pwy_conf}

	return pwy_inf_dict


def get_superpath_rxns(reports_dir):
	"""
	Accepts the path to the 'reports' directory within
	Ptools flatfile output.

	Returns dictionary of all values found in
	'pwy-evidence-list.dat' file.
	"""
	pwy_evi_rec_list = []
	pwy_evi_file = glob.glob(os.path.join(reports_dir, 'pwy-evidence-list.dat'))[0]

	with open(pwy_evi_file, 'r') as pwy_evi_in:
		data = pwy_evi_in.read()
		for line in data.split('\n'):
			if ';;;' not in line[:3]: # start of record
				pwy_evi_rec_list.append(line)
	pwy_evi_dict = {}
	for e_rec in pwy_evi_rec_list:
		if e_rec:
			parsed_sexpr = [r.value() if isinstance(r, Symbol) else str(r) for r in loads(e_rec)]
			pwy_id = parsed_sexpr[0]
			rxns = parsed_sexpr[1:]
			pwy_evi_dict[pwy_id] = {'RXNs': rxns}

	return pwy_evi_dict


def get_present_rxns(pwy_frame, pwy_evi_dict):

	pwy_expl = loads(pwy_frame.get_slot_values('EXPLANATION-CODE')[0])
	pwy_rxns = {}
	for r in pwy_expl:
		if isinstance(r, Symbol):
			r = r.value()
		elif isinstance(r, list):
			for rr in r:
				if isinstance(rr, Symbol):
					rr = rr.value()
					rr_key = rr
				elif isinstance(rr, list):
					rrr_vals = []			
					for rrr in rr:
						if isinstance(rrr, Symbol):
							rrr = rrr.value()
							rrr_vals.append(rrr)
						else:
							rrr = str(rrr)
					pwy_rxns[rr_key] = rrr_vals

	return pwy_rxns


def map_orfs2pwys(mp_outdir, pt_outdir):
	pt_id = os.path.basename(glob.glob(pt_outdir + '/*.tar.bz2')[0]).split('cyc', 1)[0]
	orf_mapfile = glob.glob(os.path.join(mp_outdir, 'results/annotation_table/*.EC_RXN_map.tsv'))[0]
	pwy_outfile = os.path.join(pt_outdir, pt_id + '_pwy.tsv')
	pwy2orf_outfile = os.path.join(pt_outdir, pt_id + '_pwy2orf.tsv')
	orf_map_df = pd.read_csv(orf_mapfile, sep='\t', header=0)
	pwy_out_df = pd.read_csv(pwy_outfile, sep='\t', header=0)
	orf_exp_list = []
	for i, row in pwy_out_df.iterrows():
		r_list = list(row)
		orfs = str(row['ORFS'])
		if ((orfs != 'nan') & (orfs != ['nan']) & (orfs != '')):
			if ',' in orfs:
				orf_list = row['ORFS'].split(',')
			else:
				orf_list = [orfs]
		else:
			orf_list = [""]
		for orf_id in orf_list:
			if orf_id:
				clean_id = orf_id.split('_', 1)[1]
				new_row = [clean_id]
				new_row.extend(r_list)
				orf_exp_list.append(new_row)
	new_cols = ['orf_id']
	new_cols.extend(pwy_out_df.columns)
	orf_exp_df = pd.DataFrame(orf_exp_list, columns=new_cols)
	pwy2orf_df = orf_map_df.merge(orf_exp_df, on='orf_id', how='outer')
	pwy2orf_df.dropna(subset=['SAMPLE'], inplace=True)
	pwy2orf_df = pwy2orf_df[['orf_id', 'SAMPLE', 'EC', 'RXN', 'PWY_COMMON_NAME', 'ref dbname',
							 'target', 'product', 'value', 'trim_target', 'PWY_NAME', 'PWY_SCORE',
							 'NUM_REACTIONS', 'NUM_COVERED_REACTIONS', 'ORF_COUNT', 'ORFS'
							 ]]
	pwy2orf_df.to_csv(pwy2orf_outfile, sep='\t', index=False)

###############################################################
# Collect inputs
arguments = docopt(__doc__, version='PGDB Workflow 1.0')

mp_dir = arguments['<mp_dir>']
sif_file = arguments['<sif_file>']
tmp_dir = arguments['<tmp_dir>']
if arguments['<taxprune>'] == None:
	taxprune = 'False'
else:
	taxprune = arguments['<taxprune>']
	
if arguments['<tag>'] == None:
	tag = 'community'
else:
	tag = arguments['<tag>']

# Build Community-level PGDB
pt_in = os.path.join(mp_dir, 'ptools')
pt_out = os.path.join(mp_dir, 'results/pgdb/community')
create_pgdb(pt_in, pt_out, sif_file, tmp_dir, taxprune, tag)
# Parse PGDB flatfiles to create PWYs TSV table
extract_pwy(pt_out)
# Map inferred pwys to ORFs and ECs/RXNs used
map_orfs2pwys(mp_dir, pt_out)


# Build MAG-level PGDBs if they exist
ms_dir = os.path.join(mp_dir, 'magsplitter/results')
if os.path.exists(ms_dir):
	mag_list = glob.glob(ms_dir + '/*')
	for pt_mag in mag_list:
		mag_id = os.path.basename(pt_mag)
		mag_tag = tag + '_' + mag_id
		pt_out = os.path.join(mp_dir, 'results/pgdb/MAGs/' + mag_id)
		create_pgdb(pt_mag, pt_out, sif_file, tmp_dir, taxprune, mag_tag)
		# Parse PGDB flatfiles to create PWYs TSV table
		extract_pwy(pt_out)
		# Map inferred pwys to ORFs and ECs/RXNs used
		map_orfs2pwys(mp_dir, pt_out)

