#!/usr/bin/env python3
## -*- python -*-
"""PGDB Workflow

Usage:
	pgdb_build_wf.py --mp_out <mp_dir> --sif <sif_file> --tmp_dir <tmp_dir>  --tag <tag> [--taxprune <taxprune> ]

Options:
	-h --help	Show this screen.
	--version	Show version.
	--mp_out=DIR	MP3 output directory.
	--sif=FILE	Path to the Ptools SIF file.
	--tmp_dir=DIR	TMP working dir for Ptools to save intermediates.
	--tag=STR	Tag for metagenome PGDB.
	--taxprune	Use taxonomic pruning when building PGDBs [True or False; default: False]
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
import argparse
from camelot_frs.camelot_frs import get_kb, get_frame, get_frame_all_children, frame_parent_of_frame_p, frame_object_p
from camelot_frs.pgdb_loader import load_pgdb, make_camelot_file
from camelot_frs.pgdb_api import genes_of_pathway
import html2text
import time


def create_pgdb(pt_inputs, pt_outputs, tprune, tag, container):

	rename_pgdb(pt_inputs, tag)
	tag_id = tag.lower()
	# Create output dir if doesn't exist
	Path(pt_outputs).mkdir(parents=True, exist_ok=True)
	if container:
		sh_tax = 'run-pathway-tools-and-copy-pgdb-taxprune.sh'
		sh_notax = 'run-pathway-tools-and-copy-pgdb.sh'
	else:
		sh_tax = 'run-pathway-tools-and-copy-pgdb-taxprune_local.sh'
		sh_notax = 'run-pathway-tools-and-copy-pgdb_local.sh'
	if tprune == True:
		pt_cmd = [sh_tax, pt_inputs, pt_outputs, tag_id]
	elif tprune == False:
		pt_cmd = [sh_notax, pt_inputs, pt_outputs, tag_id]
	os.system(' '.join(pt_cmd))
	# Uncompress PGDB to create PWYs table
	pgdb_arc = glob.glob(pt_outputs + '/*.tar.bz2')[0]
	tar_cmd = ['tar', '-xf', pgdb_arc, '-C', pt_outputs]
	tar_out = subprocess.run(tar_cmd)

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
	print(pwy_inf_data)
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
			if not frame_parent_of_frame_p(get_frame(curr_kb, 'Super-Pathways'),
										   pwy):
				
				print(pwy)
				pwy_rxn_dict = get_present_rxns(pwy)
				print(pwy_rxn_dict)
				enz_rxn_count = 0
				for rxn in pwy.get_slot_values('REACTION-LIST'):
					if 'ENZYMATIC-REACTION' in rxn.slots:
						enz_rxn_count += 1
				covered_rxn_count = len(pwy_rxn_dict['REACTIONS-PRESENT'])
				print(pwy_inf_data[pwy.frame_id])
				pscore = pwy_inf_data[pwy.frame_id]['SCORE']
				print(pwy, pscore, covered_rxn_count)

				try:
					pwy_gene_names = [ str(gene.get_slot_values('COMMON-NAME')[0]).lstrip('frame:') for gene in genes_of_pathway(pwy) ]
				except:
					pwy_gene_names = [ str(gene.get_slot_values('COMMON-NAME')) for gene in genes_of_pathway(pwy) ]
					if pwy_gene_names[0] == '[]':
						pwy_gene_names = [str(f).split(':')[1] for f in genes_of_pathway(pwy)]
				print('\t'.join([pt_id, #  curr_kb.kb_name,
								 pwy.frame_id,
								 pwy.get_slot_values('COMMON-NAME')[0],
								 pscore,
								 str(len(pwy.get_slot_values('REACTION-LIST'))),
								 str(covered_rxn_count),
								 str(len(pwy_gene_names)),
								 ','.join(pwy_gene_names)]),
					  file = report_fp)


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
	keep_dat = []
	for t_rec in trim_dat:
		if 'List of pathways pruned' in t_rec:
			k_rec = t_rec.split('List of pathways pruned')[0]
			keep_dat.append(k_rec)
	keep_dat = data #'\n'.join(keep_dat)
	pwy_inf_rec = ''
	start = False
	skip = False
	for line in keep_dat.split('\n'):
		if 'List of pathways pruned' in line:
			skip = True
		if ':::' in line:
			skip = False
		if skip == False:
			if line[:2] == ' (': # start of record
				if pwy_inf_rec != '': # add if there is something to add
					pwy_inf_rec_list.append(pwy_inf_rec)
					pwy_inf_rec = line
				else: # start a new record
					pwy_inf_rec = line
				start = True
			elif line != '':
				if ((start == True) & (line[0] == ' ')):
					pwy_inf_rec = pwy_inf_rec + line
	pwy_inf_rec_list.append(pwy_inf_rec) # add last record
	pwy_inf_dict = {}
	for p_rec in pwy_inf_rec_list:
		parsed_sexpr = [r.value() if isinstance(r, Symbol) else str(r) for r in loads(p_rec)]
		pwy_id = parsed_sexpr[0]
		pwy_pass = parsed_sexpr[1]
		pwy_conf = parsed_sexpr[2]
		pwy_score = parsed_sexpr[5]
		if pwy_pass == 'T':
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


def get_present_rxns(pwy_frame):
	
	expl_code = pwy_frame.get_slot_values('EXPLANATION-CODE')[0]
	if expl_code != 'PWY-HAS-NON-COMPUTATIONAL-EVIDENCE':
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
	else:
		pwy_rxns = {'REACTIONS-PRESENT': []}
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
				clean_id = orf_id
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
parser = argparse.ArgumentParser(description="Run Ptools on MP output and collect outputs.")
parser.add_argument("--mp_out", type=str, help="MP3 output directory.", required=True)
parser.add_argument("--tag", type=str, help="Tag for metagenome PGDB.", required=True)
parser.add_argument("--taxprune", action='store_true',
					help="Use taxonomic pruning when building PGDBs [True or False; default: False]",
					required=False
					)
parser.add_argument("--container", action='store_true', dest="container", default=False, required=False,
					help="Use when using containerized env",
					)
args = parser.parse_args()

mp_dir = args.mp_out
tag = args.tag
taxprune = args.taxprune
container = args.container


# Build Community-level PGDB
pt_in = os.path.join(mp_dir, 'ptools')
pt_out = os.path.join(mp_dir, 'results/pgdb/community')
print("Building Community-level PGDB.")
create_pgdb(pt_in, pt_out, taxprune, tag, container)
print("Completed Community-level PGDB.")
# Parse PGDB flatfiles to create PWYs TSV table
print("Extracting Community-level PGDB.")
extract_pwy(pt_out)
print("Extracting Complete.")
# Map inferred pwys to ORFs and ECs/RXNs used
print("Mapping ORFs to Inferred Pathways.")
map_orfs2pwys(mp_dir, pt_out)
print("Mapping Complete.")

# Build MAG-level PGDBs if they exist
ms_dir = os.path.join(mp_dir, 'magsplitter/results')
if os.path.exists(ms_dir):
	mag_list = glob.glob(ms_dir + '/*')
	for pt_mag in mag_list:
		if "non_binned" not in pt_mag:
			mag_id = os.path.basename(pt_mag)
			mag_tag = mag_id
			pt_out = os.path.join(mp_dir, 'results/pgdb/MAGs/' + mag_id)
			print(f"Building PGDB for {mag_id}.")
			create_pgdb(pt_mag, pt_out, taxprune, mag_tag, container)
			print("Completed PGDB.")
			# Parse PGDB flatfiles to create PWYs TSV table
			print(f"Extracting PGDB for {mag_id}.")
			extract_pwy(pt_out)
			print("Extracting Complete.")
			# Map inferred pwys to ORFs and ECs/RXNs used
			print(f"Mapping {mag_id} ORFs to Inferred Pathways.")
			map_orfs2pwys(mp_dir, pt_out)
			print("Mapping Complete.")
