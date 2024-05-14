#!/usr/bin/env python
## -*- python -*-
"""Prepping EggNOG for MetaPathways

Usage:
	eggnog_db_prep.py --en_faa <faa> --en_aliase <aliase> --out_dir <dir>

Options:
	-h --help	Show this screen.
	--version	Show version.
	--en_faa=FILE	Amino Acid FASTA file for EggNOG.
	--en_aliase=FILE	TSV containing all crossrefs for EggNOG records.
	--out_dir=DIR	Path to save all outputs.
"""


import pandas as pd
import sys
import numpy as np
from docopt import docopt
import pyfastx
import codecs
from tqdm import tqdm
import os



def get_defs(faa_file):
	fa = pyfastx.Fasta(faa_file, build_index=False)
	def_list = []
	print('Parsing FASTA...')
	for name, seq in tqdm(fa):
		def_list.append(name)
	def_df = pd.DataFrame(def_list, columns=['accesion'])
	return def_df


def load_aliase(ali_file):
	with codecs.open(ali_file, 'r', encoding='utf-8',
		errors='ignore') as ali_in:
		lines = ali_in.readlines()
	
	keep_list = ['BLAST_KEGG_EC', 'RefSeq_product',
					'BLAST_UniProt_DE', 'Ensembl_UniProt_DE'
					]
	ali_list = []
	print('Filtering Aliases...')
	for line in tqdm(lines):
		for k in keep_list:
			if k in line:
				ali_list.append(line.strip('\n').split('\t'))
	ali_df = pd.DataFrame(ali_list, columns=['accesion', 'value', 'source'])
	return ali_df


def build_names_txt(ids_df, alis_df):
	egg_dict = {}
	print('Initializing Mapping Dictionary...')
	for acc in tqdm(list(ids_df['accesion'])):
		egg_dict[acc] = ['no_product', 'no_EC']
	print('Mapping Deflines to Aliases...')
	for row in tqdm(alis_df.itertuples()):
		key = row.accesion
		source = row.source
		value = row.value
		if 'RefSeq_product' in source:
			egg_dict[key][0] = value
		elif ((('BLAST_UniProt_DE' in source) | 
				('Ensembl_UniProt_DE' in source)) &
				(egg_dict[key][0] == 'no_product')
				):
			egg_dict[key][0] = value
		if 'BLAST_KEGG_EC' in source:
			value = value.replace('EC ', '')
			egg_dict[key][1] = value
	print('Creating DataFrame...')
	egg_map = []
	for k in tqdm(egg_dict):
		v = egg_dict[k]
		egg_map.append([k, v[0], v[1]])

	map_df = pd.DataFrame(egg_map)
	return map_df





###############################################################
# Run prep
arguments = docopt(__doc__, version='EggNOG Prep 1.0')

en_faa = arguments['<faa>']
en_ali = arguments['<aliase>']
out_dir = arguments['<dir>']

if '.' in en_faa:
	en_id = os.path.basename(en_faa).split('.')[0]
else:
	en_id = os.path.basename(en_faa)

sv_file = os.path.join(out_dir, en_id + '-names.txt')

print("Prepping " + en_id + " for MetaPathways...")
print("Collecting Deflines...")
egg_ids_df = get_defs(en_faa)
print("Loading Aliases...")
egg_aliase_df = load_aliase(en_ali)
print("Building Mapping...")
egg_map_df = build_names_txt(egg_ids_df, egg_aliase_df)
print("Saving Mapping to " + sv_file + "...")
egg_map_df.to_csv(sv_file, header=False,
				  index=False, sep='\t'
				  )