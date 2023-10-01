import urllib.request
import time
import pandas as pd
import sys
from tqdm import tqdm
import pyfastx
import numpy as np

'''
uniprot_idmapping_file = 'idmapping.dat'
chunksize = 10 ** 6
df_list = []
for sub_df in tqdm(pd.read_csv(uniprot_idmapping_file, sep='\t',
						 names=['UniProtKB-AC', 'ID_type', 'ID'], chunksize=chunksize
						 )):
	# grab the KEGGs first
	sub_kegg_df = sub_df.query("ID_type == 'KEGG'")
	sub_kegg_dict = dict(zip(sub_kegg_df['UniProtKB-AC'], sub_kegg_df['ID']))
	df_list.append(sub_kegg_df)
	# Now map all other UniRefs to search out any other possible KEGGs
	k_list = list(sub_kegg_df['UniProtKB-AC'])
	uni_list = ['UniRef100', 'UniRef90', 'UniRef50']
	sub_uni_df = sub_df.query("`UniProtKB-AC` == @k_list & ID_type == @uni_list")
	sub_uni_df['new_UniProtKB-AC'] = [x.split('_', 1)[1] for x in sub_uni_df['ID']]
	sub_uni_df['new_ID_type'] = 'KEGG'
	sub_uni_df['new_ID'] = [sub_kegg_dict[x] for x in sub_uni_df['UniProtKB-AC']]
	trim_uni_df = sub_uni_df[['new_UniProtKB-AC', 'new_ID_type', 'new_ID']]
	trim_uni_df.columns = ['UniProtKB-AC', 'ID_type', 'ID']
	df_list.append(trim_uni_df)

kegg_df = pd.concat(df_list)
kegg_df.drop_duplicates(inplace=True)
kegg_df['KEGG_Prefix'] = [x.split(':')[0] for x in kegg_df['ID']]
kegg_df.to_csv('KEGG_idmapping.tsv', sep='\t', index=False)
#kegg_df = pd.read_csv('KEGG_idmapping.tsv', sep='\t', header=0)
print(kegg_df.head())
uniq_prefix_list = list(kegg_df['KEGG_Prefix'].unique())
ko2id_list = []
for prefix in tqdm(uniq_prefix_list):
	try:
		url = 'http://rest.kegg.jp/link/' + prefix + '/ko'
		fp = urllib.request.urlopen(url)
		content_list = [x.decode('utf-8').strip('\n').split('\t') for x in fp.readlines()]
		ko2id_list.extend(content_list)
	except:
		print(prefix, ' failed to download, might not exist...')

ko2id_df = pd.DataFrame(ko2id_list, columns=['KO_ID', 'ID'])
ko2id_df.to_csv('KO_idmapping.tsv', sep='\t', index=False)
print(ko2id_df.head())
kegg2ko_df = kegg_df.merge(ko2id_df, on='ID', how='left')
kegg2ko_df = kegg2ko_df.replace(np.nan, '', regex=True)
kegg2ko_df['KO_ID'] = [x.split(':')[1] if str(x) != '' else x for x in kegg2ko_df['KO_ID']]
kegg2ko_df.to_csv('KEGG2KO_idmapping.tsv', sep='\t', index=False)
'''
'''
kegg2ko_df = pd.read_csv('KEGG2KO_idmapping.tsv', sep='\t', header=0)
print(kegg2ko_df.head())
uni2ko_dict = pd.Series(kegg2ko_df['KO_ID'].values, index=kegg2ko_df['UniProtKB-AC'].values).to_dict()

ko_hier_file = 'ko00001.keg'
with open(ko_hier_file, 'r') as ko_in:
	ko_dat = ko_in.readlines()
ko_map_list = []
for i, line in enumerate(ko_dat):
	print(i, line)
	line = line.strip('\n')
	if len(line) > 1:
		if line[0] == 'A':
			a_val = line[1:].strip()
			a_code = a_val.split(' ', 1)[0]
			a_str = a_val.split(' ', 1)[1].strip()
		elif line[0] == 'B':
			b_val = line[1:].strip()
			b_code = b_val.split(' ', 1)[0]
			b_str = b_val.split(' ', 1)[1].strip()
		elif line[0] == 'C':
			c_val = line[1:].strip()
			c_code = c_val.split(' ', 1)[0]
			c_str = c_val.split(' ', 1)[1].strip()
			if 'PATH' in c_str:
				c_path = c_str.rsplit(' ', 1)[1].split(':')[1].strip(']')
			else:
				c_path = ''
		elif line[0] == 'D':
			d_val = line[1:].strip()
			if len(d_val.split(' ', 1)) == 1:
				d_code = d_val
				d_str = ''
				d_name = ''
				d_ec = ''
			else:
				d_code = d_val.split(' ', 1)[0]
				d_str = d_val.split(' ', 1)[1].strip()
				if ';' in d_str:
					d_name = d_str.split(';')[0]
				else:
					d_name = ''
				if 'EC:' in d_str:
					d_ec = d_str.rsplit(' ', 1)[1].strip('[').strip(']')
				else:
					d_ec = ''
			ko_map_list.append([i, d_code, d_name, d_ec, d_str, c_code, c_path, c_str,
				b_code, b_str, a_code, a_str
				])

ko_hier_df = pd.DataFrame(ko_map_list, columns=['unique_ID','KO', 'name', 'EC', '4_header', '3_ID',
	'Pathway', '3_header', '2_ID', '2_header', '1_ID', '1_header']
	)
ko_hier_df.to_csv('KO_hierarchy.tsv', sep='\t', index=False)

kegg2ko_df = pd.read_csv('KEGG2KO_idmapping.tsv', sep='\t', header=0)
ko_hier_df = pd.read_csv('KO_hierarchy.tsv', sep='\t', header=0)

ko_final_map_df = ko_hier_df.merge(kegg2ko_df, left_on='KO', right_on='KO_ID', how='left')[['KO', 'EC', 'UniProtKB-AC']].drop_duplicates()
ko_final_map_df.columns = ['KO', 'EC', 'UniProtKB']
ko_final_map_df.to_csv('Uniprot_KEGG_EC_map.tsv', sep='\t', index=False)

flurp
# MetaCyc
enzrxns_file = 'metacyc-26_enzrxns.dat'
with open(enzrxns_file, 'r', encoding='unicode_escape') as er_in:
    data = er_in.read().split('//')
enzrxn_list = []
for d in tqdm(data):
        d_split = d.split('\n')
        tmp_dict = {}
        skip_orf = True
        for l in d_split:
            l = l.replace('\n', '')
            attr = l.split(' - ', 1)[0]
            if l != '' and l[0] != '#':
                if attr == "ENZYME":
                    enz_id = l.split(' - ', 1)[1]
                elif attr == "REACTION":
                    rxn_id = l.split(' - ', 1)[1]
                    enzrxn_list.append([enz_id, rxn_id])
enzrxns_df = pd.DataFrame(enzrxn_list, columns=['MC', 'RXN'])
enzrxns_df.to_csv('MetaCyc-26-monomer-rxn-pairs.tsv', sep='\t', index=False)
'''
