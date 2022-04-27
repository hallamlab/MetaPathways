# This script takes a '0.pf' file as input
# The 0.pf is automatically output by MP within the ./ptools dir
# It is meant to deduplicate the ORFs based on exact match to FUNCTION
import sys
import os
import pandas as pd
pd.set_option('display.max_columns', None)


pf_file = sys.argv[1]  # the 0.pf file
orf_file = sys.argv[2] # ORF annotation file from MP3 results dir
ko_file = sys.argv[3] # Mapping file from Uniprot -> KO -> EC
mc_file = sys.argv[4] # MetaCyc Parsed Blast results table
rxn_file = sys.argv[5] # MetaCyc Monomer -> RXNS
dirpath = os.path.dirname(pf_file)  # the working dir
base = os.path.basename(pf_file)  # the original file
map_file = os.path.join(dirpath, 'orf_map.txt')  # keep an ORF map
bak_file = os.path.join(dirpath, '0.pf.bak')  # backup the orig file


# Build out mapping dictionaries for KO -> EC
orf_df = pd.read_csv(orf_file, sep='\t', header=0)
ko_df = pd.read_csv(ko_file, sep='\t', header=0)
mc_df = pd.read_csv(mc_file, sep='\t', header=0)
rxn_df = pd.read_csv(rxn_file, sep='\t', header=0)
mc_df['MC'] = [x.replace('gnl|META|', '') for x in mc_df['target']]
ec_df = ko_df[['KO', 'EC']]
kegg_df = orf_df[['# ORF_ID', 'KEGG']]
meta_df = mc_df[['#query', 'MC']]
meta_df.columns = ['ORF_ID', 'MC']
orf_ec_df = pd.merge(kegg_df, ec_df, left_on='KEGG', right_on='KO', how='left')
orf_ec_df = orf_ec_df[['# ORF_ID', 'KO', 'EC']]
orf_ec_df.columns = ['ORF_ID', 'KO', 'EC']
orf_mc_df = pd.merge(orf_ec_df, meta_df, on='ORF_ID', how='outer')
orf_rxn_df = pd.merge(orf_mc_df, rxn_df, on='MC', how='outer').drop_duplicates()
orf_map_df = orf_rxn_df.query("EC != ''")

# Iterate over the entries in 0.pf to remove duplicates and collect ORFs
funct_dict = {}
d_list = []
with open(pf_file, 'r') as pf_in:
	data = pf_in.read().split('//')
	for d in data:
		d_split = d.split('\n')
		tmp_dict = {}
		skip_orf = True
		for l in d_split:
			l = l.replace('\n', '')
			cat = l.split('\t', 1)[0]
			if l != '':
				if cat == "ID":
					pf_id = l.split('\t', 1)[1]
					orf_id = pf_id.split('_', 1)[1]
					sub_map_df = orf_map_df.query("ORF_ID == @orf_id")
					tmp_dict["ID"] = pf_id
				if sub_map_df.shape[0] != 0: # ORF needs entries
					if cat == "FUNCTION":
						l = l.replace('MULTISPECIES: ', '')  # remove multispecies tag
						if ' OS ' in l:  # if has taxa info
							split_l = l.split(' OS ', 1)
							l_func = split_l[0].split('\t', 1)[1]
							l_note = split_l[1]
							tmp_dict["FUNCTION"] = l_func
							tmp_dict["GENE-COMMENT"] = l_note
						else:
							l_func = l.split('\t', 1)[1]
							tmp_dict["FUNCTION"] = l_func
						if l_func not in funct_dict.keys():
							skip_orf = False
							funct_dict[l_func] = [pf_id]
							rxn_list = list(sub_map_df['RXN'].dropna().unique())
							ec_list = list(sub_map_df['EC'].dropna().unique())
							ko_list = list(sub_map_df['KO'].dropna().unique())
							if len(rxn_list) != 0:
								tmp_dict["RXN"] = rxn_list
							elif len(ec_list) != 0:
								tmp_dict["EC"] = ec_list
							if len(ko_list) != 0:
								if len(ko_list) > 1:
									flurp
								else:
									tmp_dict["KO"] = ko_list[0]							
						else:
							skip_orf = True
							funct_dict[l_func].append(pf_id)
					elif "PRODUCT-TYPE" in l.split('\t', 1)[0]:
						tmp_dict["PRODUCT-TYPE"] = l.split('\t', 1)[1]
					else:
						attr = l.split('\t', 1)[0]
						val = l.split('\t', 1)[1]
						tmp_dict[attr] = val					
		if skip_orf == False:
			if (('EC' in tmp_dict.keys()) or ('RXN' in tmp_dict.keys())):
				d_list.append("ID\t" + tmp_dict["ID"])
				d_list.append("NAME\t" + tmp_dict["NAME"])
				d_list.append("STARTBASE\t" + tmp_dict["STARTBASE"])
				d_list.append("ENDBASE\t" + tmp_dict["ENDBASE"])
				d_list.append("FUNCTION\t" + tmp_dict["FUNCTION"])
				d_list.append("GENE-COMMENT\t" + tmp_dict["GENE-COMMENT"])
				if 'KO' in tmp_dict.keys():
					d_list.append("DBLINK\tKO:" + tmp_dict["KO"])
				if 'RXN' in tmp_dict.keys():
					for r in tmp_dict["RXN"]:
						d_list.append("METACYC\t" + r)
				elif 'EC' in tmp_dict.keys():
					for e in tmp_dict["EC"]:
						d_list.append("EC\t" + e)
				d_list.append("PRODUCT-TYPE\t" + tmp_dict["PRODUCT-TYPE"])
				d_list.append("//")

# Move the orig file to the backup
os.rename(pf_file, bak_file)

# Save the new deduped data as the orig filename
with open(pf_file, 'w') as pf_out:
	pf_out.write('\n'.join(d_list))

# Save the complete list of ORFs that were deduped
with open(map_file, 'w') as map_out:
	for l in funct_dict.keys():
		map_out.write('\t'.join(funct_dict[l]) + '\n')

#########################################################################
################################## END ##################################
#########################################################################

