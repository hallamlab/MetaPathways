# This script takes a '0.pf' file as input
# The 0.pf is automatically output by MP within the ./ptools dir
# It is meant to deduplicate the ORFs based on exact match to FUNCTION
import sys
import os
import pandas as pd


pf_file = sys.argv[1]  # the 0.pf file
orf_file = sys.argv[2] # ORF annotation file from MP3 results dir
ko_file = sys.argv[3] # Mapping file from Uniprot -> KO -> EC
dirpath = os.path.dirname(pf_file)  # the working dir
base = os.path.basename(pf_file)  # the original file
map_file = os.path.join(dirpath, 'orf_map.txt')  # keep an ORF map
bak_file = os.path.join(dirpath, '0.pf.bak')  # backup the orig file


# Build out mapping dictionaries for KO -> EC
orf_df = pd.read_csv(orf_file, sep='\t', header=0)
ko_df = pd.read_csv(ko_file, sep='\t', header=0)
ec_df = ko_df[['KO', 'EC']]
kegg_df = orf_df[['# ORF_ID', 'KEGG']]
orf_ec_df = pd.merge(kegg_df, ec_df, left_on='KEGG', right_on='KO', how='left').drop_duplicates().fillna('')
orf_map_df = orf_ec_df[['# ORF_ID', 'KO', 'EC']]
orf_map_df.columns = ['ORF_ID', 'KO', 'EC']
orf2ko_dict = {'O_' + k: v for k,v in zip(orf_map_df['ORF_ID'], orf_map_df['KO']) if v != ''}
orf2ec_dict = {'O_' + k: v for k,v in zip(orf_map_df['ORF_ID'], orf_map_df['EC']) if v != ''}


# Iterate over the entries in 0.pf to remove duplicates and collect ORFs
ec_dict = {}
d_list = []
with open(pf_file, 'r') as pf_in:
	data = pf_in.read().split('//')
	for d in data:
		d_split = d.split('\n')
		d_tmp = []
		skip_orf = True
		for l in d_split:
			l = l.replace('\n', '')
			cat = l.split('\t', 1)[0]
			if l != '':
				if cat == "ID":
					orf_id = l.split('\t', 1)[1]
				if orf_id in orf2ec_dict.keys():  # Only write entries that have ECs
					if cat == "FUNCTION":
						l = l.replace('MULTISPECIES: ', '')  # remove multispecies tag
						if ' OS ' in l:  # if has taxa info
							split_l = l.split(' OS ', 1)
							l_func = split_l[0]
							l_note = 'GENE-COMMENT\t' + split_l[1]
							d_tmp.append(l_func)
							d_tmp.append(l_note)
						else:
							d_tmp.append(l)
						if orf_id in orf2ec_dict.keys():
							ec = orf2ec_dict[orf_id]
							if ':' in ec:
								ec = ec.split(':')[1]
							ec_l = "EC\t" + ec
							if ec not in ec_dict.keys():  # Don't add ORFs with repeat ECs
								ec_dict[ec] = [orf_id]
								d_tmp.append(ec_l)
								if orf_id in orf2ko_dict.keys():
									ko = orf2ko_dict[orf_id]
									ko_link_l = "DBLINK\tKO:" + ko
									d_tmp.append(ko_link_l)
								skip_orf = False
							else:  # If a deplicate EC, add ORF to mapping list
								ec_dict[ec].append(orf_id)
								skip_orf = True
					elif "PRODUCT-TYPE" in l.split('\t', 1)[0]:
						d_tmp.append(l)
						d_tmp.append('//')
					elif cat != "DBLINK":
						d_tmp.append(l)
					else:
						if "PRODUCT-TYPE" in l.split('\t', 1)[0]:
							d_tmp.append(l)
							d_tmp.append('//')
						else:
							d_tmp.append(l)
		if skip_orf == False:
			d_list.extend(d_tmp)


# Move the orig file to the backup
os.rename(pf_file, bak_file)

# Save the new deduped data as the orig filename
with open(pf_file, 'w') as pf_out:
	pf_out.write('\n'.join(d_list))

# Save the complete list of ORFs that were deduped
with open(map_file, 'w') as map_out:
	for l in ec_dict.keys():
		map_out.write('\t'.join(ec_dict[l]) + '\n')

#########################################################################
################################## END ##################################
#########################################################################

