import os
import sys
import pandas as pd
from tqdm import tqdm

mc_dir = sys.argv[1] # './tier12-26.5/metacyc/26.5/data/'
out_dir = sys.argv[2]
class_file = os.path.join(mc_dir, "classes.dat")
path_file = os.path.join(mc_dir, "pathways.dat")

print('Loading MetaCyc Classes...')
with open(class_file, 'rb') as class_in:
	data = class_in.read().decode(errors='replace').split('\n//')
class2types = {}
class2common = {}
for frame in data:
	lines = frame.split('\n')
	for line in lines:
		if line:
			if line[0] != '#':
				if 'UNIQUE-ID' in line:
					uniq_id = line.split(' - ')[1]
					class2common[uniq_id] = uniq_id
				elif 'TYPES' in line:
					type_id = line.split(' - ')[1]
					if uniq_id in class2types:
						class2types[uniq_id].append(type_id)
					else:
						class2types[uniq_id] = [type_id]
				elif 'COMMON-NAME' in line:
					common_name = line.split(' - ')[1]
					class2common[uniq_id] = common_name

print('Loading MetaCyc Pathways...')
with open(path_file, 'rb') as path_in:
	data = path_in.read().decode(errors='replace').split('\n//')
path2types = {}
path2common = {}
for frame in data:
	lines = frame.split('\n')
	for line in lines:
		if line:
			if line[0] != '#':
				if 'UNIQUE-ID' in line:
					uniq_id = line.split(' - ')[1]
				elif 'TYPES' in line:
					type_id = line.split(' - ')[1]
					if uniq_id in path2types:
						path2types[uniq_id].append(type_id)
					else:
						path2types[uniq_id] = [type_id]
				elif 'COMMON-NAME' in line:
					common_name = line.split(' - ')[1]
					path2common[uniq_id] = common_name

print('Building MetaCyc Ontology...')
ont_ids = {}
ont_common = {}
stop_list = set()
for k in tqdm(class2types, disable=not sys.stderr.isatty()):
	while k not in stop_list:
		if k in class2common:
			if k in ont_ids:
				mappings = ont_ids[k]
				com_mappings = ont_common[k]
				new_maps = []
				new_com_maps = []
				map_cnt = len(mappings)
				p_cnt = 0
				for i, m in enumerate(mappings):
					com_m = com_mappings[i]
					m0 = m[0]
					parents = class2types[m0]
					for p in parents:
						if p in class2common:
							new_m = [p] + m
							new_maps.append(new_m)
							new_com_m = [class2common[p]] + com_m
							new_com_maps.append(new_com_m)
						else:
							new_maps.append(m)
							new_com_maps.append(com_m)
							p_cnt += 1
					ont_ids[k] = new_maps
					ont_common[k] = new_com_maps
				if ((map_cnt == p_cnt) | ('FRAMES' in parents)):
					stop_list.add(k)
			else:
				parents = class2types[k]
				ont_ids[k] = []
				ont_common[k] = []
				for p in parents:
					if p in class2common:
						ont_ids[k].append([p, k])
						ont_common[k].append([class2common[p],
											   class2common[k]
											   ])
					else:
						stop_list.add(k)
		else:
			stop_list.add(k)

print('Saving MetaCyc Ontology...')
path_ont_list = []
h_prune_list = ['Generalized-Reactions', 'Pathways']
c_prune_list = ['All Pathways and Reactions', 'Pathways']
for path in path2types:
	for t in path2types[path]:
		for i,h in enumerate(ont_ids[t]):
			h = list(h)
			c = list(ont_common[t][i])
			for p in h_prune_list:
				if p in h:
					h.remove(p)
			for p in c_prune_list:
				if p in c:
					c.remove(p)
			common_name = path2common[path]
			path_ont_list.append([path, common_name,
								   '|'.join(h + [path]),
								   '|'.join(c + [common_name])
								   ])
path_ont_df = pd.DataFrame(path_ont_list,
							columns=[
									 'BioCyc_ID',
									 'Common_Name',
									 'MetaCyc_Ontology_IDs',
									 'MetaCyc_Ontology_Names'])
path_ont_df.to_csv(os.path.join(out_dir, 'MetaCyc_PWY_Ontology.tsv'),
					sep='\t', index=False
					)
