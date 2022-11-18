import pandas as pd

'''
# NCBI Taxonomy
ranklin_dmp = 'new_taxdump/rankedlineage.dmp'
ranklin_df = pd.read_csv(ranklin_dmp, sep='\t\|\t', header=None,
    engine='python', usecols=[0, 1]
    )
ranklin_df.columns = ['taxid', 'name']
print(ranklin_df.head())
taxlin_dmp = 'new_taxdump/taxidlineage.dmp'
taxlin_df = pd.read_csv(taxlin_dmp, sep='\t\|\t', header=None, 
    engine='python'
    )
taxlin_df.columns = ['taxid', 'lineage']
taxlin_df['lineage'] = [x.replace('\t|', '').strip() for x in taxlin_df['lineage']] 
taxlin_df['parent'] = [x.split(' ')[-1] if x.split(' ')[-1] != '' else 1 for x in taxlin_df['lineage']]
print(taxlin_df.head())
merge_df = ranklin_df.merge(taxlin_df, on='taxid', how='left')
print(merge_df.head())
final_df = merge_df[['name', 'taxid', 'parent']]
print(final_df.head())
final_df.to_csv('ncbi_taxonomy_tree.txt', sep='\t', index=False, header=False)
'''
# EggNOG Taxonomy
eggtax_file = 'e5.taxid_info.tsv'
eggtax_df = pd.read_csv(eggtax_file, sep='\t', header=0) 
taxid2name_dict = {}
parent_dict = {}
for i, line in eggtax_df.iterrows():
    names_lin = str(line['Named Lineage']).split(',')
    taxid_lin = str(line['Taxid Lineage']).split(',')
    for t_n in zip(taxid_lin, names_lin):
        t = t_n[0]
        n = t_n[1]
        taxid2name_dict[t] = n
    rev_taxid_lin = taxid_lin[::-1]
    for i, tr in enumerate(rev_taxid_lin):
        if tr != '1':
            parent_dict[tr] = rev_taxid_lin[i + 1]

final_list = []
for k in parent_dict:
    v = parent_dict[k]
    n = taxid2name_dict[k]
    final_list.append([n, k, v])
final_df = pd.DataFrame(final_list, columns=['name', 'taxid', 'parent'])
final_df.to_csv('ncbi_taxonomy_tree.txt', sep='\t', index=False, header=False)
