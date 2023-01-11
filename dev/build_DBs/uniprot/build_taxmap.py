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
