import pandas as pd



uniref_file = "UniRef90_idmapping.dat"
kegg2ec_file = "Uniprot_KEGG_EC_map.tsv"


kegg2ec_df = pd.read_csv(kegg2ec_file, sep='\t', header=0)
uniref_df = pd.read_csv(uniref_file, sep='\t', header=None)
uniref_df.columns = ['UniProtKB', 'dattype', 'UniRefID']


merge_df = kegg2ec_df.merge(uniref_df, how='left', on='UniProtKB')

print(kegg2ec_df)
print(uniref_df)
print(merge_df)

merge_df.to_csv("Uniprot_KEGG_EC_map_UniRef90.tsv", sep='\t', index=False)