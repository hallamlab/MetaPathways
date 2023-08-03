import pandas as pd
import re


expasy_dat = "enzyme.dat"
uniref_file = "UniRef90_idmapping.dat"
kegg2ec_file = "Uniprot_KEGG_EC_map.tsv"

with open(expasy_dat, 'r') as ex_handle:
    ex_dat = ex_handle.read().split('//')
xref2ec_dict = {}
for rec in ex_dat:
    if 'ID ' in rec:
        ec_id = rec.split('ID')[1].split('\n')[0].strip(' ')
        rec_lines = rec.split('\n')
        xref_list = []
        for line in rec_lines:
            if line[:2] == "DR":
                sp_line = re.split('; |, | ', line.strip('DR').strip(' ').strip(';'))
                for acc in sp_line:
                    if acc:
                        xref_list.append(acc)
                        if acc in xref2ec_dict:
                            xref2ec_dict[acc].append(ec_id)
                        else:
                            xref2ec_dict[acc] = [ec_id]
print(xref2ec_dict)

ur90_dat = pd.read_csv(uniref_file, sep='\t', names=['UniProtKB', 'dattype', 'UniRefID'])
ur90_dict = pd.Series(ur90_dat['UniRefID'].values, index=ur90_dat['UniProtKB']).to_dict()
ur902ec_list = []
for k in ur90_dict:
    v = ur90_dict[k]
    if k in xref2ec_dict:
        ecs = xref2ec_dict[k]
        for ec in ecs:
            ur902ec_list.append([k, ec, v])

ur902ec_df = pd.DataFrame(ur902ec_list, columns=['UniProtKB', 'EC', 'UniRefID'])
ur902ec_df.to_csv("Uniprot_KEGG_EC_map.tsv", sep='\t', index=False)