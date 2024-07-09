import pandas as pd
import re
import sys
import os


expasy_dat = sys.argv[1] #"enzyme.dat"
fc_dir = sys.argv[2]
xref2ec_file = os.path.join(fc_dir, "EC_map.swissprot_test.tsv")

with open(expasy_dat, 'r') as ex_handle:
    ex_dat = ex_handle.read().split('//')
xref2ec_list = []
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
                        xref2ec_list.append([acc, ec_id])
xref2ec_df = pd.DataFrame(xref2ec_list, columns=['AccID', 'EC'])
xref2ec_df.to_csv(xref2ec_file, sep='\t', index=False)