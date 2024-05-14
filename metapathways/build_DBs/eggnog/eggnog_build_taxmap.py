#!/usr/bin/env python
## -*- python -*-
"""Building EggNOG TaxMap for MetaPathways

Usage:
    eggnog_build_taxmap.py --taxid_info <tax> --out_dir <dir>

Options:
    -h --help   Show this screen.
    --version   Show version.
    --taxid_info=FILE   TaxID info file from EggNOG, example [e5.taxid_info.tsv].
    --out_dir=DIR   Path to save all outputs.
"""


import pandas as pd
import sys
import numpy as np
from docopt import docopt
import os


def build_taxmap(eggtax_file):
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
    return final_df



###############################################################
# Run taxmap
arguments = docopt(__doc__, version='EggNOG TaxMap 1.0')

en_tax = arguments['--taxid_info']
out_dir = arguments['--out_dir']

sv_file = os.path.join(out_dir, 'ncbi_taxonomy_tree.txt')

# Build EggNOG Taxonomy
tm_df = build_taxmap(en_tax)

# Save taxmap file
tm_df.to_csv(sv_file, sep='\t', index=False, header=False)

