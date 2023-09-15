import urllib.request
import time
import pandas as pd
import sys
from tqdm import tqdm
import pyfastx
import numpy as np
import os


def build_prot_dict(dat):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    prot_dict = {}
    for frame in data:
        lines = frame.split('\n')
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        enzrxn_list = []
                        compof_list = []
                        uniq_id = line.split(' - ')[1]
                    if 'CATALYZES - ' in line:
                        enzrxn_id = line.split(' - ')[1]
                        enzrxn_list.append(enzrxn_id)
                    if 'COMPONENT-OF - ' in line:
                        compof_id = line.split(' - ')[1]
                        compof_list.append(compof_id)
        if uniq_id not in prot_dict:
            prot_dict[uniq_id] = [enzrxn_list, compof_list]

    return prot_dict


def build_enz_dict(dat):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    enzrxn_dict = {}
    for frame in data:
        lines = frame.split('\n')
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        rxn_list = []
                        uniq_id = line.split(' - ')[1]
                    if 'REACTION - ' in line:
                        rxn_id = line.split(' - ')[1]
                        rxn_list.append(rxn_id)
        if uniq_id not in enzrxn_dict:
            enzrxn_dict[uniq_id] = rxn_list
        if len(rxn_list) > 1:
            print(len(rxn_list))

    return enzrxn_dict


def get_enz_comp(p_dict, k):

    enz_list = p_dict[k][0]
    com_list = p_dict[k][1]

    return enz_list, com_list


def build_cmp_dict(dat):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    cmp_dict = {}
    for frame in data:
        lines = frame.split('\n')
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        cmp_name = ''
                        uniq_id = line.split(' - ')[1]
                    if 'COMMON-NAME - ' in line:
                        cmp_name = line.split(' - ')[1]
        if uniq_id not in cmp_dict:
            cmp_dict[uniq_id] = cmp_name

    return cmp_dict


def build_rxn_dict(dat):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    rxn_dict = {}
    rxndir_dict = {}
    for frame in data:
        lines = frame.split('\n')
        rxn_dir = ''
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        rxn_name = ''
                        uniq_id = line.split(' - ')[1]
                    if 'COMMON-NAME - ' in line:
                        rxn_name = line.split(' - ')[1]
                    if 'REACTION-DIRECTION - ' in line:
                        rxn_dir = line.split(' - ')[1]
        if uniq_id not in rxn_dict:
            rxn_dict[uniq_id] = rxn_name
            rxndir_dict[uniq_id] = rxn_dir

    return rxn_dict, rxndir_dict


def build_pwycmp_dict(dat, cmp_d):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    pwy_dict = {}
    pwy2rxn_dict = {}
    for frame in data:
        lines = frame.split('\n')
        skip_super = False
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        rxn_dict = {}
                        uniq_id = line.split(' - ')[1]
                        #print(uniq_id)
                    if 'COMMON-NAME - ' in line:
                        pwy_name = line.split(' - ')[1]
                    if 'TYPES - ' in line:
                        if 'Super-Pathways' in line:
                            skip_super = True
                    if (('REACTION-LAYOUT - ' in line) & (skip_super == False)):
                        #print(line)
                        cmp_list = []
                        rxn_desc = line.split(' - ')[1]
                        rxn = rxn_desc.strip('(').split('(:')[0].strip(' ')
                        l_prime = [cmp_d[c] if c in cmp_d else c for c in line.split('(:LEFT-PRIMARIES')[1].split(')')[0].strip().split(' ')]
                        direct = line.split('(:DIRECTION')[1].split(')')[0].strip().strip(':')
                        r_prime = [cmp_d[c] if c in cmp_d else c for c in line.split('(:RIGHT-PRIMARIES')[1].split(')')[0].strip().split(' ')]
                        rxn_dict[rxn] = '|'.join(l_prime), direct, '|'.join(r_prime)
        pwy_dict[uniq_id] = pwy_name
        if uniq_id not in pwy2rxn_dict:
            pwy2rxn_dict[uniq_id] = rxn_dict

    return pwy_dict, pwy2rxn_dict


def build_rxnenz_dict(dat):
    with open(dat, 'rb') as path_in:
        data = path_in.read().decode(errors='replace').split('\n//')
    rxn2enz_dict = {}
    enz_dict = {}
    for frame in data:
        lines = frame.split('\n')
        for line in lines:
            if line:
                if line[0] != '#':
                    if 'UNIQUE-ID - ' in line:
                        rxn_list = []
                        uniq_id = line.split(' - ')[1]
                    if 'COMMON-NAME - ' in line:
                        enz_name = line.split(' - ')[1]                        
                    if 'REACTION - ' in line:
                        rxn_id = line.split(' - ')[1]
                        rxn_list.append(rxn_id)
        
        enz_dict[uniq_id] = enz_name
        for rxn in rxn_list:
            if rxn in rxn2enz_dict:
                if uniq_id not in rxn2enz_dict[rxn]:
                    rxn2enz_dict[rxn].append(uniq_id)
            else:
                rxn2enz_dict[rxn] = [uniq_id]
    return enz_dict, rxn2enz_dict



# MetaCyc
mc_dir = sys.argv[1] # './tier12-26.5/metacyc/26.5/data/'
pwy_dat = os.path.join(mc_dir, 'pathways.dat')
cmp_dat = os.path.join(mc_dir, 'compounds.dat')
rxn_dat = os.path.join(mc_dir, 'reactions.dat')
enzrxns_file = os.path.join(mc_dir, 'enzrxns.dat')
prot_fsa = os.path.join(mc_dir, 'protseq.fsa')
prot_file = os.path.join(mc_dir, 'proteins.dat')

# Collect all Monomer IDs
mono_list = []
with open(prot_fsa, 'r') as pfsa:
    data = pfsa.readlines()
    for line in data:
        if line[0] == '>':
            mono_id = line.split('|')[2].split(' ')[0]
            mono_list.append(mono_id)

# Build Protein to Catalyzes dictionary
prot2cat_dict = build_prot_dict(prot_file)
prot2enz_dict = {}
max_len = 0
for m in mono_list:
    prot2enz_dict[m] = []
    enz, com = get_enz_comp(prot2cat_dict, m)
    if enz: # Add enzrxns IDs to dict
        for e in enz:
            prot2enz_dict[m].append(e)
    if com: # Monomer is part of a complex, this can be a nested doll situation
        com_tmp = com
        com_flatlist_tmp = []
        stop = 0
        last_len = len(com_flatlist_tmp)
        while stop == 0:
            cc_tmp = []
            for c in com_tmp:
                if ((c in prot2cat_dict) &
                    (c not in com_flatlist_tmp)
                    ):
                    com_flatlist_tmp.append(c)
                    c_enz, c_com = get_enz_comp(prot2cat_dict, c)
                    for cc in c_com:
                        cc_tmp.append(cc)
            com_tmp.extend(cc_tmp)
            com_flatlist_tmp = list(set(com_flatlist_tmp))
            this_len = len(com_flatlist_tmp)
            if this_len == last_len:
                stop = 1
            else:
                last_len = this_len

        c_enz_list = []
        for c_tmp in com_flatlist_tmp:
            e = prot2cat_dict[c][0]
            c_enz_list.extend(e)
        prot2enz_dict[m].extend(c_enz_list)



# Build Enzyme to Reaction dictionary
enz2rxn_dict = build_enz_dict(enzrxns_file)
enzrxn_list = []
for m in prot2enz_dict:
    enzs = prot2enz_dict[m]
    for e in enzs:
        rxns = enz2rxn_dict[e]
        for r in rxns:
            enzrxn_list.append([m, r])

enzrxns_df = pd.DataFrame(enzrxn_list, columns=['MC', 'RXN'])
enzrxns_df.to_csv(os.path.join(mc_dir, 'MetaCyc-monomer-rxn-pairs.tsv'),
                  sep='\t', index=False
                  )

###### Pathways -> Compounds ######

# Build compounds dictionary
cmp_dict = build_cmp_dict(cmp_dat)

# Build reactions dictionary
rxn_dict, rxndir_dict = build_rxn_dict(rxn_dat)

# Build Enzyme dictionary and Enzyme to Reaction dictionary
enz_dict, rxn2enz_dict = build_rxnenz_dict(enzrxns_file)

# Build Pathway to Compounds dictionary
pwy_dict, pwy2cmp_dict = build_pwycmp_dict(pwy_dat, cmp_dict)

print(pwy_dict['METH-ACETATE-PWY'])
print(pwy2cmp_dict['METH-ACETATE-PWY'])

# Create Pathway -> Reactions -> Compounds
pwy_rxn_cmp_list = []
for p in pwy2cmp_dict:
    for rx in pwy2cmp_dict[p]:
        l, d, r = pwy2cmp_dict[p][rx]
        if rx in rxn2enz_dict:
            for e in rxn2enz_dict[rx]:
                pwy_rxn_cmp_list.append([p, pwy_dict[p], rx, rxn_dict[rx],
                    e, enz_dict[e], l, r, d, rxndir_dict[rx]
                    ])
        else:
            e = ''
            pwy_rxn_cmp_list.append([p, pwy_dict[p], rx, rxn_dict[rx],
                    e, '', l, r, d, rxndir_dict[rx]
                    ])
pwy_rxn_cmp_df = pd.DataFrame(pwy_rxn_cmp_list,
    columns=['PWY_ID', 'PWY_NAME', 'RXN', 'RXN_NAME', 'ENZ_ID', 'ENZ_NAME',
             'left_primary', 'right_primary', 'direction', 'direction_note'
             ])
pwy_rxn_cmp_df.to_csv(os.path.join(mc_dir, 'MetaCyc-PWY-RXN-CMP-map.tsv'),
                      sep='\t', index=False
                      )

