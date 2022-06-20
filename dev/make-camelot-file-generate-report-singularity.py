#!/usr/bin/env python
## -*- python -*-

"""Usage:
  make-camelot-file-generate-report.py <pgdb_flat_file_dir_path> <camelot_file_output_dir_path> <report_path>
"""

"""
Example command-line usage:

> python3 ./make-camelot-file-generate-report.py ~/data/biocyc-24.0/ecoli/24.0/data /tmp /tmp/report.tsv
"""

from docopt import docopt

from camelot_frs.camelot_frs import get_kb, get_frame, get_frame_all_children, frame_parent_of_frame_p, frame_object_p
from camelot_frs.pgdb_loader import load_pgdb, make_camelot_file
from camelot_frs.pgdb_api    import genes_of_pathway
import os
import shutil
import glob


def get_pwy_inf(reports_dir):
    """
    Accepts the path to the 'reports' directory within
    Ptools flatfile output.

    Returns dictionary of all values found in
    'pwy-inference-report_YYYY-MM-DD.txt' file.
    """
    pwy_inf_rec_list = []
    pwy_inf_file = glob.glob(os.path.join(reports_dir, 'pwy-inference-report_*.txt'))[0]


    with open(pwy_inf_file, 'r') as pwy_inf_in:

        data = pwy_inf_in.read()
        trim_dat = data.split('::: Pathway Inference Report')
        if len(trim_dat) == 3:
            keep_dat = trim_dat[2]
        else:
            keep_dat = trim_dat[1]
        keep_dat = keep_dat.split('List of pathways pruned')[0]

        pwy_inf_rec = ''
        start = False
        for line in keep_dat.split('\n'):
            if line[:2] == ' (': # start of record
                if pwy_inf_rec != '': # add if there is something to add
                    pwy_inf_rec_list.append(pwy_inf_rec)
                    pwy_inf_rec = line
                else: # start a new record
                    pwy_inf_rec = line
                start = True
            elif start == True:
                pwy_inf_rec = pwy_inf_rec + line
        pwy_inf_rec_list.append(pwy_inf_rec) # add last record
    pwy_inf_dict = {}
    for p_rec in pwy_inf_rec_list:
        split_line = p_rec.strip(' (').strip(')').split(' ')
        pwy_k = split_line[0]
        pwy_inf_dict[pwy_k] = {'CONFIDENCE': split_line[2]}

    return pwy_inf_dict




arguments = docopt(__doc__, version='run-reactionary 0.4')


## version.dat file is not in expected directory, create it
flatpath = arguments['<pgdb_flat_file_dir_path>']
verfile= os.path.join(flatpath.rsplit('/', 2)[0], 'default-version')
new_verfile = os.path.join(flatpath, 'version.dat')
shutil.copyfile(verfile, new_verfile)

## Need to create a custom sample_id since org_id is blank
sample_id = os.path.basename(flatpath.rsplit('/', 4)[0].strip('cyc'))
if sample_id == '': # just in case there is an extra '/' at the end
    sample_id = os.path.basename(
                    os.path.dirname(
                        flatpath.rsplit('/', 4)[0].strip('cyc')
                        )
                    )

## Create the .camelot file:
org_id = make_camelot_file(arguments['<pgdb_flat_file_dir_path>'],
                           arguments['<camelot_file_output_dir_path>'])
print('SAMPLE_ID:', sample_id)
print('ORGANISM_ID:', org_id)

## Load the PGDB:
load_pgdb(arguments['<camelot_file_output_dir_path>'] + '/' + org_id + '.camelot')

curr_kb = get_kb(org_id)

# Build Pathway Inference Data Dictionary from contents of ./reports/ dirextory
reportspath = os.path.dirname(flatpath) + '/reports'
pwy_inf_data = get_pwy_inf(reportspath)

## Generate the report:
headers = [ "SAMPLE",
            "PWY_NAME", 	
            "PWY_COMMON_NAME",
            "PWY_SCORE",
            "PWY_CONFIDENCE",
            "NUM_REACTIONS",
            "NUM_COVERED_REACTIONS",
    	    "ORF_COUNT",
            "ORFS" 
           ]

with open(arguments['<report_path>'],"w") as report_fp:

    print('\t'.join(headers),
          file=report_fp)

    for pwy in get_frame_all_children(get_frame(curr_kb, 'Pathways'), frame_types='instance'):
        if not frame_parent_of_frame_p(get_frame(curr_kb, 'Super-Pathways'),
                                       pwy):
            covered_rxn_count = 0
            for rxn in pwy.get_slot_values('REACTION-LIST'):
                if 'ENZYMATIC-REACTION' in rxn.slots:
                    covered_rxn_count += 1
            
            pscore = pwy.get_slot_values('SCORE')[0].strip()
            pconf = pwy_inf_data[pwy.frame_id]['CONFIDENCE']
            print(pwy, pscore, pconf)

            try:
                pwy_gene_names = [ str(gene.get_slot_values('COMMON-NAME')[0]).lstrip('frame:') for gene in genes_of_pathway(pwy) ]
            except Exception:
                pwy_genes_names = ['EcoCyc','error']
            print('\t'.join([sample_id, #  curr_kb.kb_name,
                             pwy.frame_id,
                             pwy.get_slot_values('COMMON-NAME')[0],
                             pscore,
                             pconf,
                             str(len(pwy.get_slot_values('REACTION-LIST'))),
                             str(covered_rxn_count),
                             str(len(pwy_gene_names)),
                             ','.join(pwy_gene_names)]),
                  file = report_fp)
