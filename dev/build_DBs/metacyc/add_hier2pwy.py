import pandas as pd 
import sys
"""
Adds the MetaCyc Hierarchy (ontology) to the parsed pathways output table
	Table is found in `[sample]/results/pgdb/[commmunity|MAGs]/[sample]_pwy.tsv`
	Relies on prebuilt MetaCyc Hierarchy table that is provided with MetaCyc DB build.

Inputs:
metacyc_pwy_hier: Path to .
mp_pwy_file: List of gene lengths in base pairs.

Output:
mp_hier_file: List of RPKM values for each gene.
"""

metacyc_pwy_hier = sys.argv[1]
mp_pwy_file = sys.argv[2]
mp_hier_file = sys.argv[3]


path_hier_df = pd.read_csv(metacyc_pwy_hier, sep='\t',
						header=0
						)
pwy_table_df = pd.read_csv(mp_pwy_file, sep='\t',
						header=0
						)

merge_df = pwy_table_df.merge(path_hier_df, left_on='PWY_NAME',
							  right_on='BioCyc_ID', how='left'
							  )
trim_df = merge_df[['SAMPLE', 'PWY_NAME', 'PWY_COMMON_NAME',
					'PWY_SCORE', 'PWY_CONFIDENCE', 'NUM_REACTIONS',
					'NUM_COVERED_REACTIONS', 'ORF_COUNT',
					'ORFS', 'MetaCyc_hierarchy'
					]]
trim_df.to_csv(mp_hier_file, sep='\t',
					index=False
					)