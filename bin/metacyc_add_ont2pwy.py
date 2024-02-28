import pandas as pd 
import sys
import argparse

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

def main(metacyc_pwy_hier, mp_pwy_file, mp_hier_file):
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
    print(f"Output saved to {mp_hier_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ORF and PWY tables to calculate pathway metrics.")
    parser.add_argument("metacyc_pwy_hier", help="Path to MetaCyc Ontology")
    parser.add_argument("mp_pwy_file", help="File path for the PWY table")
    parser.add_argument("mp_hier_file", help="File path for the output PWY table with added Ontology")

    args = parser.parse_args()

    main(args.metacyc_pwy_hier, args.mp_pwy_file, args.mp_hier_file)