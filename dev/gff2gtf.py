#!/usr/bin/env python3
## -*- python -*-
import argparse
import pandas as pd


def extract_gene_id(gff_row):
    for attribute in str(gff_row['attribute']).split(';'):
        if '=' not in attribute:
            continue
        key, value = attribute.strip().split('=', 1)
        if key == 'ID' and value:
            return 'gene_id "' + value + '";'
    raise ValueError(f"GFF feature is missing ID: {gff_row['seqname']}:{gff_row['start']}")


def gff_to_gtf(input_file, output_file, feature_types):
    # Read the GFF file into a pandas DataFrame
    gff_df = pd.read_csv(input_file, sep='\t', comment='#', header=None)

    # Filter out the comments and sequence region lines (if any)
    gff_df = gff_df[~gff_df[0].str.startswith('#')]

    # Filter out the rows with fewer than 9 columns (these are invalid entries)
    gff_df = gff_df[gff_df.apply(lambda row: len(row) >= 9, axis=1)]

    # Rename the columns to match GTF format
    gff_df.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    # Filter rows based on the specified feature types
    gff_df = gff_df[gff_df['feature'].isin(feature_types)]

    # Extract gene_id from the attribute column and create a new column 'gene_id'
    gff_df['gene_id'] = gff_df.apply(extract_gene_id, axis=1)

    if gff_df['gene_id'].duplicated().any():
        raise ValueError('Duplicate gene IDs in GFF; regenerate annotations with corrected RNA IDs before counting.')
    for column in ('start', 'end'):
        values = pd.to_numeric(gff_df[column], errors='raise')
        if values.isna().any() or (values < 1).any() or (values % 1 != 0).any():
            raise ValueError(f'Invalid GFF {column} coordinates')
        gff_df[column] = values.astype('int64')

    # Write the GTF output
    gff_df.to_csv(output_file, sep='\t', index=False, columns=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id'], header=False, quoting=3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GFF to GTF")
    parser.add_argument("input_file", help="Path to the input GFF file")
    parser.add_argument("output_file", help="Path to save the output GTF file")
    parser.add_argument("--feature", required=True, help="Comma-separated list of feature types to be converted (e.g., gene,exon,CDS)")
    args = parser.parse_args()

    # Convert the comma-separated string of feature types to a list
    feature_types_list = args.feature.split(',')

    gff_to_gtf(args.input_file, args.output_file, feature_types_list)