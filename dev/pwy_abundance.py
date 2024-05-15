#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse

def calculate_metrics(orfs_df, pwys_df):
    # Convert ORFS column to list of ORFs for each pathway
    pwys_df['ORFS'] = pwys_df['ORFS'].apply(lambda x: x.split(',') if pd.notnull(x) and x != '' else [])

    # Initialize new columns
    pwys_df['Average RPKM'] = np.nan
    pwys_df['Weighted Average RPKM'] = np.nan
    pwys_df['Sum TPM'] = np.nan  # Changed from 'Average TPM' to 'Sum TPM'
    pwys_df['Sum Count'] = np.nan
    pwys_df['Relative Abundance'] = np.nan  # New column for relative abundance

    total_count = orfs_df['Count'].sum()  # Calculate total count for relative abundance calculation

    for index, row in pwys_df.iterrows():
        orfs = row['ORFS']
        relevant_orfs = orfs_df[orfs_df['Gene_ID'].isin(orfs)]

        # Calculate metrics
        pwys_df.at[index, 'Average RPKM'] = relevant_orfs['RPKM'].mean()
        pwys_df.at[index, 'Weighted Average RPKM'] = (relevant_orfs['RPKM'] * relevant_orfs['Count']).sum() / relevant_orfs['Count'].sum()
        pwys_df.at[index, 'Sum TPM'] = relevant_orfs['TPM'].sum()  # Summing TPM instead of averaging
        pwys_df.at[index, 'Sum Count'] = relevant_orfs['Count'].sum()
        pwys_df.at[index, 'Relative Abundance'] = (relevant_orfs['Count'].sum() / total_count)  # Calculating relative abundance

    return pwys_df

def main(orfs_file, pwys_file, output_file):
    orfs_df = pd.read_csv(orfs_file, sep='\t')
    pwys_df = pd.read_csv(pwys_file, sep='\t')

    updated_pwys_df = calculate_metrics(orfs_df, pwys_df)
    updated_pwys_df.to_csv(output_file, sep='\t', index=False)
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process ORF and PWY tables to calculate pathway metrics.")
    parser.add_argument("orfs_file", help="File path for the ORF counts table")
    parser.add_argument("pwys_file", help="File path for the PWY table")
    parser.add_argument("output_file", help="File path for the output PWY table with added metrics")

    args = parser.parse_args()

    main(args.orfs_file, args.pwys_file, args.output_file)

