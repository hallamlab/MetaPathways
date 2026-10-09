#!/usr/bin/env python3
## -*- python -*-
import argparse
import pandas as pd
import numpy as np


def calculate_rpkm(counts, gene_lengths):
    """
    Calculate RPKM (Reads Per Kilobase Million) for each gene.

    :param counts: List of read counts for each gene.
    :param gene_lengths: List of gene lengths in base pairs.
    :return: List of RPKM values for each gene.
    """
    
    counts = np.array(counts)
    gene_lengths = np.array(gene_lengths)
    total_reads = sum(counts)
    tot_per_million = total_reads / 1e6
    rpkm_values = counts / (gene_lengths/1000 * tot_per_million)
    
    return rpkm_values


def calculate_tpm(counts, gene_lengths):
    """
    Calculate TPM (Transcripts Per Million) for each gene.

    :param counts: List of read counts for each gene.
    :param gene_lengths: List of gene lengths in base pairs.
    :return: List of TPM values for each gene.
    """

    counts = np.array(counts)
    gene_lengths = np.array(gene_lengths)
    tpm_values = []
    counts_per_gene = counts / (gene_lengths / 1000.0)
    counts_per_million = sum(counts_per_gene) / 1e6
    tpm_values = counts_per_gene / counts_per_million

    return tpm_values


def abundance_table(counts_file, gtf_file, gene_lengths_file=None):
    """One abundance row per featureCounts gene ID, using its union feature length."""
    counts_df = pd.read_csv(counts_file, sep="\t", index_col=0, comment="#")
    if counts_df.index.has_duplicates:
        raise ValueError("featureCounts contains duplicate gene IDs")
    if 'Length' in counts_df:
        lengths = pd.to_numeric(counts_df['Length'], errors='raise')
    elif gene_lengths_file:
        legacy = pd.read_csv(gene_lengths_file, sep="\t", index_col=0, header=None)
        if legacy.index.has_duplicates or set(legacy.index) != set(counts_df.index):
            raise ValueError("Legacy gene lengths must contain each counted gene ID exactly once")
        lengths = pd.to_numeric(legacy.iloc[:, 0].reindex(counts_df.index), errors='raise')
    else:
        raise ValueError("Missing featureCounts Length column")
    counts = pd.to_numeric(counts_df.iloc[:, -1], errors='raise')
    if not np.isfinite(lengths).all() or (lengths <= 0).any():
        raise ValueError("Gene lengths must be finite and positive")
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    results = pd.DataFrame({
        'Gene_ID': counts_df.index, 'Count': counts.to_numpy(),
        'Length': lengths.to_numpy(),
        'RPKM': calculate_rpkm(counts, lengths) if counts.sum() else np.zeros(len(counts)),
        'TPM': calculate_tpm(counts, lengths) if counts.sum() else np.zeros(len(counts)),
    })
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None, names=columns,
                      dtype=str, keep_default_na=False)
    gtf['Gene_ID'] = gtf['attributes'].str.extract(r'(?:^|;)\s*gene_id "([^"]+)"')
    if gtf['Gene_ID'].isna().any():
        raise ValueError("GTF row is missing gene_id")
    missing = set(counts_df.index) - set(gtf['Gene_ID'])
    uncounted = set(gtf['Gene_ID']) - set(counts_df.index)
    if missing or uncounted:
        raise ValueError(f'Gene IDs differ between counts and GTF: missing from GTF={sorted(missing)[:5]}, missing from counts={sorted(uncounted)[:5]}')
    if gtf['Gene_ID'].duplicated().any():
        duplicates = gtf.loc[gtf['Gene_ID'].duplicated(), 'Gene_ID'].unique()[:5]
        raise ValueError(f"GTF contains duplicate gene IDs: {list(duplicates)}. Regenerate annotations with corrected RNA IDs; do not pool distinct loci.")
    return results.merge(gtf, on='Gene_ID', how='left', validate='one_to_one')



def main():
    parser = argparse.ArgumentParser(description="Calculate RPKM and TPM values for gene expression data.")
    parser.add_argument("--output", type=str, help="Output file name", required=True)
    parser.add_argument("--counts-file", type=str, help="Path to the file containing read counts in TSV format", required=True)
    parser.add_argument("--gene-lengths-file", type=str, help="Legacy length table (optional); featureCounts Length is preferred", required=False)
    parser.add_argument("--gtf-file", type=str, help="Path to the GTF file", required=True)

    args = parser.parse_args()

    merged_df = abundance_table(args.counts_file, args.gtf_file, args.gene_lengths_file)

    # Save the merged results to a tab-separated file using pandas
    merged_df.to_csv(args.output, sep="\t", index=False)
    print(f"Merged results saved to '{args.output}'")


if __name__ == "__main__":
    main()