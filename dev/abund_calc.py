#!/usr/bin/env python3
## -*- python -*-
import argparse
import pandas as pd


def calculate_rpkm(counts, gene_lengths):
    """
    Calculate RPKM (Reads Per Kilobase Million) for each gene.

    :param counts: List of read counts for each gene.
    :param gene_lengths: List of gene lengths in base pairs.
    :return: List of RPKM values for each gene.
    """
    rpkm_values = []
    total_reads = sum(counts)
    for count, length in zip(counts, gene_lengths):
        rpkm = (count / (length / 1000.0)) / (total_reads / 1e6)
        rpkm_values.append(rpkm)
    return rpkm_values


def calculate_tpm(counts, gene_lengths):
    """
    Calculate TPM (Transcripts Per Million) for each gene.

    :param counts: List of read counts for each gene.
    :param gene_lengths: List of gene lengths in base pairs.
    :return: List of TPM values for each gene.
    """
    tpm_values = []
    total_reads = sum(counts)
    counts_per_gene = [count / (length / 1000.0) for count, length in zip(counts, gene_lengths)]
    sum_counts_per_gene = sum(counts_per_gene)
    for count in counts_per_gene:
        tpm = (count / sum_counts_per_gene) * 1e6
        tpm_values.append(tpm)
    return tpm_values


def clean_gene_id(merge_row):
    gen_id = merge_row["Gene_ID"]
    feature = merge_row["feature"]
    if feature == "CDS":
        gene_id = gen_id.split("|")[1]
    elif feature == "rRNA":
        gene_id = gen_id.split("|")[0].split('-')[-1] + '-' + gen_id.split("|")[1]
    elif ((feature == 'tRNA') | (feature == 'pseudogene')):
        gene_id = gen_id.split("|")[1].split('-', 2)[2]
    return gene_id
    

def main():
    parser = argparse.ArgumentParser(description="Calculate RPKM and TPM values for gene expression data.")
    parser.add_argument("--output", type=str, help="Output file name", required=True)
    parser.add_argument("--counts-file", type=str, help="Path to the file containing read counts in TSV format", required=True)
    parser.add_argument("--gene-lengths-file", type=str, help="Path to the file containing gene lengths in TSV format", required=True)
    parser.add_argument("--gtf-file", type=str, help="Path to the GTF file", required=True)

    args = parser.parse_args()

    # Read counts and gene lengths from input files using pandas
    counts_df = pd.read_csv(args.counts_file, sep="\t", index_col=0, comment="#")
    gene_lengths_df = pd.read_csv(args.gene_lengths_file, sep="\t", index_col=0, header=None)
    gene_lengths_df.columns = ['counts']
    # Make sure the input DataFrames have the same index (gene IDs)
    if not counts_df.index.equals(gene_lengths_df.index):
        raise ValueError("Gene IDs in the input files do not match.")

    # Get the read counts from the last column
    counts = counts_df.iloc[:, -1].tolist()
    # Get the gene lengths as lists
    gene_lengths = gene_lengths_df.iloc[:, 0].tolist()

    # Calculate RPKM and TPM values
    rpkm_values = calculate_rpkm(counts, gene_lengths)
    tpm_values = calculate_tpm(counts, gene_lengths)

    # Create a pandas DataFrame with the results
    data = {
        "Gene_ID": counts_df.index,
        "Count": counts,
        "RPKM": rpkm_values,
        "TPM": tpm_values
    }
    results_df = pd.DataFrame(data)

    # Read the GTF file using pandas with specified column names
    gtf_columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes", "gene_id"]
    gtf_df = pd.read_csv(args.gtf_file, sep="\t", comment="#", header=None, names=gtf_columns)

    # Extract gene ID information from the GTF file
    attributes = gtf_df["attributes"].str.split(';', expand=True)
    gtf_df["gene_id"] = attributes[attributes[0].str.contains('gene_id')][0].str.extract(r'gene_id "([^"]+)"')

    # Rename the columns of the GTF DataFrame
    gtf_df.rename(columns={0: "seqname", 1: "source", 2: "feature", 3: "start", 4: "end", 5: "score",
                           6: "strand", 7: "frame", 8: "attributes"}, inplace=True)


    # Merge the results table with the GTF information based on gene ID
    merged_df = pd.merge(results_df, gtf_df, left_on="Gene_ID", right_on="gene_id", how="left")
    merged_df.drop(columns=["gene_id"], inplace=True)

    # Clean up the Gene_ID now that all the merging is done
    merged_df['Gene_ID'] = merged_df.apply(clean_gene_id, axis=1)

    # Save the merged results to a tab-separated file using pandas
    merged_df.to_csv(args.output, sep="\t", index=False)
    print(f"Merged results saved to '{args.output}'")


if __name__ == "__main__":
    main()