import os
import pandas as pd
import argparse
from difflib import SequenceMatcher
from tqdm import tqdm
import multiprocessing
from functools import partial
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy

def rename_column(df, old_substring, new_substring):
    # Find columns containing the old_substring in the name and rename them to new_substring
    for col in df.columns:
        if old_substring in col:
            df.rename(columns={col: new_substring}, inplace=True)
    return df

def rename_columns(df, rename_rules):
    for old_substring, new_substring in rename_rules.items():
        df = rename_column(df, old_substring, new_substring)
    return df

def group_and_sum(df):
    # Group by "Filename" and "uniref90(product)" and explicitly sum the "ReadCount" and "TPM" columns
    sum_columns = ['ReadCount', 'TPM']
    grouped_df = df.groupby(['Filename', 'uniref90(product)'])[sum_columns].sum().reset_index()

    # Pivot the table to have 'Filename' as columns, 'uniref90(product)' as rows, and 'TPM' as values
    pivoted_df = grouped_df.pivot(index='uniref90(product)', columns='Filename', values='TPM').reset_index()

    # Rename the columns to remove 'Filename' from the column names
    pivoted_df.columns.name = None

    return pivoted_df

def compute_similarity_row(index_values, i):
    str1 = index_values[i]
    row_similarities = []
    for j in range(len(index_values)):
        str2 = index_values[j]
        similarity_ratio = SequenceMatcher(None, str1, str2).ratio()
        row_similarities.append(similarity_ratio)
    return row_similarities

def compute_similarity_matrix_parallel(df, num_threads=1):
    index_values = df['uniref90(product)']
    similarity_matrix = []

    # Use partial to create a picklable function with predefined arguments
    compute_similarity_row_partial = partial(compute_similarity_row, index_values)

    # Use multiprocessing Pool for parallel computation
    with multiprocessing.Pool(processes=num_threads) as pool:
        with tqdm(total=len(index_values), desc='Computing Similarity Matrix') as pbar:
            for row_similarities in pool.imap(compute_similarity_row_partial, range(len(index_values))):
                similarity_matrix.append(row_similarities)
                pbar.update(1)

    return pd.DataFrame(similarity_matrix, columns=index_values, index=index_values)

def cluster_similarity_matrix(similarity_matrix, similarity_threshold=0.97):
    clustered_indices = []
    visited_indices = set()

    for idx, row in similarity_matrix.iterrows():
        if idx not in visited_indices:
            cluster_indices = row[row >= similarity_threshold].index.tolist()
            clustered_indices.append(cluster_indices)
            visited_indices.update(cluster_indices)

    return clustered_indices

def select_representatives(clustered_indices):
    representatives_dict = {}
    for cluster in clustered_indices:
        representative = max(cluster, key=len)
        representatives_dict[representative] = cluster
    return representatives_dict

def sum_tpm_by_representatives(filtered_df, representatives_dict):
    filtered_df.set_index('uniref90(product)', inplace=True)
    rows = []
    for representative, members in representatives_dict.items():
        subset_df = filtered_df[filtered_df.index.isin(members)]
        summed_values = subset_df.sum()
        rows.append([representative] + summed_values.tolist())

    columns = ['uniref90(product)'] + list(filtered_df.columns)
    summed_data = pd.DataFrame(rows, columns=columns)

    return summed_data

def concatenate_and_add_filename(tsv_files):
    combined_data = pd.DataFrame()
    with tqdm(total=len(tsv_files), desc='Concatenating Files') as pbar:
        for file in tsv_files:
            df = pd.read_csv(file, sep='\t')
            filename = os.path.basename(file)
            # Remove ".merge" from the filename
            filename_without_ext = os.path.splitext(filename)[0].replace(".merge", "")
            df['Filename'] = filename_without_ext

            # Define the renaming rules
            rename_rules = {
                'TPM': 'TPM',
                'Read Count': 'ReadCount',
                'RPKM': 'RPKM'
                # Add more substrings and their corresponding replacements here if needed
            }

            # Rename columns based on the defined rules
            df = rename_columns(df, rename_rules)

            combined_data = pd.concat([combined_data, df], ignore_index=True)
            pbar.update(1)

    grouped_and_pivoted_data = group_and_sum(combined_data)

    return grouped_and_pivoted_data

def compute_similarity_matrix(grouped_and_pivoted_data):
    index_values = grouped_and_pivoted_data['uniref90(product)']
    similarity_matrix = []

    with tqdm(total=len(index_values), desc='Computing Similarity Matrix') as pbar:
        for i in range(len(index_values)):
            str1 = index_values[i]
            row_similarities = []
            for j in range(len(index_values)):
                str2 = index_values[j]
                similarity_ratio = SequenceMatcher(None, str1, str2).ratio()
                row_similarities.append(similarity_ratio)
            similarity_matrix.append(row_similarities)
            pbar.update(1)

    return pd.DataFrame(similarity_matrix, columns=index_values, index=index_values)

def filter_by_tpm_and_uniref_threshold(grouped_and_pivoted_data, tpm_threshold=0.001):
    # Define strings to exclude from 'uniref90(product)'
    excluded_strings = ["uncharacterized", "hypothetical",
                        "(fragment)", "unannotated", "unknown",
                        "domain-containing", "ATP", "ABC",
                        "transporter", "non-specific"]
    #Exclude rows with 'uniref90(product)' containing the specified strings
    exclude_df = grouped_and_pivoted_data[~grouped_and_pivoted_data['uniref90(product)'
                                         ].str.lower().str.contains('|'.join(excluded_strings),
                                            case=False)
                                         ]

    # Calculate the total TPM for each column
    total_tpm_per_column = exclude_df.drop(columns=['uniref90(product)']).sum()

    # Filter rows that have at least 0.1% of the total TPM for at least one column
    valid_rows = exclude_df[exclude_df.drop(columns=['uniref90(product)'])
         .apply(lambda row: any(row >= total_tpm_per_column * tpm_threshold), axis=1)
         ]

    valid_rows.reset_index(inplace=True, drop=True)

    return valid_rows

def file_exists(file_path):
    return os.path.exists(file_path)

def file_is_empty(file_path):
    if not os.path.exists(file_path):
        return True
    return os.path.isfile(file_path) and os.stat(file_path).st_size == 0

def save_as_tsv(df, file_path):
    df.to_csv(file_path, sep='\t', index=False)
    print(f"Data saved to {file_path}")

def save_heatmap_to_pdf(summed_data, output_file):
    summed_data.set_index('uniref90(product)', inplace=True)

    # Create a heatmap using seaborn
    plt.figure(figsize=(24, 36))

    # Plot the heatmap
    sns.heatmap(summed_data, cmap='YlGnBu', annot=True, fmt=".1f", linewidths=0.5, cbar_kws={"label": "TPM Sum"})

    plt.title('Summed TPM Heatmap', fontsize=16)  # Set the title
    plt.xlabel('Filename', fontsize=12)  # Set the x-axis label
    plt.ylabel('uniref90(product)', fontsize=12)  # Set the y-axis label
    plt.tight_layout()

    # Save the heatmap to PDF
    plt.savefig(output_file, format='pdf', bbox_inches="tight")
    plt.close()
    
    print(f"Heatmap saved to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenate TSV files, add a filename column, and group the data.')
    parser.add_argument('files', nargs='+', help='List of TSV files to concatenate')
    parser.add_argument('--pivot-output', help='Output file path for the pivot table in TSV format')
    parser.add_argument('--similarity-output', help='Output file path for the similarity matrix in TSV format')
    parser.add_argument('--representatives-output', help='Output file path for representatives in TSV format')
    parser.add_argument('--summed-tpm-output', help='Output file path for summed TPM data in TSV format')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for multiprocessing')
    parser.add_argument('--force', action='store_true', help='Force recomputation and overwrite existing files')
    args = parser.parse_args()

    tsv_files = args.files
    pivot_output_file = args.pivot_output
    similarity_output_file = args.similarity_output
    representatives_output_file = args.representatives_output
    summed_tpm_output_file = args.summed_tpm_output
    num_threads = args.threads
    force_recomputation = args.force

    try:
        grouped_and_pivoted_data = None
        similarity_matrix = None
        representatives_dict = None
        summed_data = None

        # Check if pivot table file exists and load if available
        if not force_recomputation and not file_is_empty(pivot_output_file):
            grouped_and_pivoted_data = pd.read_csv(pivot_output_file, sep='\t')
            print(f"Using existing pivot table from {pivot_output_file}")
        else:
            grouped_and_pivoted_data = concatenate_and_add_filename(tsv_files)
            save_as_tsv(grouped_and_pivoted_data, pivot_output_file)

        print("Grouped and Pivoted Data:")
        print(grouped_and_pivoted_data)

        # Filter rows based on TPM threshold
        filtered_data = filter_by_tpm_and_uniref_threshold(grouped_and_pivoted_data)
        print("Filtered Data:")
        print(filtered_data)

        # Check if similarity matrix file exists and load if available
        if not force_recomputation and not file_is_empty(similarity_output_file):
            similarity_matrix = pd.read_csv(similarity_output_file, sep='\t', index_col=0)
            print(f"Using existing similarity matrix from {similarity_output_file}")
        else:
            if num_threads > 1:
                similarity_matrix = compute_similarity_matrix_parallel(filtered_data, num_threads)
            else:
                similarity_matrix = compute_similarity_matrix(filtered_data)
            save_as_tsv(similarity_matrix, similarity_output_file)

        print("Similarity Matrix:")
        print(similarity_matrix)

        # Check if representatives file exists and load if available
        if not force_recomputation and not file_is_empty(representatives_output_file):
            with open(representatives_output_file, 'r') as f:
                representatives_dict = {}
                for line in f:
                    representative, *members = line.strip().split('\t')
                    representatives_dict[representative] = members
            print(f"Using existing representatives from {representatives_output_file}")
        else:
            # Cluster the similarity matrix and select representatives
            clustered_indices = cluster_similarity_matrix(similarity_matrix)
            representatives_dict = select_representatives(clustered_indices)
            with open(representatives_output_file, 'w') as f:
                for representative, members in representatives_dict.items():
                    members_str = '\t'.join(members)
                    f.write(f"{representative}\t{members_str}\n")
            print(f"Representatives saved to {representatives_output_file}")

        print("Representatives selectied")
        #print(representatives_dict)

        # Check if summed TPM data file exists and load if available
        if not force_recomputation and not file_is_empty(summed_tpm_output_file):
            summed_data = pd.read_csv(summed_tpm_output_file, sep='\t')
            print(f"Using existing summed TPM data from {summed_tpm_output_file}")
        else:
            # Sum the TPM values for each representative
            summed_data = sum_tpm_by_representatives(filtered_data, representatives_dict)
            save_as_tsv(summed_data, summed_tpm_output_file)

        print("Summed TPM Data:")
        print(summed_data)

        # Save the summed TPM data as a clustered heatmap to PDF
        save_heatmap_to_pdf(summed_data, 'summed_tpm_heatmap.pdf')

    except Exception as e:
        print("An error occurred:", str(e))
