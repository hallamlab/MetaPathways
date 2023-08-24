import pandas as pd
import argparse

def merge_tables(table1, table2, table3, table4):
    # Read the tables into pandas DataFrames using tab separator
    func_df = pd.read_csv(table1, sep='\t')
    annos_df = pd.read_csv(table2, sep='\t')
    counts_df = pd.read_csv(table3, sep='\t')
    mag_df = pd.read_csv(table4, sep='\t', names=['Contig_Name', 'MAG_ID'])

    func_df['orf_id'] = [x.split('_', 1)[1] for x in func_df['ORF_ID']]

    counts_df['contig_id'] = [x.rsplit('_', 1)[0] for x in counts_df['Contig']]
    counts_df['orf_id'] = [x.split('_', 1)[1] for x in counts_df['Contig']]

    # Merge the tables using the specified column
    merged1_df = func_df.merge(annos_df, on='orf_id', how='left')
    merged2_df = merged1_df.merge(counts_df, on='orf_id', how='left')
    merged3_df = merged2_df.merge(mag_df, on='Contig_Name', how='left')

    return merged3_df

def main():
    parser = argparse.ArgumentParser(description='Merge three tab-separated tables on a shared column.')
    parser.add_argument('func_annos', type=str, help='Path to the *.functional_and_taxonomic_table.txt')
    parser.add_argument('orf_annos', type=str, help='Path to the *.2.txt')
    parser.add_argument('orf_counts', type=str, help='Path to the *.orf_counts.txt')
    parser.add_argument('map_map', type=str, help='Path to the MS_contig_map.tsv')
    parser.add_argument('output_file', type=str, help='Path to the output merged table file')

    args = parser.parse_args()

    try:
        # Call the merge_tables function to merge the tables
        result_df = merge_tables(args.func_annos, args.orf_annos, args.orf_counts, args.map_map)
        result_df.to_csv(args.output_file, sep='\t', index=False)

        # Display the merged result
        print(result_df)

    except FileNotFoundError as e:
        print(f"Error: {e.filename} not found.")
    except pd.errors.EmptyDataError:
        print("Error: One or more input tables are empty.")
    except pd.errors.MergeError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
