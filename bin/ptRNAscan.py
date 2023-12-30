#!/usr/bin/env python

import subprocess
import threading
import argparse
import signal
import os
import pyfastx
import shutil

# Global list to keep track of subprocesses
processes = []


def run_tRNAscan(input_file, output_files, thread_id, args):
    #print(f"Thread {thread_id}: Processing {input_file}")
    cmd = ["tRNAscan-SE", input_file]

    for flag, file_path in output_files.items():
        cmd.extend([flag, file_path])

    if args.bacterial:
        cmd.append("-B")
    if args.archaeal:
        cmd.append("-A")
    if args.mito:
        cmd.extend(["-M", args.mito])
    if args.general:
        cmd.append("-O")
    if args.genomic:
        cmd.append("-G")

    if args.eukaryotic:
        cmd.append("-E")
    if args.infernal:
        cmd.append("-I")
    if args.max:
        cmd.append("--max")
    if args.legacy:
        cmd.append("-L")
    if args.cove:
        cmd.append("--cove")
    if args.nopseudo:
        cmd.append("--nopseudo")
    if args.quiet:
        cmd.append("--quiet")
    # ... Include other flags as needed

    process = subprocess.Popen(cmd)
    processes.append(process)
    process.wait()
    #print(f"Thread {thread_id}: Completed processing {input_file}")

def cat_files(file_parts, stats_output_file):
    with open(stats_output_file, 'w') as outfile:
        for infile_path in file_parts:
            with open(infile_path) as infile:
                outfile.writelines(infile)

def filter_and_split_file(input_file, num_chunks, tmp_dir, min_length):
    fasta = pyfastx.Fasta(input_file, build_index=False)
    total_sequences = sum(1 for rec in fasta if len(rec[1]) >= min_length)
    lines_per_file = total_sequences // num_chunks
    extra_lines = total_sequences % num_chunks

    # Ensure at least one sequence per file if the total number is too small
    if lines_per_file <= 0:
        lines_per_file = 1
        num_chunks = total_sequences

    file_parts = []
    count = 0
    file_index = 0
    smallfile = None

    for name, seq in fasta:
        if len(seq) >= min_length:
            if count % lines_per_file == 0 and file_index < num_chunks:
                if smallfile:
                    smallfile.close()
                part_filename = os.path.join(tmp_dir, f'{os.path.basename(input_file)}_{file_index}.fasta')
                smallfile = open(part_filename, "w")
                file_parts.append(part_filename)
                file_index += 1
            smallfile.write(f">{name}\n{seq}\n")
            count += 1

            # When the regular files have reached their capacity, increase the limit to account for extra lines
            if file_index == num_chunks and extra_lines > 0:
                lines_per_file += 1
                extra_lines -= 1

            #print(f"Writing sequence: {name}, length: {len(seq)}")  # Logging statement

    if smallfile:
        smallfile.close()

    #print(f"Total sequences written: {count}")  # Logging total count
    #print(f"Total chunks created: {len(file_parts)}")  # Logging total count

    return file_parts

def combine_output_files(file_parts, output_file):
    total_lines_written = 0
    header_written = False
    with open(output_file, 'w') as outfile:
        for i, infile_path in enumerate(file_parts):
            lines_from_file = 0
            header_section = True
            header = []
            with open(infile_path) as infile:
                for line in infile:
                    if header_section:
                        header.append(line)
                        if '--------' in line:  # Dashed line marks the end of the header
                            header_section = False
                            if not header_written:
                                outfile.write(''.join(header))
                                header_written = True
                    elif not header_section:
                        outfile.write(line)
                        lines_from_file += 1
            
            #print(f"Processed {lines_from_file} lines from {infile_path}")
            total_lines_written += lines_from_file

    #print(f"Total lines written to final output: {total_lines_written}")

def cleanup_processes():
    #print("Cleaning up running processes...")
    for process in processes:
        if process.poll() is None:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)

def signal_handler(signum, frame):
    #print("\nSignal received, cleaning up...")
    cleanup_processes()
    exit(1)

def count_contigs(input_file, min_length):
    fasta = pyfastx.Fasta(input_file, build_index=False)
    return sum(1 for rec in fasta if len(rec[1]) >= min_length)

def main(input_file, num_threads, tmp_dir, min_length, args):
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Count the number of contigs
    num_contigs = count_contigs(input_file, min_length)

    # Adjust the number of threads if necessary
    if num_contigs < num_threads:
        num_threads = num_contigs

    #print("Filtering and splitting input file...")
    file_parts = filter_and_split_file(input_file, num_threads, tmp_dir, min_length)

    #print(f"Starting {num_threads} tRNAscan-SE threads...")
    threads = []
    output_types = ["--output", "--struct", "--isospecific", "--stats", "--bed", "--gff", "--fasta", "--log"]
    output_files = {otype: getattr(args, otype.split("--")[1]) for otype in output_types if getattr(args, otype.split("--")[1])}
    
    for i, part_file in enumerate(file_parts):
        part_output_files = {flag: os.path.join(tmp_dir, f"{os.path.basename(input_file)}_{flag.strip('-')}_{i}.txt") for flag in output_files}
        thread = threading.Thread(target=run_tRNAscan, args=(part_file, part_output_files, i, args))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()
    
    # Combining files for each output type
    for otype, ofile in output_files.items():
        part_files = [os.path.join(tmp_dir, f"{os.path.basename(input_file)}_{otype.strip('-')}_{i}.txt") for i in range(num_threads)]
        combine_function = cat_files if otype != "--output" else combine_output_files
        combine_function(part_files, ofile)

    #print(f"Processing complete. Outputs saved.")
    
    try:
        shutil.rmtree(tmp_dir)
        #print(f"Successfully removed temporary directory: {tmp_dir}")
    except Exception as e:
        print(f"Error removing temporary directory: {e}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run tRNA-Scan SE in parallel.")
    # Add arguments for various output types
    parser.add_argument("--input_file", help="Input FASTA file.")
    parser.add_argument("--num_threads", type=int, help="Number of threads to use.")
    parser.add_argument("--tmp_dir", help="Temporary directory for intermediate files.")
    parser.add_argument("--min_length", type=int, default=0, help="Minimum sequence length to consider.")

    parser.add_argument("-o", "--output", help="File for final tRNA results")
    parser.add_argument("--struct", help="File for tRNA secondary structures")
    parser.add_argument("--isospecific", help="File for isotype-specific model results")
    parser.add_argument("-m", "--stats", help="File for statistics summary")
    parser.add_argument("--bed", help="File for BED format results")
    parser.add_argument("-j", "--gff", help="File for GFF3 format results")
    parser.add_argument("-a", "--fasta", help="File for FASTA format of predicted tRNA sequences")
    parser.add_argument("--log", help="File for log of program progress")
    parser.add_argument('-B', '--bacterial', action='store_true', help='Bacterial tRNA model')
    parser.add_argument('-A', '--archaeal', action='store_true', help='Archaeal tRNA model')
    parser.add_argument('-M', '--mito', nargs='?', const='mammal', default=False, help='Mitochondrial tRNA model')
    parser.add_argument('-O', '--general', action='store_true', help='General tRNA model')
    parser.add_argument('-G','--genomic', action='store_true', help='Genomic tRNA model')
    parser.add_argument('-E', '--eukaryotic', action='store_true', help='Eukaryotic tRNA model')
    parser.add_argument('-I', '--infernal', action='store_true', help='Use Infernal for search')
    parser.add_argument('--max', action='store_true', help='Maximum sensitivity mode')
    parser.add_argument('-L', '--legacy', action='store_true', help='Use legacy search method')
    parser.add_argument('--cove', action='store_true', dest='cove', help='Search using COVE analysis only')
    parser.add_argument('--nopseudo', action='store_true', help='Disable pseudogene checking')
    parser.add_argument('--quiet', action='store_true', help='Disable stdout')
    
    args = parser.parse_args()
    
    # Create temporary and output directories if they don't exist
    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    # Call the main function with parsed arguments
    main(args.input_file, args.num_threads, args.tmp_dir, args.min_length, args)