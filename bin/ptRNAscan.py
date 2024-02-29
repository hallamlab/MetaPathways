#!/usr/bin/env python

import subprocess
import threading
import argparse
import signal
import os
import pyfastx
import shutil
import math
import queue
import sys

# Global list to keep track of subprocesses
processes = []
# Queue for file parts
file_parts_queue = queue.Queue()

def run_tRNAscan(thread_id, args):
    while not file_parts_queue.empty():
        try:
            input_file, output_files = file_parts_queue.get_nowait()
        except queue.Empty:
            break

        print(f"Thread {thread_id}: Processing {input_file}")
        cmd = ["tRNAscan-SE", input_file]

        # Add output file flags and paths
        for flag, file_path in output_files.items():
            cmd.extend([flag, file_path])

        # Include other flags as needed
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

        process = subprocess.Popen(cmd, preexec_fn=os.setsid)
        processes.append(process)
        process.wait()
        print(f"Thread {thread_id}: Completed processing {input_file}")
        file_parts_queue.task_done()

def filter_and_split_file(input_file, tmp_dir, min_length, num_threads, max_contigs_per_file=None):
    fasta = pyfastx.Fasta(input_file, build_index=False)
    contigs = [(name, seq) for name, seq in fasta if len(seq) >= min_length]

    num_files = min(math.ceil(len(contigs) / max_contigs_per_file) if max_contigs_per_file else len(contigs), num_threads)

    file_parts = []
    for i in range(num_files):
        part_filename = os.path.join(tmp_dir, f'{os.path.basename(input_file)}_{i}.fasta')
        with open(part_filename, 'w') as smallfile:
            for name, seq in contigs[i::num_files]:
                smallfile.write(f">{name}\n{seq}\n")
        file_parts.append(part_filename)

    return file_parts

def cleanup_processes():
    for process in processes:
        if process.poll() is None:
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)

def signal_handler(signum, frame):
    cleanup_processes()
    sys.exit(1)

def main(input_file, num_threads, tmp_dir, min_length, args):
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    print("Filtering and splitting input file...")
    file_parts = filter_and_split_file(input_file, tmp_dir, min_length, num_threads, args.split_by_seq)

    # Prepare output file paths for each part
    output_types = ["--output", "--struct", "--isospecific", "--stats", "--bed", "--gff", "--fasta", "--log"]
    for part_file in file_parts:
        output_files = {flag: getattr(args, flag.strip("--")) for flag in output_types if getattr(args, flag.strip("--"), None)}
        file_parts_queue.put((part_file, output_files))

    actual_threads = min(len(file_parts), num_threads)
    print(f"Starting {actual_threads} tRNAscan-SE threads...")
    threads = []
    for i in range(actual_threads):
        thread = threading.Thread(target=run_tRNAscan, args=(i, args))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    print(f"Processing complete. Outputs saved.")

    try:
        shutil.rmtree(tmp_dir)
        print(f"Successfully removed temporary directory: {tmp_dir}")
    except Exception as e:
        print(f"Error removing temporary directory: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run tRNA-Scan SE in parallel.")
    parser.add_argument("--input_file", required=True, help="Input FASTA file.")
    parser.add_argument("--num_threads", type=int, default=4, help="Number of threads to use.")
    parser.add_argument("--tmp_dir", required=True, help="Temporary directory for intermediate files.")
    parser.add_argument("--min_length", type=int, default=0, help="Minimum sequence length to consider.")
    parser.add_argument("--split_by_seq", type=int, help="Maximum number of sequences per file.", default=1000)

    # Define other arguments here
    parser.add_argument('-B', '--bacterial', action='store_true', help='Bacterial tRNA model')
    parser.add_argument('-A', '--archaeal', action='store_true', help='Archaeal tRNA model')
    parser.add_argument('-M', '--mito', help='Mitochondrial tRNA model')
    parser.add_argument('-O', '--general', action='store_true', help='General tRNA model')
    parser.add_argument('-G', '--genomic', action='store_true', help='Genomic tRNA model')
    parser.add_argument('-E', '--eukaryotic', action='store_true', help='Eukaryotic tRNA model')
    parser.add_argument('-I', '--infernal', action='store_true', help='Use Infernal for search')
    parser.add_argument('--max', action='store_true', help='Maximum sensitivity mode')
    parser.add_argument('-L', '--legacy', action='store_true', help='Use legacy search method')
    parser.add_argument('--cove', action='store_true', help='Search using COVE analysis only')
    parser.add_argument('--nopseudo', action='store_true', help='Disable pseudogene checking')
    parser.add_argument('--quiet', action='store_true', help='Disable stdout')

    # Add output file arguments as needed

    args = parser.parse_args()

    main(args.input_file, args.num_threads, args.tmp_dir, args.min_length, args)
