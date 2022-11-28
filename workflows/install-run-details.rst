Create and Run the SRAtoolkit Snakemake Workflow
*******************

This page contains instructions on how to install and run
the SRAtoolkit Snakemake workflow.

Create SRAtoolkit Conda Env
===========================

Assuming that you have all prerequisites satisfied, open a terminal in the MetaPathways repo:
::
	cd ./workflows/
	conda env create -f sratk_env.yml


Run the SRAtoolkit test
========================

This test will attempt to download 3 ecoli SRA accessions using `prefetch`.
It will then use `fasterq-dump` to extract the FASTQ files from the downloaded
`.sra` archives. Directories and files are auto-generated and removed as the
process completes. The final outputs are located in `SRA_gz` which will contain
all the extracted FASTQs that have been gzipped to save space.
::
	snakemake -p -s Snakefile.sra_download -j 3 --use-conda --conda-frontend conda
