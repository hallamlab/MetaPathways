# Steps to setting up the EggNOG DB for MetaPathways

1. Download required files from EggNOG. As of 01/11/2023, the FTP link is: `http://eggnog5.embl.de/download/eggnog_5.0/`
	Files to download:
		   e5.proteomes.faa
		   e5.sequence_aliases.tsv
		   e5.taxid_info.tsv

2. Create the `e5-names.txt` for EggNOG.
	Command: `python dev/build_DBs/eggnog/eggnog_db_prep.py --en_faa [path to]/e5.proteomes.faa --en_aliase [path to]/e5.sequence_aliases.tsv --out_dir [output dir]`
	Note: By Default anything before the first '.' is used as the DB name by MP. You might want to rename the proteomes to something else depending on your needs or RefDB naming conventions. example: e5.proteomes.faa -> eggnog-v5.faa

3. Build the EggNOG taxmap.
	Command: `python dev/build_DBs/eggnog/eggnog_build_taxmap.py --taxid_info [path to]/e5.taxid_info.tsv --out_dir [output dir]`
	Note: This creates a file named 'ncbi_taxonomy_tree.txt', which is used by MPs LCA function to assign taxonomy to functional annotations. Don't run this if you don't plan to use EggNOG for LCA. Instead there is another taxmap script for general use available.
