# This describes how to build the current MetaPathways reference database
# NOTE: You must have a working username and password for downloading the MetaCyc flatfiles

# Build MetaCyc, SwissProt, and SILVA
snakemake -p -s Snakefile --use-conda --cores all --until stage_fast_lite --config metacyc_username=<username> metacyc_password=<password> ref_db_dir=</path/to/refDB>

# Build SwissProt and SILVA, but skip MetaCyc
snakemake -p -s Snakefile --use-conda --cores all --until stage_fast_noMeta ref_db_dir=</path/to/refDB>

# Build Larger RefDB that includes: MetaCyc, SwissProt, CAZymes, UniRef90, and SILVA
# Warning, this is huge! Do not build unless you really really want it :)
snakemake -p -s Snakefile --use-conda --cores all --until stage_fast_full --config metacyc_username=<username> metacyc_password=<password> ref_db_dir=</path/to/refDB>

