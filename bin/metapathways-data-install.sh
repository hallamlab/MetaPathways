#!/bin/bash

data_dir="$1"
command="$2"
metapathways_pkg_dir=`python3 -c "import metapathways; print(metapathways.__path__[0])"`

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --cores 1 \
	  --config ref_db_dir="$data_dir" \
	           pkg_dir="$metapathways_pkg_dir" \
	  --snakefile "$metapathways_pkg_dir/Snakefile" \
	  -- \
	  "$command"

conda deactivate
