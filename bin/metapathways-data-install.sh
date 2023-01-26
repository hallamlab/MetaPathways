#!/bin/bash

data_dir="$1"
command="$2"

metapathways_pkg_dir=$(python -c "import site; print(site.getsitepackages()[0])")/metapathways

# eval "$(conda shell.bash hook)"
# conda activate metapathways

snakemake --cores 1 --scheduler greedy \
	  --config ref_db_dir="$data_dir"  \
	           pkg_dir="$metapathways_pkg_dir" \
	  --snakefile "$metapathways_pkg_dir/Snakefile" \
	  -- \
	  "$command"

# conda deactivate
