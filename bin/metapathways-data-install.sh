#!/bin/bash

print_help() {
  echo "Usage: $0 <data_dir> <db_type> [nthreads] [metacyc_user] [metacyc_pswd]"
  echo "Options:"
  echo "	data_dir - location to save database."
  echo "	db_type - version of database to build."
  echo "		  acceptable values are:"
  echo "			stage_fast_full"
  echo "			stage_fast_lite"
  echo "			stage_fast_noMeta"
  echo "			stage_blast_full"
  echo "			stage_blast_lite"
  echo "	nthreads - number of threads to use. [all]"
  echo "	metacyc_user - username for MetaCyc data download."
  echo "	metacyc_pswd - password for MetaCyc data download. "  
  echo "  -h, --help       Display this help message."
}

while getopts ":h" opt; do
  case $opt in
    h)
      print_help
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      print_help
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

if [ "$#" -lt 2 ] || [ "$#" -gt 5 ]; then
  echo "Error: Invalid number of arguments."
  print_help
  exit 1
fi

data_dir="$1"
db_type="$2"
nthreads="${3:-all}"
metacyc_user="${4:-}"
metacyc_pswd="${5:-}"

metapathways_pkg_dir=$(python -c "import site; print(site.getsitepackages()[0])")/metapathways/build_DBs

snakemake \
  -p \
  -s "$metapathways_pkg_dir/Snakefile" \
  --cores "$nthreads" \
  --until "$db_type" \
  --config ref_db_dir="$data_dir" \
  metacyc_username="$metacyc_user" \
  metacyc_password="$metacyc_pswd" \
  script_path="$metapathways_pkg_dir"

