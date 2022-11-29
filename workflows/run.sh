# requires conda installed locally on sockeye
# requires srtk env active

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export XDG_CACHE_HOME=$SCRIPT_DIR/cache

snakemake -p -s Snakefile.sra_download \
    --cluster 'qsub -l {resource_list} -N {params.name} -A {account} -m {mail_opt} -M {email} -V -o {params.pbs_out} -e {params.pbs_err}' \
    -j 10 --keep-going \
    --use-conda --conda-frontend conda

rm -r ./cache/*
