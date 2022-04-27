#!/bin/sh

inputpath=$1
sifpath=$2
srcpath=$3
mpdpath=$4
mappath=$5
tmppath=$6

orf_file=$(ls ${inputpath}/results/annotation_table/*.ORF_annotation_table.txt)
sampleid=$(basename ${orf_file} | cut -d'.' -f1)
metacyc_file=$(ls ${inputpath}/blast_results/*metacyc*.FASTout.parsed.txt)
python3 ${mpdpath}/format_pf.py ${inputpath}/ptools/0.pf ${orf_file} ${mappath}/KO2EC_mapping.tsv ${metacyc_file} ${mappath}/MetaCyc-monomer-rxn-pairs.tsv

sh ${mpdpath}/build-pgdb-with-ptools-in-singularity.sh ${inputpath} ${sifpath} ${srcpath} ${mpdpath} ${tmppath}

rm -rf ${tmppath}/*

tar -xf ${inputpath}/results/pgdb/*.tar.bz2 -C ${inputpath}/results/pgdb/

rm -rf ${inputpath}/results/pgdb/*.tar.bz2

python3 ${mpdpath}/make-camelot-file-generate-report-singularity.py ${inputpath}/results/pgdb/1.0/data ${inputpath}/results/pgdb ${inputpath}/results/pgdb/${sampleid}_pwy.tsv




