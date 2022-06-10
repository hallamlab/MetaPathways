#!/bin/sh

inputpath=$1
sifpath=$2
srcpath=$3
mpdpath=$4
tmppath=$5

sampleid=$(basename ${inputpath} | cut -d'.' -f1)

sh ${mpdpath}/build-pgdb-with-ptools-in-singularity.sh ${inputpath} ${sifpath} ${srcpath} ${mpdpath} ${tmppath}

rm -rf ${tmppath}/*

tar -xf ${inputpath}/results/pgdb/*.tar.bz2 -C ${inputpath}/results/pgdb/

rm -rf ${inputpath}/results/pgdb/*.tar.bz2

python3 ${mpdpath}/make-camelot-file-generate-report-singularity.py ${inputpath}/results/pgdb/1.0/data ${inputpath}/results/pgdb ${inputpath}/results/pgdb/${sampleid}_pwy.tsv




