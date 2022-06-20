#!/bin/sh
### The purpose of this scripts is to run the Ptools Singularity container on MP3 output
### to create ePGDBs for a sample.

inputpath=$1 # MP3 output directory
sifpath=$2 # path to the Ptools SIF file
srcpath=$3 # path to `ptools-container` repository (uses files within this dir)
mpdpath=$4 # path to `metapathways/dev` dir within MP3 repo (contains convenience scripts)
tmppath=$5 # tmp working dir for Ptools to save intermediates

sampleid=$(basename ${inputpath} | cut -d'.' -f1) # creates sample ID from MP3 output dir

#sh ${mpdpath}/build-pgdb-with-ptools-in-singularity.sh ${inputpath} ${sifpath} ${srcpath} ${mpdpath} ${tmppath} # runs Ptools using `${inputpath}/ptools` file as inputs

#rm -rf ${tmppath}/* # clears tmp dir after ePGDB has been dumped into MP3 output dir

#tar -xf ${inputpath}/results/pgdb/*.tar.bz2 -C ${inputpath}/results/pgdb/ # uncompresses ePGDB tar

#rm -rf ${inputpath}/results/pgdb/*.tar.bz2 # removes tar archive

python3 ${mpdpath}/make-camelot-file-generate-report-singularity.py ${inputpath}/results/pgdb/1.0/data ${inputpath}/results/pgdb ${inputpath}/results/pgdb/${sampleid}_pwy.tsv # creates a TSV file of the pathways that were predicted by PATHOLOGIC within Ptools




