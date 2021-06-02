#!/bin/bash

# $1 is the sample folder containing the 0.pf file 
# $2 is the output folder where pgdb.tar would be available
inpath=$1
outpath=$2

rm -rf ${inpath}/tmp

docker run \
    -v ${inpath}:/data \
    -v ${outpath}:/output \
    quay.io/mcglock/pt_v24.5 \
    run-pathway-tools-and-copy-pgdb.sh \
    -patho /data -no-taxonomic-pruning -no-web-cel-overview -tip
