#!/bin/bash

# $1 is the sample folder containing the 0.pf file 
# $2 is the output folder where pgdb.tar would be available
inpath="$1"
outpath="$2"

rm -rf ${inpath}/tmp

docker run \
    -v ${inpath}:/data \
    -v ${outpath}:/output \
    -v ${PWD}:/code \
    --network=host \
    quay.io/hallamlab/ptools-container:24.5-v1 \
    /code/run-pathway-tools-and-copy-pgdb.sh \
    -patho /data -no-taxonomic-pruning -no-web-cel-overview -tip -no-patch-download
