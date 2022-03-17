#!/bin/bash

# $1 is the parent directory from MP3 that contains the ./ptools dir
#    it is assumed to contain ./results/pgdb as well
inpath="$1"
srcpath="$2"
tmppath="$3"
mpdpath="$4"

singularity exec ${srcpath}/ptools-dev.sif ${mpdpath}/run-pathway-tools-and-copy-pgdb-singularity.sh $inpath $srcpath $tmppath
