#!/bin/bash

# $1 is the parent directory from MP3 that contains the ./ptools dir
#    it is assumed to contain ./results/pgdb as well
inpath="$1"
sifpath="$2"
srcpath="$3"
mpdpath="$4"

singularity exec -B ${TMPDIR}:/data $sifpath ${mpdpath}/run-pathway-tools-and-copy-pgdb-singularity.sh $inpath $srcpath
