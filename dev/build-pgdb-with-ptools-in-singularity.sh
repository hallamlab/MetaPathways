#!/bin/bash

# $1 is the parent directory from MP3 that contains the ./ptools dir
#    it is assumed to contain ./results/pgdb as well
inpath="$1"

singularity exec /home/ryan/Desktop/test_ptools_sif/ptools-dev.sif /home/ryan/Projects/BCB2/metapathways/dev/run-pathway-tools-and-copy-pgdb-singularity.sh $inpath
