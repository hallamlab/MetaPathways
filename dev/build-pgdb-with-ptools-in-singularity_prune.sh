#!/bin/bash

# $1 is the parent directory from MP3 that contains the ./ptools dir
#    it is assumed to contain ./results/pgdb as well
inpath="$1"
sifpath="$2"
srcpath="$3"
mpdpath="$4"
tmppath="$5"

singularity exec -B ${inpath}:/pt_inputs,${mpdpath}:/mpdevpath,${srcpath}:/pt_src,${tmppath}:/data $sifpath /mpdevpath/run-pathway-tools-and-copy-pgdb-singularity_prune.sh /pt_inputs /pt_src
