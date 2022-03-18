#!/bin/bash

ppath=$1  # path to MP3 main output diretory
ptoolpath=${ppath}/ptools
pgdbpath=${ppath}/results/pgdb
srcpath=$2  # path to ptools-container directory
tmppath=$3


#rm -rf ${tmppath}/data
#mkdir -p ${tmppath}/data/ptools-local/pgdbs/user
#mkdir -p ${tmppath}/data/blastdb
#cp ${srcpath}/ptools-init.dat ${tmppath}/data/ptools-local/ptools-init.dat
#cp ${srcpath}/.ncbirc $HOME/

mkdir -p /data/ptools-local/pgdbs/user
mkdir -p /data/blastdb
cp ${srcpath}/ptools-init.dat /data/ptools-local/ptools-init.dat

Xvfb :${DISPLAY#*:} &

## This builds the PGDB:
/opt/pathway-tools/pathway-tools -patho ${ptoolpath} -no-taxonomic-pruning -no-web-cel-overview -tip -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile

## Get the Org ID of the just-built PGDB:
org_id=$(basename ${ppath})  # `awk -F"\t" '$1 == "ID" { print $2 }' ${tmppath}/data/ptools-local/pgdbs/user/*cyc/1.0/input/organism.dat`

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-web-cel-overview -tip -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"

sub_id=$(echo ${org_id} | cut -d'_' -f 2-)
tar -cjf ${pgdbpath}/${org_id}cyc.tar.bz2 -C ${tmppath}/data/ptools-local/pgdbs/user/*_${sub_id} .



