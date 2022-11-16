#!/bin/bash

ptoolpath=$1 # path to ../ptools for community-level
pgdbpath=$2  # path to ../results/pgdb
srcpath=$3   # path to ptools-container directory

mkdir -p /data/ptools-local/pgdbs/user
mkdir -p /data/blastdb
cp ${srcpath}/ptools-init.dat /data/ptools-local/ptools-init.dat

Xvfb :${DISPLAY#*:} &

## This builds the PGDB:
/opt/pathway-tools/pathway-tools -patho ${ptoolpath} -no-taxonomic-pruning -no-web-cel-overview -no-patch-download -tip -no-cel-overview -disable-metadata-saving -nologfile

## Get the Org ID of the just-built PGDB:
org_id=$(awk -F"\t" '$1 == "ID" { print $2 }' /data/ptools-local/pgdbs/user/*cyc/1.0/input/organism.dat)

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-web-cel-overview -tip -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"

sub_id=$(echo ${org_id} | tr '[:upper:]' '[:lower:]')
subpath=/data/ptools-local/pgdbs/user/${sub_id}cyc
tar -cjf ${pgdbpath}/${org_id}cyc.tar.bz2 -C ${subpath} .



