#!/bin/bash

ptoolpath=$1 # path to ../ptools input
pgdbpath=$2  # path to ../results/pgdb output
tag_id=$3 # organism/tag id

rm -rf /data/ptools-local # clear any previous runs
mkdir -p /data/ptools-local/pgdbs/user # create the required dirs and files
mkdir -p /data/blastdb
cp /opt/sandbox/ptools-init.dat /data/ptools-local/ptools-init.dat

Xvfb :${DISPLAY#*:} &
sleep 10

## This builds the PGDB:
/opt/pathway-tools/pathway-tools -patho ${ptoolpath} -no-taxonomic-pruning -no-web-cel-overview -no-patch-download -tip -no-cel-overview -disable-metadata-saving -nologfile
sleep 10

## Get the Org ID of the just-built PGDB:
org_id=$tag_id
sleep 10

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-web-cel-overview -tip -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"
sleep 10

sub_id=$(echo ${org_id} | tr '[:upper:]' '[:lower:]')
subpath=/data/ptools-local/pgdbs/user/${sub_id}cyc
tar -cjf ${pgdbpath}/${org_id}cyc.tar.bz2 -C ${subpath} .
sleep 10

rm -rf /data/ptools-local # clear present run
ps -ef | grep Xvfb | grep -v grep | awk '{print $2}' | xargs kill

