#!/bin/bash

ppath=$1
ptoolpath=${ppath}/ptools
pgdbpath=${ppath}/results/pgdb

mkdir -p /scratch/st-shallam-1/mcglock/tmp/data/ptools-local/pgdbs/user
mkdir -p /scratch/st-shallam-1/mcglock/tmp/data/blastdb
cp /home/ryan/Projects/BCB2/ptools-container/ptools-init.dat /scratch/st-shallam-1/mcglock/tmp/data/ptools-local/ptools-init.dat
cp /home/ryan/Projects/BCB2/ptools-container/.ncbirc $HOME/
export DATA_LOADERS=/scratch/st-shallam-1/mcglock/tmp/data/blastdb

## This builds the PGDB:
/opt/pathway-tools/pathway-tools -patho ${ptoolpath} -no-taxonomic-pruning -no-web-cel-overview -tip -no-patch-download

## Get the Org ID of the just-built PGDB:
org_id=`awk -F"\t" '$1 == "ID" { print $2 }' /scratch/st-shallam-1/mcglock/tmp/data/ptools-local/pgdbs/user/*cyc/1.0/input/organism.dat`

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-patch-download \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"

tar -cjf ${pgdbpath}/${org_id}cyc.tar.bz2  -C /scratch/st-shallam-1/mcglock/tmp/data/ptools-local/pgdbs/user .

#rm -rf /tmp/data
#rm -rf $HOME/.ncbirc
