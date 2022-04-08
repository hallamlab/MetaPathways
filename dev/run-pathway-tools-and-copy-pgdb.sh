#!/bin/bash

## This starts a virtual frame buffer that some of the libraries that PTools is reliant on needs to be able to access:
Xvfb $DISPLAY &

## This builds the PGDB:
/opt/pathway-tools/pathway-tools "$@"

## Get the Org ID of the just-built PGDB:
org_id=`awk -F"\t" '$1 == "ID" { print $2 }' /opt/data/ptools-local/pgdbs/user/*cyc/1.0/input/organism.dat`

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-patch-download \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"

tar -cjf /output/pgdb.tar.bz2  -C /opt/data/ptools-local/pgdbs/user ${org_id}cyc

