#!/bin/bash

Xvfb :${DISPLAY#*:} &

## This builds the PGDB:
/opt/pathway-tools/pathway-tools "$@"

## Get the Org ID of the just-built PGDB:
org_id=`awk -F"\t" '$1 == "ID" { print $2 }' /opt/data/ptools-local/pgdbs/user/*cyc/1.0/input/organism.dat`

## This gets PTools to dump out the flat-files of the PGDB:
/opt/pathway-tools/pathway-tools \
    -no-patch-download \
    -eval "(progn (with-organism (:org-id '$org_id) (dump-frames-to-attribute-value-files (org-data-dir)))(exit))"

tar -cjf /output/${org_id}cyc.tar.bz2  -C /opt/data/ptools-local/pgdbs/user .


