#!/bin/bash

Xvfb $DISPLAY &

/opt/pathway-tools/pathway-tools "$@"

tar -cjf /output/pgdb.tar  -C /opt/data/ptools-local/pgdbs/user .

