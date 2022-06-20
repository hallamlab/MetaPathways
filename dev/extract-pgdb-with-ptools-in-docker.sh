# $1 is the sample

#userpgdbfolders=/home/kishori/pgdb_output/output/user
#socketfolder=/home/kishori/pgdb_output/tmp

userpgdbfolders=$1
socketfolder=$2

docker run \
    -v ${userpgdbfolders}:/opt/data/ptools-local/pgdbs/user \
    -v ${socketfolder}:/tmp \
    -v ${PWD}:/tools \
    quay.io/hallamlab/ptools-container:24.5-v1 \
    /tools/run-pathway-tools-to-extract.sh 
