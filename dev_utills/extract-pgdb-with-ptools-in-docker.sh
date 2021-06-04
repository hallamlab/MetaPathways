# $1 is the sample

#userpgdbfolders=/home/kishori/pgdb_output/output/user
#socketfolder=/home/kishori/pgdb_output/tmp

userpgdbfolders=$1
socketfolder=$2

#docker run \
#    -v ${inpath}:/data quay.io/mcglock/pt_v24.5 \
#    run-pathway-tools.sh \
#    -patho /data -no-taxonomic-pruning -no-web-cel-overview -tip

docker run \
    -v ${userpgdbfolders}:/opt/data/ptools-local/pgdbs/user \
    -v ${socketfolder}:/tmp \
    -v ${PWD}:/tools \
    quay.io/mcglock/pt_v24.5 \
    /tools/run-pathway-tools-to-extract.sh 
