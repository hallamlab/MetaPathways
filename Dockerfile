FROM continuumio/miniconda

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y wget \
    				   	   	      ncbi-blast+ \
						      zlib1g-dev \
						      


ENV METAPATHWAYS_DB=/tmp/mp_db_dir
ENV MP_DB_URI=https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip

RUN cd $METAPATHWAYS_DB && wget $MP_DB_URI

RUN source ~/dev/metapathways_engcyc/MetaPathwaysrc