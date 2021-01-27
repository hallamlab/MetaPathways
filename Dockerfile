FROM continuumio/miniconda

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Definitions:

#ENV METAPATHWAYS_DB=/tmp/mp_db_dir
#ENV MP_DB_URI=https://www.dropbox.com/s/ye3kpve041e0r39/MetaPathways_DBs.zip

ENV PYTHONPATH=/root/mp_repo:/root/mp_repo/libs


### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make #python3 #zlib1g-dev

<<<<<<< HEAD
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install python3-pip
#RUN DEBIAN_FRONTEND=noninteractive apt-get -y install wget
RUN DEBIAN_FRONTEND=noninteractive pip3 install metapathways
=======
#RUN mkdir $METAPATHWAYS_DB						      
#RUN cd $METAPATHWAYS_DB
#RUN wget $MP_DB_URI

#RUN git clone --recurse-submodules https://taltman1@bitbucket.org/BCB2/metapathways_engcyc.git
>>>>>>> c796fa62db2b2b6e220829d0ec94e8c01cc90e40

### Copying the repo files into the Docker image:
COPY executables     /root/mp_repo/executables/
COPY resources       /root/mp_repo/resources/
COPY Makefile        /root/mp_repo/
COPY libs            /root/mp_repo/libs/
COPY MetaPathways.py /root/mp_repo/
COPY MetaPathwaysrc  /root/mp_repo/


#RUN cd /root/mp_repo && make pre-docker-builds
#RUN cd /root/mp_repo && make METAPATHWAYS_DB_FETCH
RUN touch /root/mp_repo/executables/linux/FGS+
RUN touch /root/mp_repo/executables/linux/ptools
RUN mkdir /root/pgdb_dir

### EntryPoint source:
##RUN . /root/mp_repo/MetaPathwaysrc
#RUN time ./MetaPathways.py -i regtests/input/B1.fasta -o test/B1_MPout/ -p test/mp_param.txt -c test/mp_config.txt
