FROM continuumio/miniconda3

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Definitions:

ENV PYTHONPATH=/root/mp_repo:/root/mp_repo/libs


### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update -y 
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make \
    				   	   	      python3 \
						      zlib1g-dev \
						      python3-pip \
						      wget
RUN DEBIAN_FRONTEND=noninteractive pip3 install metapathways

### Copying the repo files into the Docker image:
COPY executables     /root/mp_repo/executables/
COPY resources       /root/mp_repo/resources/
COPY Makefile        /root/mp_repo/
COPY libs            /root/mp_repo/libs/
#COPY MetaPathways.py /root/mp_repo/
#COPY MetaPathwaysrc  /root/mp_repo/


RUN touch /root/mp_repo/executables/linux/FGS+
RUN touch /root/mp_repo/executables/linux/ptools
RUN mkdir /root/pgdb_dir


## Set up Conda:
RUN make -f /root/mp_repo/Makefile conda-install-deps 

### EntryPoint source:
##RUN . /root/mp_repo/MetaPathwaysrc
#RUN time ./MetaPathways.py -i regtests/input/B1.fasta -o test/B1_MPout/ -p test/mp_param.txt -c test/mp_config.txt
