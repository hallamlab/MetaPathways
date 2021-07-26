FROM continuumio/miniconda3

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /opt

### Definitions:

ENV PYTHONPATH=/opt/mp_repo:/opt/mp_repo/libs


### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update -y 
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make \
    				   	   	      python3 \
						      zlib1g-dev \
						      python3-pip \
						      wget

# Install MetaPathways:
RUN pip3 install metapathways

## COPY over Makefile:
COPY Makefile        /opt/mp_repo/

## Set up Conda:
RUN make -C mp_repo conda-install-deps 


### Copying the repo files into the Docker image:
COPY extensions	      /opt/mp_repo/extensions/
COPY resources       /opt/mp_repo/resources/


RUN mkdir /opt/pgdb_dir

## Compile & Install Extensions:
RUN make -C mp_repo extensions-build
RUN make -C mp_repo extensions-install

## Make things work for Singularity by relaxing the permissions:
RUN chmod -R 755 /opt/mp_repo
RUN chmod -R 755 /opt/conda

### EntryPoint source:
## TODO