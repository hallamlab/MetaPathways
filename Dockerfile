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

## Create the mp_repo directory, and copy over the Makefile
COPY Makefile        /opt/mp_repo/

## Set up Conda:
RUN make -C mp_repo conda-install-deps 

# Install MetaPathways:
RUN pip3 install git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways

### Copying the repo files into the Docker image:
COPY extensions	     /opt/mp_repo/extensions/
COPY resources       /opt/mp_repo/resources/
COPY tests           /opt/mp_repo/tests/

RUN mkdir /opt/pgdb_dir

## Compile & Install Extensions:
RUN make -C mp_repo extensions-build
RUN make -C mp_repo extensions-install


## Copy over Snakemake file & config file:
COPY Snakefile /opt/mp_repo/
COPY snakemake_config.yaml /opt/mp_repo


## Make things work for Singularity by relaxing the permissions:
RUN chmod -R 755 /opt/mp_repo
RUN chmod -R 755 /opt/conda

### EntryPoint source:
## TODO