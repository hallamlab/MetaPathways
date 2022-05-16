FROM ubuntu:focal

LABEL base.image="ubuntu:focal"
LABEL dockerfile.version="0"
LABEL software="MetaPathways + deps"
LABEL software.version="3.5"
LABEL description="[add description]"
LABEL website="[add website]"
LABEL license="[add license]"
LABEL maintainer="Ryan J. McLaughlin"
LABEL maintainer.email="mclaughlinr2@gmail.com"

Workdir /opt

### Definitions:

ENV PYTHONPATH=/opt/mp_repo:/opt/mp_repo/libs

# Install base utilities
RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make \
							python3 \
							zlib1g-dev \
							python3-pip \
							wget \
							git

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

## Create the mp_repo directory, and copy over the Makefile
COPY Makefile        /opt/mp_repo/

## Set up Conda:
RUN make -C /opt/mp_repo conda-install-deps 

# Install MetaPathways:
RUN pip3 install git+https://bitbucket.org/BCB2/metapathways.git@prokka_test#egg=MetaPathways

### Copying the repo files into the Docker image:
COPY resources       /opt/mp_repo/resources/
COPY tests           /opt/mp_repo/tests/
COPY extensions	     /opt/mp_repo/extensions/

RUN mkdir /opt/pgdb_dir

## Compile & Install Extensions:
#RUN make -C mp_repo extensions-build
RUN make -C /opt/mp_repo extensions-install
COPY extensions/FAST/fastal /usr/local/bin/fastal
COPY extensions/FAST/fastdb /usr/local/bin/fastdb
RUN chmod 755 /usr/local/bin/fastal
RUN chmod 755 /usr/local/bin/fastdb

## Copy over Snakemake file & config file:
COPY Snakefile /opt/mp_repo/
COPY snakemake_config.yaml /opt/mp_repo

## Hack so that Prokka can't run tbl2asn
COPY resources/tbl2asn.dummy /opt/conda/bin/tbl2asn
RUN chmod +x /opt/conda/bin/tbl2asn

## Make things work for Singularity by relaxing the permissions:
RUN chmod -R 755 /opt/mp_repo
RUN chmod -R 755 /opt/conda

