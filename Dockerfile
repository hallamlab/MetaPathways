FROM continuumio/miniconda3
MAINTAINER Tomer Altman, Altman Analytics LLC

### EXAMPLES ###

### Build dev branch
# sudo docker build --network=host -t metapathways:dev .

### Build test branch
# sudo docker build --build-arg git_branch=test --network=host -t metapathways:test .

################

Workdir /opt

### Definitions:

ENV PYTHONPATH=/opt/mp_repo:/opt/mp_repo/libs

ARG git_branch=dev

### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update -y 
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y python3 \
						      python3-pip \
						      wget

## Create the mp_repo directory, and copy over the Makefile
RUN mkdir /opt/mp_repo/


# Install MetaPathways:
RUN pip3 install git+https://bitbucket.org/BCB2/metapathways.git@${git_branch}#egg=MetaPathways

## Set up Conda:
## We do some umask munging to avoid having to use chmod later on,
## as it is painfully slow on large directores in Docker.
RUN old_umask=`umask` && \
    umask 0000 && \
    metapathways-install-deps.sh && \
    umask $old_umask

RUN mkdir /opt/pgdb_dir

## Make things work for Singularity by relaxing the permissions:
RUN chmod -R 755 /opt/mp_repo
#RUN chmod -R 755 /opt/conda

### EntryPoint source:
## TODOFROM continuumio/miniconda3
MAINTAINER Tomer Altman, Altman Analytics LLC

### EXAMPLES ###

### Build dev branch
# sudo docker build --network=host -t metapathways:dev .

### Build test branch
# sudo docker build --build-arg git_branch=test --network=host -t metapathways:test .

################

Workdir /opt

### Definitions:

ENV PYTHONPATH=/opt/mp_repo:/opt/mp_repo/libs

ARG git_branch=dev

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
## We do some umask munging to avoid having to use chmod later on,
## as it is painfully slow on large directores in Docker.
RUN old_umask=`umask` && \
    umask 0000 && \
    make -C mp_repo conda-install-deps && \
    umask $old_umask

# Install MetaPathways:
RUN pip3 install git+https://bitbucket.org/BCB2/metapathways.git@${git_branch}#egg=MetaPathways

RUN mkdir /opt/pgdb_dir

## Make things work for Singularity by relaxing the permissions:
RUN chmod -R 755 /opt/mp_repo
#RUN chmod -R 755 /opt/conda

### EntryPoint source:
## TODO
