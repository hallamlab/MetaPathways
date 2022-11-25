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
						      zlib1g-dev \
						      liblzma-dev \
						      libbz2-dev \
						      wget
RUN DEBIAN_FRONTEND=noninteractive apt-get install libtinfo6


## Create the mp_repo directory, and copy over the Makefile
RUN mkdir /opt/mp_repo/
RUN mkdir /opt/pgdb_dir

# Create the environment:
RUN conda create -n metapathways python=3.10

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "metapathways", "/bin/bash", "-c"]

# Install MetaPathways and dependencies
RUN pip3 install git+https://bitbucket.org/BCB2/metapathways.git@${git_branch}#egg=MetaPathways
RUN metapathways-install-deps.sh
# Demonstrate the environment is activated:
RUN echo "Make sure MetaPathways is installed:"
RUN MetaPathways -h

## We do some umask munging to avoid having to use chmod later on,
## as it is painfully slow on large directores in Docker.
RUN old_umask=`umask` && \
    umask 0000 && \
    umask $old_umask

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "metapathways"]





