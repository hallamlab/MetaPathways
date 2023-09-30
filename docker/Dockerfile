FROM condaforge/mambaforge as build-env
# https://mamba.readthedocs.io/en/latest/user_guide/mamba.html
MAINTAINER Tomer Altman, Altman Analytics LLC

### EXAMPLES ###

### Build dev branch
# sudo docker build --network=host -t metapathways:dev .

### Build test branch
# sudo docker build --build-arg git_branch=test --network=host -t metapathways:test .

################

# Create the environment:
# packages from `metapathways-install-deps.sh` are included in the yml
COPY ./docker/conda_env.yml /opt/
RUN mamba env create --no-default-packages -f /opt/conda_env.yml

# move to clean execution environment
# which includes dependencies like wget & libtinfo6
FROM ubuntu
COPY --from=build-env /opt/conda/envs/metapathways /opt/conda/envs/metapathways
ENV PATH /opt/conda/envs/metapathways/bin:$PATH

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# Install MetaPathways and dependencies
# use the local install so local builds are up to date
# !! remember to run `make create-package` at least once !!
COPY ./dist/*.gz /opt/pip_dist/
RUN pip3 install /opt/pip_dist/*.gz \
    && echo "Make sure MetaPathways is installed:" \
    && MetaPathways -h

## We do some umask munging to avoid having to use chmod later on,
## as it is painfully slow on large directores in Docker.
RUN old_umask=`umask` && \
    umask 0000 && \
    umask $old_umask

# singularity doesn't use the -s flag, and that causes warnings
ENTRYPOINT ["/tini", "-s", "--"]
