FROM continuumio/miniconda

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Install apt dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make python3 #zlib1g-dev

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install python3-pip
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install wget
RUN DEBIAN_FRONTEND=noninteractive pip3 install metapathways

RUN DEBIAN_FRONTEND=noninteractive apt install -y zlib1g-dev
RUN DEBIAN_FRONTEND=noninteractive  apt-get -y install ncbi-blast+

### Install binaries from C/C++ code
RUN wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/c_cpp_sources.1.0.tar.gz

RUN tar -zxvf c_cpp_sources.1.0.tar.gz && cd c_cpp_sources && make && make install

