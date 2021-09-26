FROM ubuntu:bionic

MAINTAINER Ryan J. McLaughlin, University of British Columbia

Workdir /root

### Install apt dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y make python3 #zlib1g-dev

RUN DEBIAN_FRONTEND=noninteractive apt-get -y install python3-pip
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install wget
#RUN DEBIAN_FRONTEND=noninteractive pip3 install metapathways

RUN DEBIAN_FRONTEND=noninteractive apt install -y zlib1g-dev

RUN DEBIAN_FRONTEND=noninteractive  apt-get -y install ncbi-blast+
RUN DEBIAN_FRONTEND=noninteractive  apt-get -y install samtools
RUN DEBIAN_FRONTEND=noninteractive  apt-get -y install bwa
RUN DEBIAN_FRONTEND=noninteractive  apt-get -y install prodigal

### Add source code for various executables MP3 requires
ADD ./ $HOME/root/metapathways_engcyc/

RUN cd /root/metapathways_engcyc/ && pip3 install .
RUN cd /root/metapathways_engcyc/extensions/FAST && make && cp fast* /usr/local/bin/.
RUN cd /root/metapathways_engcyc/extensions/metacount && make && cp metacount /usr/local/bin/.
RUN cd /root/metapathways_engcyc/extensions/trnascan && make && cp trnascan-1.4 /usr/local/bin/.