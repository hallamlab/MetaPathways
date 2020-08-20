FROM continuumio/miniconda

MAINTAINER Tomer Altman, Altman Analytics LLC

Workdir /root

### Install apt dependencies

RUN DEBIAN_FRONTEND=noninteractive apt-get update
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y wget emboss samtools parallel bcftools tabix make infernal bowtie2 hmmer gawk


