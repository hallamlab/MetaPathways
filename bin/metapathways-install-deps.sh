#!/bin/sh

echo "Installing MetaPathways dependencies using Conda:"
conda install --yes -c conda-forge mamba
mamba install --yes -c conda-forge curl
mamba install --yes -c bioconda blast prodigal bwa samtools barrnap trnascan-se
mamba create --yes -c conda-forge -c bioconda -n snakemake snakemake
echo "Installation of MetaPathways dependencies complete!"
