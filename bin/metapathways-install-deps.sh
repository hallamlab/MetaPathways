#!/bin/sh

echo "Installing MetaPathways dependencies using Conda:"
conda install --yes --channel conda-forge mamba
#mamba install --yes --channel conda-forge curl
#mamba install --yes --channel bioconda blast prodigal bwa samtools barrnap trnascan-se
mamba install --yes \
      --channel conda-forge \
      --channel bioconda \
      --name metapathways \
      curl blast prodigal bwa samtools barrnap trnascan-se snakemake
echo "Installation of MetaPathways dependencies complete!"
