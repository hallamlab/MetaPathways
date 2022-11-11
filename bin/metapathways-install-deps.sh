#!/bin/sh

echo "Installing MetaPathways dependencies using Conda:"
conda install --yes --channel conda-forge mamba
conda install --yes --channel conda-forge --channel bioconda coverm
mamba install --yes --channel conda-forge curl
mamba install --yes --channel bioconda blast prodigal bwa samtools barrnap trnascan-se
echo "Installation of MetaPathways dependencies complete!"
