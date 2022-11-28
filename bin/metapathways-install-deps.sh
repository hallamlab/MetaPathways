#!/bin/sh

echo "Installing MetaPathways dependencies using Conda:"
conda install --yes --channel conda-forge mamba
mamba install --yes --channel conda-forge --channel bioconda curl ncurses pysam pybedtools coverm blast prodigal bwa samtools barrnap trnascan-se
echo "Installation of MetaPathways dependencies complete!"
