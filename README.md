# MetaPathways

Functional and taxonomic annotation of environmental genomes, with community- and population-level pathway inference.

## Quickstart

**Linux x86-64 · Conda/Mamba · included K12 example.** Install and run in a writable environment with internet access:

```bash
mamba create -n metapathways_env --override-channels \
  -c hallamlab -c conda-forge -c bioconda metapathways=3.5.1
mamba activate metapathways_env
metapathways --version

mkdir -p ~/metapathways-review
cd ~/metapathways-review
metapathways build_db --test
metapathways run --test
```

**Check the result:** outputs are in `test/k12_test/`. Inspect
`metapathways_steps_log.txt` for successful stages and `errors_warnings_log.txt`
for problems; do not rely on the exit code alone. Key outputs are
`results/annotation_table/k12_test.functional_and_taxonomic_table.txt`,
`genbank/k12_test.gbk`, and `results/rpkm/k12_test.contig_counts.tsv`.

The example includes small K12 FASTA/FASTQ and SwissProt/SILVA fixtures. Database
preparation downloads ExPASy enzyme records and NCBI taxonomy and writes indexes
inside the installed package. This is an installation test; it does not reproduce
the manuscript's CAMI2 benchmark or run Pathway Tools.

[Docker / Apptainer instructions](docker/README.quay.md) ·
[Release downloads](https://github.com/hallamlab/MetaPathways/releases) ·
[Full usage](https://metapathways.readthedocs.io/en/latest/usage.html) ·
[Benchmark provenance](docs/src/reproducibility.rst)

[![Version 3.5.1](https://img.shields.io/badge/Version-3.5.1-blue.svg)](https://github.com/hallamlab/MetaPathways/releases)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Install this source revision

For development on Linux x86-64 with Python 3.11:

```bash
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
mamba env create -f docker/conda_base.yml
mamba activate metapathways
mamba install -c conda-forge pip wheel 'setuptools>=83,<85'
python -m pip install --no-deps --no-build-isolation .
```

Tagged releases validate packages and containers through the
[release workflow](docs/releasing.md). Record the release build or container digest
used for a review. Optional MAG splitting and Pathway Tools commands require
additional dependencies; see the [full usage documentation](https://metapathways.readthedocs.io/en/latest/usage.html).

## Inputs, outputs, and reproducibility

`metapathways run` accepts nucleotide FASTA or amino-acid FASTA (`--input_format
fasta-amino`). Paired or interleaved FASTQ reads can be supplied for abundance
estimation. The current CLI does not accept GFF or GenBank as primary inputs.

Results include annotation tables, annotated GFF/GenBank files, predicted
sequences, abundance tables when reads are supplied, run statistics, and Pathway
Tools input files. Reference databases are prepared separately with
`metapathways build_db`; see [reproducibility notes](docs/src/reproducibility.rst)
for dependency versions, output locations, and the manuscript benchmark provenance.
The small included example is an installation test, not a reproduction of the
manuscript's CAMI2 performance benchmark.

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of information into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. MetaPathways v1.0 is a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

MetaPathways has since advanced as a modular tool, deepening our understanding of microbial metabolism at various biological levels. With this release, we have addressed previous challenges in modularity and database management. v3.5 enhances user accessibility through streamlined installation via package indexes or containers, refined modules, and interface upgrades. It boasts updated algorithm support for sequence feature prediction, annotation, metabolic inference, and coverage metrics. Tested on mock community data, Metapathways v3.5 demonstrates improved performance and usability. With automated installation and database management, this open-source tool makes advanced metagenomic analysis more accessible. Metapathways v3.5 represents a significant step forward in automated, comprehensive metagenomic analysis, facilitating a deeper exploration of microbial interactions and metabolic functions in environmental genomics.

## Team and repository

**Current Team:** Ryan J. McLaughlin, Tony X. Liu, Tomer Altman, Aditi N. Nallan, Aria S. Hahn, Julia Anstett, Connor Morgan-Lang, Kishori M. Konwar, and Steven J. Hallam

**Previous Team Members:** Niels W. Hanson and Shang-Ju Wu

The canonical source is [hallamlab/MetaPathways](https://github.com/hallamlab/MetaPathways). Earlier code is preserved separately in [MetaPathways-legacy](https://github.com/hallamlab/MetaPathways-legacy).

## [Documentation](https://metapathways.readthedocs.io/en/latest/)

## Support

[Technical questions, bug reports, and general inquires can be made here.](https://github.com/hallamlab/MetaPathways/issues)

## Citation

If you use MetaPathways in your research, please cite the following article:

> Ryan J. McLaughlin, Tony X. Liu, Tomer Altman, Aditi N. Nallan, Aria S. Hahn, Julia Anstett, Connor Morgan-Lang, Kishori M. Konwar, Steven J. Hallam. *MetaPathways v3.5: Modularity and Scalability Improvements for Pathway Inference from Environmental Genomes* bioRxiv (2024): 2024-06. [doi: https://doi.org/10.1101/2024.06.04.597460](https://doi.org/10.1101/2024.06.04.597460)


## Maintainer releases

Use [the release controller and CI workflow](docs/releasing.md) to set a version,
build and test packages, and publish downloadable GitHub releases and optional
Anaconda.org packages. The source version is declared in `metapathways/_version.py`.
