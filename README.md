[![Generic badge](https://img.shields.io/badge/Codebase-MetaPathways-blue.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/) [![Generic badge](https://img.shields.io/badge/Version-3.5.0-blue.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/) [![Python 3.10](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/) [![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

# MetaPathways

The canonical MetaPathways v3.5 source repository is [hallamlab/MetaPathways](https://github.com/hallamlab/MetaPathways). The earlier GitHub codebase is maintained separately as [MetaPathways-legacy](https://github.com/hallamlab/MetaPathways-legacy); it is not part of this repository.

A master-worker model for environmental Pathway/Genome Database construction on grids and clouds

**Current Team:** Ryan J. McLaughlin, Tony X. Liu, Tomer Altman, Aditi N. Nallan, Aria S. Hahn, Julia Anstett, Connor Morgan-Lang, Kishori M. Konwar, and Steven J. Hallam

**Previous Team Members:** Niels W. Hanson and Shang-Ju Wu

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of information into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. MetaPathways v1.0 is a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

MetaPathways has since advanced as a modular tool, deepening our understanding of microbial metabolism at various biological levels. With this release, we have addressed previous challenges in modularity and database management. v3.5 enhances user accessibility through streamlined installation via package indexes or containers, refined modules, and interface upgrades. It boasts updated algorithm support for sequence feature prediction, annotation, metabolic inference, and coverage metrics. Tested on mock community data, Metapathways v3.5 demonstrates improved performance and usability. With automated installation and database management, this open-source tool makes advanced metagenomic analysis more accessible. Metapathways v3.5 represents a significant step forward in automated, comprehensive metagenomic analysis, facilitating a deeper exploration of microbial interactions and metabolic functions in environmental genomics.

## Quickstart

### Installation with Conda/Mamba
```
mamba create -n metapathways_env -c hallamlab -c bioconda -c conda-forge metapathways
```

### Installation with containers
```
apptainer pull docker://quay.io/hallamlab/metapathways

docker pull quay.io/hallamlab/metapathways
```

### Install this source revision

The Conda package and container above are distributed separately from GitHub releases.
To run this checkout on Linux with Python 3.10:

```bash
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
mamba env create -f docker/conda_base.yml
mamba activate metapathways
python -m pip install --no-deps .
```

Optional MAG splitting and Pathway Tools commands also need MAGSplitter,
camelot-frs, and a separately installed/licensed Pathway Tools. See the
[full usage documentation](https://metapathways.readthedocs.io/en/latest/usage.html).

### Run the included example

```bash
metapathways build_db --test
metapathways run --test
```

The example uses the bundled K12 FASTA/paired FASTQ inputs and small SwissProt/SILVA
reference fixtures. Database preparation requires internet access to download
ExPASy enzyme records and NCBI taxonomy, and writes indexes into the installed
package's test database directory. Run it in a writable environment. Results go to
`./test/`; inspect the sample's `metapathways_steps_log.txt` and
`errors_warnings_log.txt` as well as the command output.

### Inputs, outputs, and reproducibility

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
