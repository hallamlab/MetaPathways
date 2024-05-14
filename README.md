[![Generic badge](https://img.shields.io/badge/Codebase-MetaPathways-<COLOR>.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/) [![Generic badge](https://img.shields.io/badge/Version-3.5-<COLOR>.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/)[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/) 

# MetaPathways

A master-worker model for environmental Pathway/Genome Database construction on grids and clouds

**Current Team:** Tomer Altman, Julia Anstett, Aria Hahn, Kishori M. Konwar, Tony Liu, Ryan McLaughlin, Aditi Nagaraj, and Steven J. Hallam

**Previous Team Members:** Niels W. Hanson and Shang-Ju Wu

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of information into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. MetaPathways v1.0 is a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

MetaPathways has since advanced as a modular tool, deepening our understanding of microbial metabolism at various biological levels. With this release, we have addressed previous challenges in modularity and database management. v3.5 enhances user accessibility through streamlined installation via package indexes or containers, refined modules, and interface upgrades. It boasts updated algorithm support for sequence feature prediction, annotation, metabolic inference, and coverage metrics. Tested on mock community data, Metapathways v3.5 demonstrates improved performance and usability. With automated installation and database management, this open-source tool makes advanced metagenomic analysis more accessible. Metapathways v3.5 represents a significant step forward in automated, comprehensive metagenomic analysis, facilitating a deeper exploration of microbial interactions and metabolic functions in environmental genomics.

## Quickstart

### Installation with Anacodna
```
conda install -c hallamlab -c bioconda -c conda-forge metapathways
```

### Installation with containers
```
singularity pull docker://quay.io/hallamlab/metapathways

docker pull quay.io/hallamlab/metapathways
```

### Execution
```{bash}
metapathways run --help
```

## [Documentation](https://metapathways.readthedocs.io/en/dev/)

## Support

[Technical questions, bug reports, and general inquires can be made here.](https://bitbucket.org/BCB2/metapathways/issues?status=new&status=open)

## Citation

If you use MetaPathways in your research, please cite the following article:

> Kishori M. Konwar, Niels W. Hanson, Antoine P Pagé, and Steven J. Hallam. *MetaPathways: a modular pipeline for constructing pathway/genome databases from environmental sequence information.* BMC Bioinformatics, 2013;14(1):202. [doi:10.1186/1471-2105-14-202](https://doi.org/10.1186/1471-2105-14-202)

