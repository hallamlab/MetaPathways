[![Generic badge](https://img.shields.io/badge/Codebase-MetaPathways-<COLOR>.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/) [![Generic badge](https://img.shields.io/badge/Version-2.5-<COLOR>.svg)](https://hallam.microbiology.ubc.ca/MetaPathways/) [![Generic badge](https://img.shields.io/badge/Python Version->2.7-<COLOR>.svg)](https://www.python.org/) [![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/) 

# MetaPathways

A master-worker model for environmental Pathway/Genome Database construction on grids and clouds

**Current Team:** Tomer Altman, Aria Hahn, Kishori M. Konwar, Ryan McLaughlin, and Steven J. Hallam

**Previous Team Members:** Niels W. Hanson and Shang-Ju Wu

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of infor- mation into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. We recently introduced MetaPathways v1.0, a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

Here we present MetaPathways v3.0, incorporating a graphical user interface (GUI) and refined task management methods. The MetaPathways GUI provides an intuitive display for setup and process monitoring and supports interactive data visualization and sub-setting via a custom Knowledge Engine data structure. A master-worker model is adopted for task management allowing users to scavenge computational results from a number of worker grids in an ad hoc, asynchronous, distributed network that dramatically increases fault tolerance. This model facilitates the use of EC2 instances extending ePGDB construction to the Amazon Elastic Cloud.

## Installation

MetaPathways v3.0 requires Python 3.0 or greater. For full functionality, you should also install [Pathway Tools](http://bioinformatics.ai.sri.com/ptools/), developed by SRI International.

Please see the [MetaPathways v3.0 documentation](https://metapathways.readthedocs.io/en/dev/) for installation details.


## Citation

If you use MetaPathways in your research, please cite the following article:

> Niels W. Hanson, Kishori M. Konwar, Shang-Ju Wu, Steven J. Hallam. *MetaPathways v2.0: A master-worker model for environmental Pathway/Genome Database construction on grids and clouds.* Proceedings of the 2014 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2014), Honolulu, HI, USA, May 21-24, 2014. [doi:10.1109/CIBCB.2014.6845516](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6845516)

