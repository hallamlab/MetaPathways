# MetaPathways

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of information into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. MetaPathways v1.0 is a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

MetaPathways has since advanced as a modular tool, deepening our understanding of microbial metabolism at various biological levels. With this release, we have addressed previous challenges in modularity and database management. v3.5 enhances user accessibility through streamlined installation via package indexes or containers, refined modules, and interface upgrades. It boasts updated algorithm support for sequence feature prediction, annotation, metabolic inference, and coverage metrics. Tested on mock community data, Metapathways v3.5 demonstrates improved performance and usability. With automated installation and database management, this open-source tool makes advanced metagenomic analysis more accessible. Metapathways v3.5 represents a significant step forward in automated, comprehensive metagenomic analysis, facilitating a deeper exploration of microbial interactions and metabolic functions in environmental genomics.

**[Full user guide](https://hallamlab-metapathways.readthedocs.io/en/latest/index.html)** · [Workflow test](https://hallamlab-metapathways.readthedocs.io/en/latest/test.html) · [Issues and feature requests](https://github.com/hallamlab/MetaPathways/issues)

## Quick start

On Linux x86-64, install MetaPathways and its workflow dependencies with Mamba:

```bash
mamba create -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.1
conda activate metapathways
```

Prefer a container or a source checkout? Follow the [Docker / Apptainer guide](https://hallamlab-metapathways.readthedocs.io/en/latest/containers.html) or [GitHub installation guide](https://hallamlab-metapathways.readthedocs.io/en/latest/installation.html).

## Try the bundled three-sample dataset

```bash
metapathways prepare_test -o ~/mp-test
cd ~/mp-test
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-test/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
metapathways report -o all --serve --no-browser --port 8765
```

Open the report server's printed URL; for remote runs, follow the [SSH viewing guide](https://hallamlab-metapathways.readthedocs.io/en/latest/reports-tutorial.html#view-a-remote-report-through-ssh). The bundled 2.4 MiB CAMI II subset ([Meyer et al., 2022](#cami-references)) exercises annotation, read abundance, genome splitting, reports and exploration. Database preparation downloads supporting reference records. Pathway inference requires your own license and is skipped in this test. See the [test walkthrough](https://hallamlab-metapathways.readthedocs.io/en/latest/test.html) for validation and expected results.

## Run your own data

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways run -i /path/to/assembly.fasta -o results -d ~/MPDB --threads 8
```

**For PGDBs, follow the [Pathway Tools installation guide](https://hallamlab-metapathways.readthedocs.io/en/latest/pathway-tools.html) first.** Then follow the [complete workflow guide](https://hallamlab-metapathways.readthedocs.io/en/latest/analysis.html) to include reads and genome maps. Test references are for testing only.

## Workflow

[![MetaPathways appnote-style workflow: sequence processing, annotation, optional pathways and abundance, reports and explorer.](docs/assets/workflow-main.svg?v=mp-tools-20261008)](https://hallamlab-metapathways.readthedocs.io/en/latest/workflow.html)

The **[full user guide](https://hallamlab-metapathways.readthedocs.io/)** covers inputs, databases, Pathway Tools, local and Slurm resources, all commands, reporting and troubleshooting. [Detailed workflow and tool citations](https://hallamlab-metapathways.readthedocs.io/en/latest/detailed-workflow.html).

## Team, support and citation

Current team: Ryan J. McLaughlin, Tony X. Liu, Tomer Altman, Aditi N. Nallan, Aria S. Hahn, Julia Anstett, Connor Morgan-Lang, Kishori M. Konwar and Steven J. Hallam. Previous contributors include Niels W. Hanson and Shang-Ju Wu.

Source: [hallamlab/MetaPathways](https://github.com/hallamlab/MetaPathways). Historical code: [MetaPathways-legacy](https://github.com/hallamlab/MetaPathways-legacy). Questions and bug reports: [GitHub issues](https://github.com/hallamlab/MetaPathways/issues). License: [MIT](LICENSE), with bundled third-party license notices retained.

Please cite:

> McLaughlin RJ, Liu TX, Altman T, Nallan AN, Hahn AS, Anstett J, Morgan-Lang C, Konwar KM, Hallam SJ. *MetaPathways v3.5: Modularity and Scalability Improvements for Pathway Inference from Environmental Genomes*. bioRxiv (2024). [doi:10.1101/2024.06.04.597460](https://doi.org/10.1101/2024.06.04.597460).

## CAMI references

- **CAMI:** Sczyrba, A., Hofmann, P., Belmann, P., et al. (2017). *Critical Assessment of Metagenome Interpretation—a benchmark of metagenomics software*. **Nature Methods 14**(11), 1063–1071. [DOI: 10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458). [CAMI project website](https://cami-challenge.org/).
- **CAMI II:** Meyer, F., Fritz, A., Deng, Z.-L., et al. (2022). *Critical Assessment of Metagenome Interpretation: the second round of challenges*. **Nature Methods 19**(4), 429–440. [DOI: 10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4). [CAMI project website](https://cami-challenge.org/).
- **Source dataset for the MP test subset:** CAMI II multi-sample human microbiome dataset. [Dataset DOI: 10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). The bundled inputs are selected and cropped subsets of this collection; their exact transformations and file hashes are recorded in the bundle provenance.
