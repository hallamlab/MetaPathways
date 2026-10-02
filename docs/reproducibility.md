# Reproducibility and validation

[User guide](index.md) · [Release process](releasing.md)

## Record the actual revision and inputs

MetaPathways is distributed under the [MIT license](https://github.com/hallamlab/MetaPathways/blob/HEAD/LICENSE); bundled third-party sources retain their notices. GitHub tags, Conda packages and container images are separate publications. A tag alone does not update the other distributions. Record the source Git commit, package build or immutable container digest used for each analysis.

Preserve:

- Assembly/read identifiers, layout (single, paired or interleaved), checksums and sample naming.
- Original contig-to-MAG maps, including the binning version/method and QC decisions.
- MP commands, workflow parameters, tool versions and environment exports.
- Reference releases, acquisition dates and file checksums; public downloads can change.
- The Pathway Tools installer/image checksum, version and taxonomic-pruning setting.
- Final outputs, run logs, Nextflow trace/report/timeline and report `sources`/`schema.json`.

Example environment records:

```bash
mamba list --explicit > conda-explicit.txt
metapathways version > metapathways-version.txt
```

For a source installation, also record `git rev-parse HEAD` from the checkout and `python -m pip freeze` from its activated environment.

The environment specification in `docker/conda_base.yml` is not an exact lock. MAGSplitter and Camelot revisions are pinned in `requirements-workflow.txt`; the release package includes both helpers from those exact sources. Record the MP artifact checksum along with the environment. The source revision's setup checks used Nextflow 26.04.6, Python 3.11 and Apptainer 1.5.4, with Pathway Tools 29.5 startup validation. This does not establish end-to-end biological equivalence for the orchestration migration.

## Installation example versus biological validation

`prepare_test -o WORKSPACE` copies the bundled three-sample inputs and reference FASTAs. From that workspace, `build_db --test -d MPDB` formats the fixtures and downloads current ExPASy/NCBI support records. The reviewer test is neither offline nor fully pinned; record the downloaded reference dates. Its default workflow skips Pathway Tools.

The test suite includes synthetic/unit checks for read-layout arguments, task ordering/resources, resume/cleanup, container isolation interfaces, report joins and export behavior. Real synthetic Nextflow tasks exercise scheduling and logging; concurrent Pathway Tools startup checks exercise private container state. These checks do not substitute for real annotations, PGDB construction, or execution on an actual Slurm cluster. The small reviewer walkthrough and any user-run benchmarks must be evaluated through their retained results and task records. Do not infer full biological validation from unit-test success. See the [reviewer protocol](reviewer-test.md) and [benchmark measurement guide](benchmarking.md).

Run the software checks from the checkout:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python -m unittest discover -s tests/release -p 'test_*.py'
python scripts/check_docs.py
```

The report HTTP tests require permission to bind a temporary loopback socket. Biological integration checks require the full Conda environment, appropriate references and an explicit decision to run them. Reports created from existing outputs perform indexing only.

## Historical benchmark provenance

The repository's previous reproducibility notes describe the [v1 preprint](https://doi.org/10.1101/2024.06.04.597460), section 3.1, as reporting 15 CAMI2 Human Microbiome metagenomes and 622 MAGs from five body sites, SwissProt 2023_05, MetaCyc 27.1, MetaBAT2 2.15 and Pathway Tools 27.0, with 16 cores on Ubuntu 20.04.6. These are historical settings, not the defaults or a newly verified reproduction of this source revision.

The [CAMI2 study](https://doi.org/10.1038/s41592-022-01431-4) identifies the human collection at [PUBLISSO, DOI 10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). A collection DOI does not identify the exact assembly/read/bin selections used by a particular benchmark. The K12 installation fixture does not substitute for those records.

A fresh benchmark needs a complete sample/input manifest, contig maps, reference/software versions, resource limits, commands and retained task traces. Historical output affected by the paired-read mapping bug must have mapping corrected before abundance is used. Expected MAG Pathway Tools failures should be recorded as such; they are not evidence that no pathways exist.

Only measured trace fields should support resource plots. A report rebuilt over old outputs cannot recreate missing peak RAM, CPU time or stage runtime. Reused tasks, failed tasks, incomplete output and successful fresh tasks must be distinguished when making supplementary tables or benchmark claims. Report database row counts are data-accounting units, not performance measurements.

## Documentation and release records

Documentation is maintained in this README-led GitHub tree. CI checks local Markdown links and generated CLI reference freshness; no Sphinx/Read the Docs deployment is required. Git history retains the former RST documentation.

No new version DOI or archive has been assigned by this feature work. Publish and archive a validated release before claiming a corresponding version DOI in a manuscript. Use [maintainer release instructions](releasing.md) for package/container publication.
