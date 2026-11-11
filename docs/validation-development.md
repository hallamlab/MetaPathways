# Validation and test-fixture maintenance

[Development](development.md) · [PR testing](pr-testing.md)

For contributors validating software changes or regenerating bundled fixtures.

## Installation example versus biological validation

`prepare_test -o WORKSPACE` copies the bundled three-sample inputs and reference FASTAs. From that workspace, `build_db --test -d MPDB` formats the fixtures and downloads current ExPASy/NCBI support records. The workflow test is neither offline nor fully pinned; record the downloaded reference dates. Its default workflow skips Pathway Tools.

The test suite includes synthetic/unit checks for read-layout arguments, task ordering/resources, resume/cleanup, container isolation interfaces, report joins and export behavior. Real synthetic Nextflow tasks exercise scheduling and logging; concurrent Pathway Tools startup checks exercise private container state. These checks do not substitute for real annotations, PGDB construction, or execution on an actual Slurm cluster. The small test walkthrough and any user-run benchmarks must be evaluated through their retained results and task records. Do not infer full biological validation from unit-test success. See the [test protocol](test.md) and [benchmark measurement guide](benchmarking.md).

Run the software checks from the checkout:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python -m unittest discover -s tests/release -p 'test_*.py'
python scripts/check_docs.py
```

The report HTTP tests require permission to bind a temporary loopback socket. Biological integration checks require the full Conda environment, appropriate references and an explicit decision to run them. Reports created from existing outputs perform indexing only.

## Regenerate the CAMI test fixture

Maintainers can reproduce the data from an MP manifest pointing at the original local CAMI downloads:

```bash
python scripts/prepare_cami_test.py \
  --source-manifest /path/to/full-cami-manifest.tsv \
  --output /path/to/new-cami-test \
  --minimap2 /path/to/minimap2
```

Run from the source checkout. The output directory must not already exist. Source assemblies must have a sibling `gsa_mapping.tsv`, and the selected sample read inputs must be the original interleaved FASTQs. Python and minimap2 are required. Temporary read subsets and alignments are cleaned up after preparation. Reproduction requires the original large CAMI files; tests use the bundled small files directly.

## Historical benchmark provenance

The repository's previous reproducibility notes describe the [v1 preprint](https://doi.org/10.1101/2024.06.04.597460), section 3.1, as reporting 15 CAMI II human-microbiome metagenomes ([Meyer et al., 2022](#cami-references)) and 622 MAGs from five body sites, SwissProt 2023_05, MetaCyc 27.1, MetaBAT2 2.15 and Pathway Tools 27.0, with 16 cores on Ubuntu 20.04.6. These are historical settings, not the defaults or a newly verified reproduction of this source revision.

The [CAMI2 study](https://doi.org/10.1038/s41592-022-01431-4) identifies the human collection at [PUBLISSO, DOI 10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). A collection DOI does not identify the exact assembly/read/bin selections used by a particular benchmark. The K12 installation fixture does not substitute for those records.

A fresh benchmark needs a complete sample/input manifest, contig maps, reference/software versions, resource limits, commands and retained task traces. Historical output affected by the paired-read mapping bug must have mapping corrected before abundance is used. Expected MAG Pathway Tools failures should be recorded as such; they are not evidence that no pathways exist.

Only measured trace fields should support resource plots. A report rebuilt over old outputs cannot recreate missing peak RAM, CPU time or stage runtime. Reused tasks, failed tasks, incomplete output and successful fresh tasks must be distinguished when making supplementary tables or benchmark claims. Report database row counts are data-accounting units, not performance measurements.

## Recorded three-sample validation

The 2026-10-02 local SwissProt test run completed all 51 tasks successfully. The [audit record](validation/reviewer-2026-10-02.json) confirms all 393 CDS records were retained, paired read identifiers matched, all 421 gene/RNA abundance records and 9 contig records matched their source measurements, and taxonomy remained tied to each reference database and target. Each sample produced three CAMI genome bins. The 28 report placeholders were RNA features, not missing CDS annotations. Skin_28 had no rRNA queries, explaining its two BLAST empty-query warnings. This run used `--skip_ptools`; it does not establish successful PGDB construction or validate full-scale/HPC performance.

```{include} includes/cami-references.md
```

## Historical software validation notes

The environment specification in `docker/conda_base.yml` is not an exact lock. MAGSplitter and Camelot revisions are pinned in `requirements-workflow.txt`; the release package includes both helpers from those exact sources. Record the MP artifact checksum along with the environment. The source revision's setup checks used Nextflow 26.04.6, Python 3.11 and Apptainer 1.5.4, with Pathway Tools 29.5 startup validation. This does not establish end-to-end biological equivalence for the orchestration migration.

## Documentation and release records

The Markdown guides and generated CLI references are checked in CI and built for
Read the Docs. Git history retains the former RST documentation. Consult the
[release process](releasing.md) for publication and archival DOI requirements;
a source tag alone does not establish package, container or archive publication.
