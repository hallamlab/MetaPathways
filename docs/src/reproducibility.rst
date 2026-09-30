Reproducibility
***************

Source and license
==================

The canonical v3.5 source is `hallamlab/MetaPathways
<https://github.com/hallamlab/MetaPathways>`_. MetaPathways is distributed under
the MIT license; bundled third-party source retains its own license notices.
The historical GitHub codebase is separate from v3.5.

The Conda package and Quay container are distributed separately. A GitHub tag
alone does not update either distribution. Record the package version or
container digest used for an analysis, alongside the Git commit when installing
from source.

Installation test
=================

In a writable environment with the dependencies from ``docker/conda_base.yml``
and this checkout installed, run:

.. code-block:: bash

   metapathways build_db --test
   metapathways run --test

The package includes the compressed K12 nucleotide FASTA, paired FASTQ reads,
and small SwissProt/SILVA sequence fixtures in ``metapathways/regtests``.
The first command formats these fixtures and downloads the current ExPASy enzyme
records and NCBI taxonomy. It needs internet access and writes into the installed
package's test database directory. It is not an offline or fully pinned benchmark.
The second command writes results to ``./test/``. Inspect command output,
``metapathways_steps_log.txt`` and ``errors_warnings_log.txt`` in the sample output
directory; do not rely on the process exit code alone.

Keep generated databases and outputs outside version control. Preserve their
checksums and download dates when recording an analysis. The raw sequence fixtures
are versioned with the source; they must not be removed as generated data.

Supported inputs and outputs
============================

The current ``run`` CLI accepts nucleotide FASTA (``--input_format fasta``) or
amino-acid FASTA (``--input_format fasta-amino``). Compressed FASTA is used by the
included test. FASTQ reads supplied with ``-1``/``-2`` or ``--interleaved`` enable
abundance calculations. GFF and GenBank are outputs, not supported primary input
formats of this CLI.

Outputs under each sample directory include:

* ``preprocessed/`` and ``orf_prediction/``: processed nucleotide/protein sequences
  and predicted features.
* ``genbank/``: annotated GFF and GenBank files when their stages are enabled.
* ``results/annotation_table/``: tab-separated functional/taxonomic annotations.
* ``results/rpkm/``: abundance/count tables when read mapping is enabled.
* ``results/rRNA/`` and ``results/tRNA/``: RNA annotations and statistics.
* ``run_statistics/`` and the sample log files: processing counts and step status.
* ``ptools/``: Pathway Tools input files. Creating an ePGDB is a separate
  ``metapathways ptools`` command requiring a licensed Pathway Tools installation
  and the optional PGDB dependencies.

Dependency versions
===================

``docker/conda_base.yml`` specifies the runtime environment; optional development,
MAG splitting and PGDB dependencies are listed in ``docker/conda_dev.yml`` and
``conda_recipe/meta_template.yaml``. These environment specifications are not
complete version locks. The migration validation environment on Linux x86-64 used:

.. list-table:: Validation environment
   :header-rows: 1

   * - Dependency
     - Version
   * - Python
     - 3.10.19
   * - Snakemake
     - 7.32.4
   * - Prodigal / PProdigal
     - 2.6.3 / 1.0.1
   * - BLAST+
     - 2.17.0
   * - BWA / SAMtools
     - 0.7.19 / 1.22.1
   * - CoverM
     - 0.7.0
   * - tRNAscan-SE / Barrnap
     - 2.0.12 / 0.9
   * - pandas / pybedtools / pysam / pyfastx
     - 2.3.3 / 0.12.0 / 0.23.3 / 2.2.0

FAST and metacount executables are included in the package. Their native source,
Makefiles, and relevant third-party tests are retained under ``extensions/``.
The bundled executables target Linux; native Windows is not supported. Record an
explicit environment export and checksums of reference files for reproducible
analyses.

Manuscript benchmark provenance
===============================

The `v1 preprint <https://doi.org/10.1101/2024.06.04.597460>`_, section 3.1,
reports 15 CAMI2 Human Microbiome metagenomes and 622 MAGs across five body sites.
It specifies SwissProt release 2023_05, MetaCyc 27.1, MetaBAT2 2.15, and Pathway
Tools 27.0, with 16 cores on Ubuntu 20.04.6. These historical reference/software
versions differ from the moving downloads used by the installation test.

The `CAMI2 study <https://doi.org/10.1038/s41592-022-01431-4>`_ identifies the
human dataset at `PUBLISSO, DOI 10.4126/FRL01-006425518
<https://doi.org/10.4126/FRL01-006425518>`_. This identifies the source collection;
it does not establish the exact input files used for the MetaPathways results.

Before claiming reproduction of the manuscript benchmark, the authors still need
to provide the exact 15-sample input/assembly manifest with checksums, commands
for preparing the 622 MAGs and contig mappings, reference acquisition/version
records, the pipeline revision and environment, benchmark commands, and the
result tables underlying the reported metrics. These artifacts were not found in
the current source checkout. The bundled K12 example does not substitute for
this benchmark archive.

No Zenodo DOI has been assigned by this migration. After the validated GitHub
release is published and archived through Zenodo, cite the actual version DOI
in the manuscript and data-availability statement.

Automated checks
================

The GitHub Actions smoke workflow builds the source/wheel distributions, checks
installed CLI startup on Python 3.10, and builds documentation. The full
bioinformatics example above is a separate integration check requiring the
Conda dependencies and public database downloads. It is not run by the lightweight
smoke workflow. The old Python 3.6--3.8 CircleCI/tox configuration was retired
because it referenced a removed test suite and requirements file.
