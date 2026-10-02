# Getting started: from a terminal to your first result

[Home and reading order](../README.md#start-here) · Next: [reviewer walkthrough](reviewer-test.md)

## What you will do

Install MP on a Linux computer, activate its software environment, run a small included example, and open a report. You do not need a Pathway Tools license for the first example. Add pathway inference after the example works.

If you already have a working MP installation, go directly to the [reviewer walkthrough](reviewer-test.md). If your lab manages the software centrally, ask for the activation command and MPDB path rather than installing a second copy.

## Before typing commands

A **terminal** accepts commands. A **working directory** is the folder relative paths refer to. `pwd` displays it; `ls` lists its contents; `cd DIRECTORY` changes it. `~` means your home directory. `/data/project/sample.fasta` is an absolute path; `sample.fasta` means a file in the working directory. Linux distinguishes `SampleA` from `samplea`.

Copy commands inside the code blocks, without Markdown backticks or a shell prompt. A backslash at the end of a line continues the command on the next line: do not put spaces after it. Values such as `/path/to/MPDB` are placeholders; replace them with real paths. Use simple project and sample names without spaces or shell punctuation, because some legacy tool wrappers construct shell commands.

An **environment** is a collection of compatible software. Activating it makes your terminal find the intended MP and tool executables. Activate it again in each new terminal. It does not require setting an MP-specific environment variable.

On a remote server, the analysis commands run in the SSH terminal on that server. Your web browser normally runs on your own computer. The [SSH report instructions](reports-tutorial.md#view-a-remote-report-through-ssh) connect the two. For a long run, use your site's persistent terminal/session practice so disconnecting SSH does not interrupt your controller.

## Install on Linux x86-64

The supported source-install path uses Linux x86-64, Conda/Mamba, and Git. Windows and macOS users can connect to a Linux server. This guide does not claim a tested native Windows, Apple Silicon, or WSL installation.

If `conda --version` and `mamba --version` do not work, follow the official [Miniforge installation instructions](https://github.com/conda-forge/miniforge). Reopen your terminal after shell initialization. Ask your administrator about Git and cluster installation policies if needed.

Use a writable source checkout. The following commands put it in `~/src/MetaPathways`, which the reviewer walkthrough also uses:

```bash
mkdir -p ~/src
cd ~/src
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
git checkout feat/nextflow-controller-db-build
mamba env create -f docker/conda_base.yml
conda activate metapathways
mamba install --yes --override-channels --strict-channel-priority -c conda-forge -c bioconda pip wheel git
python -m pip install --no-deps --no-build-isolation .
python -m pip install --no-deps -r requirements-workflow.txt
```

This is the development branch documented here. For a release, use its tag and matching documentation. The environment installs the workflow and biological dependencies; MAGSplitter handles genome-bin splitting and Camelot handles pathway extraction. The installation can require substantial downloads. Do not use `sudo pip` to install MP into the operating system's Python.

This installs a fixed copy of MP into the environment. Rerun the MP `pip install` command after updating the checkout. Record `git rev-parse HEAD` alongside `metapathways version`: candidate revisions can share version `3.5.2`. Developers may opt into editable mode with `-e .`.

Confirm the installation:

```bash
which python
which metapathways
metapathways version
metapathways analysis_wf --help
nextflow -version
java -version
apptainer --version
magsplitter --help
```

`which` should point into your activated environment, or to an intentional site-managed executable. `analysis_wf --help` should list manifest and resource options. A help command does not run an analysis. If a command is missing, first confirm you activated the correct environment. On hosts that restrict user namespaces, Apptainer setup may require an administrator; see [Pathway Tools build prerequisites](pathway-tools.md#host-prerequisites).

## First installation check

Follow the **[three-sample reviewer walkthrough](reviewer-test.md)**. It uses the same `analysis_wf` command as a real analysis, with small bundled references, assemblies, paired reads and genome maps. It exercises annotation, read abundance, genome splitting and the reports without requiring Pathway Tools. Database preparation downloads enzyme and taxonomy support records, so internet access is required.

The older `run --test` K12 example remains available for compatibility, but it does not exercise the complete workflow and is not the reviewer acceptance test.

## Your first real assembly

The tiny test references are only for the example. Build or obtain a production MPDB before interpreting your own data:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways run -i /path/to/sample.fasta -o ~/mp-results -d ~/MPDB \
  --threads 8 --max_cpus 8 --max_memory '32 GB'
```

The assembly goes through quality control, gene/RNA prediction, reference searches, and annotation. Without reads, there is no measured read abundance. Without a separate `ptools` step, there are no inferred PGDB pathways. A single `analysis_wf` command can combine these steps once the inputs and licensed container are ready.

Continue with the [reviewer walkthrough](reviewer-test.md), [command cookbook](commands.md), and [input organization](inputs.md).

## Terms used in the guides

| Term | Meaning in MP |
| --- | --- |
| Assembly / contig | DNA sequences assembled from reads; each FASTA record is a contig |
| FASTA / FASTQ | Sequence file / read file with per-base quality scores |
| Paired / interleaved | Two mate files / one file with alternating mate records |
| ORF | Predicted protein-coding region; RNA features are also retained in relevant outputs |
| MAG / genome bin | A group of contigs assigned to a genome; MP consumes assignments, not a binning algorithm |
| MPDB | MP's reference sequences, indexes, taxonomy and functional mapping tables |
| MetaCyc | Reference pathways/reactions and associated proteins; not a sample's pathway predictions |
| PGDB | Pathway/Genome Database inferred for one community or genome bin |
| SIF | An Apptainer image containing the licensed Pathway Tools installation |
| Nextflow / task | Scheduler behind MP / one scheduled unit of work |
| Manifest | A TSV listing sample IDs and exact input paths |
| Receipt | MP's record used to decide whether a completed task can be reused |
| TSV / CSV | Table with tab-separated / comma-separated fields |

## Wheel-build cleanup errors on shared filesystems

If pip fails with `Directory not empty` while removing `build/bdist.../wheel`, the installation did not complete. This is a build-directory cleanup failure, not an analysis or Slurm error. Update the checkout, then build from a clean snapshot in node-local `/tmp` using your activated MP environment:

```bash
git pull --ff-only
mp_build_dir=$(mktemp -d /tmp/metapathways-build.XXXXXX)
git archive HEAD | tar -x -C "$mp_build_dir"
python -m pip install --no-cache-dir --no-deps --no-build-isolation "$mp_build_dir"
```

Only continue to the workflow after pip reports successful installation. The snapshot uses committed files from the checkout; it excludes stale build products and uncommitted edits. Your existing MPDB, SIF, inputs and outputs stay in place. The temporary build directory can be removed after a successful installation:

```bash
rm -rf -- "$mp_build_dir"
```
