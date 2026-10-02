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
python -m pip install --no-deps --no-build-isolation -e .
python -m pip install \
  'git+https://github.com/hallamlab/MAGSplitter.git@main' \
  'git+https://bitbucket.org/tomeraltman/camelot-frs@dev#egg=camelot-frs'
```

This is the development branch documented here. For a release, use its tag and matching documentation. The environment installs the workflow and biological dependencies; MAGSplitter handles genome-bin splitting and Camelot handles pathway extraction. The installation can require substantial downloads. Do not use `sudo pip` to install MP into the operating system's Python.

`-e .` installs MP in editable mode: Python modules are read from this checkout. Some helper scripts are copied into the environment, so rerun the MP `pip install` command after updating the checkout. Package version `3.5.1` alone does not distinguish this feature branch from older code.

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

```bash
mkdir -p ~/mp-reviewer
cd ~/mp-reviewer
metapathways build_db --test
metapathways run --test
```

The database command formats the bundled small reference fixtures and downloads support records for enzymes and taxonomy; it needs internet access and a writable installed fixture directory. The run command uses the bundled K12 assembly and paired reads, forces one thread for this example, and writes `~/mp-reviewer/test/k12_test/`. It does not invoke Pathway Tools. Do not pass your production data to `--test`: that flag selects its own inputs and output location.

Expect annotation files under `test/k12_test/results/annotation_table/`, abundance under `test/k12_test/results/rpkm/`, and `test/reports/MP_run_report.html`. Inspect task status and logs if any expected output is absent. A report file existing is not by itself proof that all stages succeeded.

```bash
metapathways report -o test --serve --no-rebuild
```

Keep this terminal running while browsing. Ctrl-C stops the report server, not the already completed analysis. On a remote server, use `--no-browser --port 8765` and the SSH tunnel instructions.

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
