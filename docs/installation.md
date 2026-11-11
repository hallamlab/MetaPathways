# Installation and first run

For MetaCyc references or pathway inference, [build/register the licensed SIF](pathway-tools.md#build-and-register-once)
first. Then use [database setup](databases.md) to select your references.
The bundled installation test skips this licensed setup.

Install MetaPathways, try the included three-sample dataset, then use your own data. Linux x86-64 is supported. Choose one installation method below; **Mamba is recommended**.

You need a Linux x86-64 host and an internet connection. If Mamba is not installed, follow the [Miniforge installation instructions](https://github.com/conda-forge/miniforge) or your site’s environment setup. Choose one installation route; the other routes below are alternatives.

## 1. Conda package with Mamba (preferred)

```bash
mamba create -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.8
conda activate metapathways
```

This installs MP and its workflow dependencies, including MAGSplitter and Camelot. Prepare the included data, build the small reference database, and run all three samples:

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

Open the URL printed by the report server. On a remote server, use an [SSH tunnel](reports-tutorial.md#view-a-remote-report-through-ssh). The **2.4 MiB input dataset is included** in the package and is derived from CAMI II ([Meyer et al., 2022](#cami-references)); database preparation downloads enzyme and taxonomy support records. The test covers annotation, paired-read abundance, genome splitting, reports and exploration. Pathway inference is skipped because it requires your own Pathway Tools license.

`prepare_test` creates `cami-test/` and `MPDB/` inside the requested workspace.
Run the subsequent commands from that workspace (`cd ~/mp-test` above), so
`-d MPDB` selects its reference database. `build_db --test` also copies missing
bundled reference seeds into its destination, preserving existing files.

### Before adding licensed Pathway Tools

If you want MetaCyc or pathway inference, check container readiness in the active
Mamba environment before obtaining/running the licensed installer:

```bash
metapathways build_pt --check -o ~/mp-containers
```

This downloads nothing and needs no sudo. MP automatically tests the Mamba
Apptainer and available system/site installations, then remembers a working choice
for this host. On a compatible host, MP handles the
public container dependencies automatically. If host policy blocks the check,
the log identifies what your administrator needs to investigate; Mamba cannot
grant host permissions. Follow the [Pathway Tools setup](pathway-tools.md#host-prerequisites)
before continuing with licensed features. The public test above is independent.

## 2. Quay: Docker or Apptainer

```bash
# Docker
docker pull quay.io/hallamlab/metapathways:4.0.8

# Or Apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:4.0.8
```

The image includes the same workflow dependencies and test data. Follow the [Docker three-sample test](containers.md#docker-three-sample-test) or [Apptainer three-sample test](containers.md#apptainer-three-sample-test) to run the commands with your working directory mounted for persistent results. Licensed Pathway Tools is a separate image.

## 3. Local installation from GitHub

```bash
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
mamba env create -f docker/conda_base.yml
mamba run -n metapathways python -m pip install .
conda activate metapathways
```

Then run the **same three-sample commands under option 1**, starting with `metapathways prepare_test -o ~/mp-test`. MP installs its Python workflow helpers automatically. The data comes from the installed package; the test does not depend on your checkout location.

## Upgrade to 4.0.8

This update automatically selects a working Mamba or site Apptainer, checks
building and execution separately, and supports older build options. The cached
public dependency image avoids APT during licensed builds. Run `build_pt --check`
after updating if you plan to use licensed features. Existing working SIFs do
not need rebuilding for this change.

This release separates image and database setup: `build_pt` builds only the SIF;
use `build_db` for MetaCyc and public references. Existing MPDB additions require
`--add`, and replacing selected references requires `--update`. See the
[database guide](databases.md) before modifying an existing MPDB.

For an existing Conda installation, update the environment before submitting jobs:

```bash
mamba install -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.8
conda activate metapathways
metapathways version
```

For containers, pull the new versioned image above and use it for subsequent runs.
Updating a Git checkout does not update an installed Conda package or SIF.
Stop active workflows before replacing their software environment.

The FAST staging and search-walltime options are opt-in; existing commands retain
their defaults. See [Slurm guidance](resources.md#fast-databases-on-shared-storage).
Previously compacted outputs use implementation-specific completion markers and
require a new output directory when upgrading. Noncompacted runs may recompute
tasks whose tracked command or input/output evidence differs. Externally generated
FAST results are not automatically accepted as workflow successes.

After the test finishes, use the [test walkthrough](test.md#3-check-outputs-and-explore-tables) to check its outputs. Stop the report server with Ctrl-C before entering more commands in the same terminal.

## Try your own data

Build a production reference database, then annotate an assembly:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways run -i /path/to/assembly.fasta -o results -d ~/MPDB --threads 8
```

For assemblies with reads and genome maps, follow the [complete workflow](inputs.md). **If you want PGDBs, complete the [Pathway Tools installation guide](pathway-tools.md) before starting that workflow.** The small test references are for testing only.

```{include} includes/cami-references.md
```
