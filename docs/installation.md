# Installation and first run

Install MetaPathways, try the included three-sample dataset, then use your own data. Linux x86-64 is supported. Choose one installation method below; **Mamba is recommended**.

## 1. Conda package with Mamba (preferred)

```bash
mamba create -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.0
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

## 2. Quay: Docker or Apptainer

```bash
# Docker
docker pull quay.io/hallamlab/metapathways:4.0.0

# Or Apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:4.0.0
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

## Try your own data

Build a production reference database, then annotate an assembly:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways run -i /path/to/assembly.fasta -o results -d ~/MPDB --threads 8
```

For assemblies with reads and genome maps, follow the [complete workflow](inputs.md). **If you want PGDBs, complete the [Pathway Tools installation guide](pathway-tools.md) before starting that workflow.** The small test references are for testing only.

```{include} includes/cami-references.md
```
