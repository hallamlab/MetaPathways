# Quay containers: Docker and Apptainer

[Home and quick start](installation.md) · [Reviewer dataset](reviewer-test.md) · [Licensed Pathway Tools](pathway-tools.md)

Use the versioned `quay.io/hallamlab/metapathways:3.5.2` image on Linux x86-64. It includes MP, its workflow dependencies (including MAGSplitter and Camelot), and the three-sample reviewer dataset. Production references and licensed Pathway Tools are supplied separately.

## Docker three-sample test

Create a working directory and open a shell in the image. The mount keeps inputs, references and results on your host. Run these commands in your host terminal:

```bash
mkdir -p ~/mp-reviewer-docker
cd ~/mp-reviewer-docker
docker pull quay.io/hallamlab/metapathways:3.5.2
docker run --rm -it --network host --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work quay.io/hallamlab/metapathways:3.5.2 bash
```

Inside the container, run:

```bash
metapathways prepare_test -o .
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
metapathways report -o all --serve --no-browser --port 8765
```

Open the printed URL in your browser. On a remote Linux host, use the [SSH tunnel instructions](reports-tutorial.md#view-a-remote-report-through-ssh). Host networking makes the server's loopback address available on that host. Keep the server running while browsing; Ctrl-C stops it. Type `exit` to close the container shell. Files in the working directory remain on the host and belong to your user.

To return later, repeat the `docker run` command from the same host directory, then run `metapathways report -o all --serve --no-browser --port 8765`.

## Apptainer three-sample test

Create a working directory and pull the image on an internet-connected host:

```bash
mkdir -p ~/mp-reviewer-apptainer
cd ~/mp-reviewer-apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:3.5.2
apptainer exec --bind "$PWD:/work" --pwd /work metapathways.sif bash
```

Inside the container, run:

```bash
metapathways prepare_test -o .
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
metapathways report -o all --serve --no-browser --port 8765
```

Open the printed URL, using an [SSH tunnel](reports-tutorial.md#view-a-remote-report-through-ssh) if remote. Ctrl-C stops the report server; `exit` closes the shell. Your inputs, references and results remain in the host working directory. The SIF stays read-only. Pulling this MP image does not require building a custom image or using `--fakeroot`.

## From the test to your own data

The reference build downloads enzyme and taxonomy support files. Complete reference preparation on an internet-connected host before using offline compute nodes. The tiny reviewer references are for testing only. Build production references with `metapathways build_db -d MPDB --func swissprot -a fast` in a separate project directory.

Mount all inputs, references and outputs into the container and use their container-visible paths in manifests. For Slurm, use the [Mamba installation](installation.md#1-conda-package-with-mamba-preferred) on shared storage so the controller and compute jobs can use the same environment. Follow the [resource and Slurm guide](resources.md#resources-and-slurm) for submission limits.

The public MP image does not include Pathway Tools or MetaCyc. Licensed users should follow the [Pathway Tools guide](pathway-tools.md) to build a separate SIF from their own installer. That image can be copied to another compatible host and selected with `--image`.
