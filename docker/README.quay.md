# Quay containers: Docker and Apptainer

[Home and quick start](../README.md#quick-start) · [Reviewer dataset](../docs/reviewer-test.md) · [Licensed Pathway Tools](../docs/pathway-tools.md)

Use a versioned `quay.io/hallamlab/metapathways` image on Linux x86-64. These examples use MP **3.5.2**.

The public MP image contains annotation software and the small reviewer fixtures, not production reference databases or licensed Pathway Tools. Current core release images omit MAGSplitter and Camelot; the one-time helper layer below supplies their pinned revisions. It does not rebuild MP or include any licensed material. The same three-sample test is used for both container runtimes.

## Docker three-sample test

Docker must already work on the host. Create a new directory; all reference downloads and results will stay here.

```bash
mkdir -p ~/mp-reviewer-docker
cd ~/mp-reviewer-docker
docker pull quay.io/hallamlab/metapathways:3.5.2

cat > Dockerfile <<'DOCKERFILE'
FROM quay.io/hallamlab/metapathways:3.5.2
USER root
RUN micromamba install --yes --prefix /opt/conda --override-channels -c conda-forge pip git && \
    /opt/conda/bin/python -m pip install --no-deps \
      'magsplitter @ git+https://github.com/hallamlab/MAGSplitter.git@a83bcd0c6479d3fda15025262b981e348c082bc9' \
      'camelot-frs @ git+https://bitbucket.org/tomeraltman/camelot-frs.git@30a774c7fe7bb8a88b5cd6cadfac1c6d258efc96' && \
    micromamba clean --all --yes
DOCKERFILE

docker build -t metapathways-workflow:3.5.2 .
```

Copy the included inputs and reference seeds out of the image. Determine its installed test-database path so the writable copy can be mounted there; `build_db --test` uses that installed path.

```bash
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work metapathways-workflow:3.5.2 \
  python -c 'from pathlib import Path; import shutil, metapathways; p=Path(metapathways.__file__).parent / "regtests"; shutil.copytree(p / "cami_reviewer", "cami-reviewer"); shutil.copytree(p / "test_db", "MPDB")'

mp_test_db=$(docker run --rm metapathways-workflow:3.5.2 \
  python -c 'from pathlib import Path; import metapathways; print(Path(metapathways.__file__).parent / "regtests/test_db")')

mp_container() {
  docker run --rm --network host --env NXF_HOME=/work/.nextflow --user "$(id -u):$(id -g)" \
    -v "$PWD:/work" -v "$PWD/MPDB:$mp_test_db" -w /work \
    metapathways-workflow:3.5.2 "$@"
}
```

`mp_container` is a shell convenience for running a command inside this image. Keep using this terminal and working directory. Linux host networking lets the report server bind the host's loopback address, so the same local/SSH browser instructions apply. Output ownership follows your user ID. The wrapper puts Nextflow’s launcher cache in the writable work directory.

```bash
mp_container metapathways build_db --test
mp_container metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
mp_container metapathways report -o all --serve --no-browser --port 8765
```

Open the printed report URL. Keep the terminal running while browsing; Ctrl-C stops only the report server after the analysis has finished. Reports and data remain in the mounted host directory when a container exits.

## Apptainer three-sample test

Use a host where Apptainer can pull images and build with `--fakeroot`. Some HPC sites require an administrator-provided build service or require building images elsewhere; a finished SIF can then be copied to shared storage. Prepare the image and reference support files on an internet-connected host, not an offline compute node.

```bash
mkdir -p ~/mp-reviewer-apptainer
cd ~/mp-reviewer-apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:3.5.2

cat > workflow.def <<'DEFINITION'
Bootstrap: localimage
From: metapathways.sif

%post
    /bin/micromamba install --yes --prefix /opt/conda --override-channels -c conda-forge pip git
    /opt/conda/bin/python -m pip install --no-deps \
      'magsplitter @ git+https://github.com/hallamlab/MAGSplitter.git@a83bcd0c6479d3fda15025262b981e348c082bc9' \
      'camelot-frs @ git+https://bitbucket.org/tomeraltman/camelot-frs.git@30a774c7fe7bb8a88b5cd6cadfac1c6d258efc96'
    /bin/micromamba clean --all --yes
DEFINITION

apptainer build --fakeroot metapathways-workflow.sif workflow.def
```

Copy the bundled inputs and test references to writable storage, then mount the reference copy over the image's read-only test-database directory:

```bash
apptainer exec --bind "$PWD:/work" --pwd /work metapathways-workflow.sif \
  python -c 'from pathlib import Path; import shutil, metapathways; p=Path(metapathways.__file__).parent / "regtests"; shutil.copytree(p / "cami_reviewer", "cami-reviewer"); shutil.copytree(p / "test_db", "MPDB")'

mp_test_db=$(apptainer exec metapathways-workflow.sif \
  python -c 'from pathlib import Path; import metapathways; print(Path(metapathways.__file__).parent / "regtests/test_db")')

mp_container() {
  apptainer exec --bind "$PWD:/work,$PWD/MPDB:$mp_test_db" --pwd /work \
    metapathways-workflow.sif "$@"
}

mp_container metapathways build_db --test
mp_container metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
mp_container metapathways report -o all --serve --no-browser --port 8765
```

Keep this terminal and working directory when using `mp_container`. The image remains read-only; database preparation and the workflow write through the bind mounts. For remote browsing, follow [SSH report access](../docs/reports-tutorial.md#view-a-remote-report-through-ssh).

## What the test covers

All three samples run annotation, paired-read abundance, genome splitting and report/explorer generation. `--skip_ptools` deliberately excludes licensed pathway inference. `all/reports/` contains the HTML report, EDA portal and indexed tables. Review task statuses as described in the [reviewer walkthrough](../docs/reviewer-test.md). Reference support downloads are substantially larger than the roughly 2.4 MiB input bundle.

For production HPC workflows, the [Mamba installation](../README.md#1-conda-package-with-mamba-preferred) plus the separate licensed Pathway Tools SIF is the recommended route. Launching Slurm from inside the public MP image requires site-specific scheduler access; running a Pathway Tools SIF from inside another container requires nested-container support. The examples above test the local executor inside Docker/Apptainer and do not claim either setup is universally supported.
