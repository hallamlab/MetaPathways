# Quay containers: Docker and Apptainer

[Home and quick start](installation.md) · [Test dataset](test.md) · [Licensed Pathway Tools](pathway-tools.md)

Use the versioned `quay.io/hallamlab/metapathways:4.0.1` image on Linux x86-64. It includes MP, its workflow dependencies (including MAGSplitter and Camelot), and the three-sample test dataset derived from CAMI II ([Meyer et al., 2022](#cami-references)). Production references and licensed Pathway Tools are supplied separately.

## Docker three-sample test

Create a working directory and open a shell in the image. The mount keeps inputs, references and results on your host. Run these commands in your host terminal:

```bash
mkdir -p ~/mp-test-docker
cd ~/mp-test-docker
docker pull quay.io/hallamlab/metapathways:4.0.1
docker run --rm -it --network host --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" -w /work quay.io/hallamlab/metapathways:4.0.1 bash
```

Inside the container, run:

```bash
metapathways prepare_test -o .
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-test/all.tsv -o all -d MPDB \
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
mkdir -p ~/mp-test-apptainer
cd ~/mp-test-apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:4.0.1
apptainer exec --bind "$PWD:/work" --pwd /work metapathways.sif bash
```

Inside the container, run:

```bash
metapathways prepare_test -o .
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-test/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
metapathways report -o all --serve --no-browser --port 8765
```

Open the printed URL, using an [SSH tunnel](reports-tutorial.md#view-a-remote-report-through-ssh) if remote. Ctrl-C stops the report server; `exit` closes the shell. Your inputs, references and results remain in the host working directory. The SIF stays read-only. Pulling this MP image does not require building a custom image or using `--fakeroot`.

## From the test to your own data

The reference build downloads enzyme and taxonomy support files. Complete reference preparation on an internet-connected host before using offline compute nodes. The tiny test references are for testing only. Build production references with `metapathways build_db -d MPDB --func swissprot -a fast` in a separate project directory.

Mount all inputs, references and outputs into the container and use their container-visible paths in manifests. For Slurm, use the [Mamba installation](installation.md#1-conda-package-with-mamba-preferred) on shared storage so the controller and compute jobs can use the same environment. Follow the [resource and Slurm guide](resources.md#resources-and-slurm) for submission limits.

The public MP image does not include Pathway Tools or MetaCyc. Licensed users should follow the [Pathway Tools guide](pathway-tools.md) to build a separate SIF from their own installer. That image can be copied to another compatible host and selected with `--image`.

```{include} includes/cami-references.md
```

## Nested Pathway Tools containers

MP's release container includes an unlicensed Ubuntu dependency base for
`build_pt`. When present, MP uses it automatically instead of installing system
packages inside the nested build. Conda installations without the bundled base
continue to build from Ubuntu using host Apptainer.

A local nested build passed installation, official patch startup, SIF creation,
BLAST validation and Pathway Tools readiness. Complete nested PGDB inference and
portability across HPC installations require separate validation. Use the Mamba
route when a host restricts nested user namespaces or container execution.

### What the nested checks established

- Java/Nextflow stalled on a FUSE filesystem request when executing from the
  outer SIF on the test host. Launching the outer image with `--unsquash`
  resolved this; extraction requires additional temporary space and startup time.
- Outer bind mounts were inherited by the inner build and targeted directories
  absent from its writable root. Clearing inherited bind variables inside the
  outer container let the build proceed. Explicit inner runtime mounts are still
  required; removing all mounts indiscriminately is not a general solution.
- The inner build fell back to a root-mapped namespace without full UID mappings
  or a fakeroot helper. Package installation failed when APT tried to switch to
  its `_apt` user. Preparing dependencies outside this restricted namespace
  removed that step; the licensed installer then completed successfully.
- The prepared-base build completed in approximately six minutes, including
  official patch startup, SIF packaging, BLAST database/search validation and
  the Pathway Tools readiness marker. This duration excludes building the base.
- SRI returned HTTP 503 twice during the experiment. The successful test reused
  107 unmodified official patch files from an existing private SIF, checking
  their hashes against its retained manifest. A fresh vendor download was
  therefore not validated in that successful attempt.

### Build behavior and provenance

The inner build subprocess clears inherited bind variables, while leaving the
parent environment unchanged. The built image metadata records its installer,
recipe, image and official-patch hashes, plus a package-manifest hash for the
bundled dependency base. Validation remains blocking: a failed build is never
registered. Transient vendor HTTP errors receive at most three attempts; missing
patches still stop the build.

Use `--unsquash` on the outer Apptainer invocation if Java stalls on its FUSE
mount. This is a host launch setting; code inside an already-mounted container
cannot remount that outer image. Keep configuration inside the mounted project
when building from a container, and pass the container-visible SIF path with
`--image` when running it. A `/work` registration is meaningful inside that mount,
not from an unrelated host shell.

Maintainers can repeat the licensed checks on a target host from a source checkout:

```bash
python scripts/check_nested_ptools.py \
  --mp-image /absolute/path/metapathways.sif \
  --pt-image /absolute/path/pathway-tools.sif \
  --installer /absolute/path/pathway-tools-29.5-linux-64-tier1-install \
  --output /absolute/path/new-nested-test
```

Activate the host Apptainer environment first so its helper tools are on `PATH`.
Omit `--installer` to check only execution of an existing licensed image. The
output directory must be new; it contains private logs and a machine-readable
result. Each check has a configurable timeout (`--timeout`, default 900 seconds).
The test does not upload licensed files. Its configuration directory is isolated
under the requested output, so it does not replace the normal MP image registry.
Use `--unsquash` when testing the extraction route described above; the default
continues to test the mounted SIF.
A passed startup check verifies Pathway Tools readiness and BLAST integration;
it does not replace a complete PGDB inference test.
