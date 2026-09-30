# MetaPathways

MetaPathways is a pipeline for processing and annotating assembled metagenomic sequences.

- **Source and documentation:** https://github.com/hallamlab/MetaPathways
- **Releases and Apptainer downloads:** https://github.com/hallamlab/MetaPathways/releases
- **Issues:** https://github.com/hallamlab/MetaPathways/issues
- **License:** MIT
- **Platform:** Linux x86-64 (amd64), Python 3.10

Release containers contain the validated MetaPathways Conda package and its pinned dependency versions. They include the core annotation pipeline; optional MAGSplitter, camelot-frs, and licensed Pathway Tools require separate installation.

## Docker

Use a version tag for reproducible work (replace VERSION with a published release, such as 3.5.0):

```bash
docker pull quay.io/hallamlab/metapathways:VERSION
docker run --rm quay.io/hallamlab/metapathways:VERSION metapathways version
docker run --rm -v "$PWD:/work" -w /work --user "$(id -u):$(id -g)" quay.io/hallamlab/metapathways:VERSION metapathways run --help
```

Mount your input files, configuration, reference databases, and output directory when running a pipeline. The `latest` tag tracks stable releases; release candidates do not update it.

## Apptainer / Singularity

Download the `.sif` file and container checksums from the corresponding GitHub release, or convert the Docker image:

```bash
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:VERSION
apptainer exec metapathways.sif metapathways version
apptainer exec --bind "$PWD:/work" --pwd /work metapathways.sif metapathways run --help
```

The image contains software, not production reference databases. Follow the source repository's documentation to configure your data and databases.
