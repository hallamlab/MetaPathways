# Reproducibility and validation

[User guide](index.md) · [Advanced user guide](advanced.md)

## Record the actual revision and inputs

MetaPathways is distributed under the {download}`MIT license <../LICENSE>`; bundled third-party sources retain their notices. GitHub tags, Conda packages and container images are separate publications. A tag alone does not update the other distributions. Record the source Git commit, package build or immutable container digest used for each analysis.

Preserve:

- Assembly/read identifiers, layout (single, paired or interleaved), checksums and sample naming.
- Original contig-to-MAG maps, including the binning version/method and QC decisions.
- MP commands, workflow parameters, tool versions and environment exports.
- Reference releases, acquisition dates and file checksums; public downloads can change.
- The Pathway Tools installer/image checksum, version and taxonomic-pruning setting.
- Final outputs, run logs, Nextflow trace/report/timeline and report `sources`/`schema.json`.

Example environment records:

```bash
mamba list --explicit > conda-explicit.txt
metapathways version > metapathways-version.txt
```

For a source installation, also record `git rev-parse HEAD` from the checkout and `python -m pip freeze` from its activated environment.

## What the installation test establishes

The bundled three-sample walkthrough checks that your installed workflow runs on
small inputs. It skips licensed Pathway Tools by default and downloads current
support records. Record the download dates; a successful example run is not an
accuracy benchmark for your own biological dataset.

For measured comparisons, use the [benchmarking guide](benchmarking.md), retain
the original traces, and distinguish fresh runs from resumed work. Reports built
from existing outputs cannot reconstruct measurements that were not recorded.
