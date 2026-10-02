# MetaPathways 3.5.2 release readiness

[Release instructions](releasing.md) · [Tester checklist](pr-testing.md)

Status updated 2026-10-02. This is a maintainer record of pre-release testing. Version remains 3.5.2. No release tag or registry publication is authorized by these checks.

## Registry and container status

Public metadata checked on 2026-10-02: [Anaconda](https://api.anaconda.org/package/hallamlab/metapathways) lists 3.5.1 as latest; the [latest GitHub release](https://api.github.com/repos/hallamlab/MetaPathways/releases/latest) is `v3.5.1-build1`; [Quay tags](https://quay.io/api/v1/repository/hallamlab/metapathways/tag/?limit=5) list 3.5.1 images. The 3.5.2 installation examples are release instructions, not evidence of publication. Until publication, testers install the reviewed checkout.

The current Conda recipe and core release image omit the MAGSplitter/Camelot VCS helpers. The user guides explicitly install pinned helpers for the complete three-sample test. Container instructions include a helper layer and writable test-reference mounts. These Docker/Apptainer walkthroughs still require end-to-end validation against the published 3.5.2 artifacts; documentation checks alone do not establish that validation. The source/Mamba three-sample test has been run successfully.

## Completed checks

- 125 runtime/regression tests and 17 release-control tests pass. Tests cover independent per-database taxonomy, SwissProt taxon parsing, database order, report joins, abundance exports, workflow scheduling, checkpoint behavior, sequence staging and release safeguards.
- The user completed the three-sample CAMI reviewer workflow locally with SwissProt and `--skip_ptools`: 51/51 successful tasks. [Detailed result audit](validation/reviewer-2026-10-02.json).
- All CDS records, mapped gene/RNA abundance rows, contig measurements and nine genome bins were reconciled. Report SQLite integrity and foreign keys pass. Expected RNA warnings and unavailable PGDBs are explained in the audit.
- CLI documentation and local documentation links pass validation.
- The earlier 3.5.2 candidate was built and installed successfully as a Conda package. Artifacts are specific to their source commit; rebuild from the final feature commit to include subsequent taxonomy and explorer changes.
- Source archives, installed native binaries and bundled three-sample fixture assets have been checked. Public packaging excludes licensed Pathway Tools installers, SIFs and MetaCyc databases.
- Feature pushes run smoke checks. Publication helpers require reviewed `dev`/`master`; no feature-branch publication is permitted.

## Before merge and public release

1. Build/install the final committed candidate and record artifact checksums alongside its source commit.
2. Complete the full local benchmark with the licensed Pathway Tools image and MPDB, retaining task traces and logs.
3. Complete the HPC Slurm test and benchmark with shared-storage paths, recording requested resources separately from measured consumption.
4. Have the tester review the feature branch/PR before merging to `dev` or `master`.
5. Run the release integration and optional container checks, verify registry credentials and metadata, and publish only after explicit release approval. Anaconda, Quay and Zenodo uploads are not part of the local test or feature push.
