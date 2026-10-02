# Nextflow feature release readiness

[Release instructions](releasing.md) · [Tester checklist](pr-testing.md)

This records the local preparation audit on 2026-10-02 for `feat/nextflow-controller-db-build`. It is not a published-release certificate. The selected next version is 3.5.2, prepared locally at the user’s request. No release tag has been published.

## Prepared and locally checked

- PR/push smoke checks cover `master`, `main`, `dev`, and feature-branch pushes.
- Release publishing requires a clean `dev`/`master` checkout; a feature branch cannot use the publish helper.
- Manual Release dispatch with blank publication inputs builds/test artifacts without public uploads. `test_containers` additionally validates Docker/SIF and scans the image.
- Conda recipe includes Nextflow and Apptainer, expanded CLI checks, and an exact source SHA256. Rendering passed; full package solve/build remains a CI check.
- Core CI integration has explicit small-fixture resource budgets.
- Source archive and Linux x86-64 CPython 3.11 wheel built successfully from an isolated copy of the working tree. The native-binary wheel is no longer marked `py3-none-any`.
- Wheel installed in a temporary environment; all eight CLI entry/help checks passed outside the source checkout. This environment shared runtime dependencies with the existing environment; it was not a fresh Conda solve.
- Installed workflow modules, native binaries, report assets, all three CAMI manifests, and fixture SHA256 hashes passed inspection.
- 102 runtime/unit tests passed across the initial run and the report-test rerun with localhost permission; 17 release-control tests passed. These include fake local Git remotes; no actual publication occurred.
- CFF schema validation passed; release preparation updates its version. Documentation links, CLI examples and shell syntax passed.

## Required before merge/release

1. Commit and push the reviewed feature work, open a PR targeting `dev` or `master`, and record its exact commit. No branch push, PR creation or publication was performed by this local audit.
2. Run PR smoke checks and the nonpublishing Release workflow on that commit. Confirm a fresh Conda solve/build/install, K12 integration, Docker integration, SIF startup and security scan. These full builds were not run locally during the benchmark.
3. Obtain independent tester sign-off using the CAMI single/pair commands, checkpoint reuse, explorer exports, and any applicable Pathway Tools/Slurm checks. MAGSplitter and Camelot remain separately installed dependencies; the public container does not certify the complete licensed workflow.
4. Verify account-side Anaconda/Quay credentials and permissions, GitHub branch protection, and the Zenodo repository connection. None of these settings were verified locally. Review software authors and third-party fixture redistribution terms before public release.
5. Use the selected version 3.5.2, merge after testing, and publish only through the documented release process. Confirm actual uploaded artifacts, checksums/digests and Zenodo record contents/DOI afterward.

The source/archive audit artifacts are local temporary files under `/data/tmp/mp-package-audit-xpb6ytrx`; retain CI artifacts for the actual reviewed commit instead of treating this uncommitted snapshot as the release.
