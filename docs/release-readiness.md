# MetaPathways 4.0.0 release readiness

[Release instructions](releasing.md) · [Tester checklist](pr-testing.md)

The production branch is `main`; `dev` is the integration branch. Installation
examples use version 4.0.0. A version update or branch push does not publish a
package, image, GitHub release or DOI.

## Validation evidence

The corrected 49-sample HPC benchmark completed under version 3.5.2 with 6,519
successful uncached tasks, 49 community PGDBs and 5,392 population PGDBs. The
4.0.0 release must retain its own tested source commit and package/container
validation; changing the version does not relabel that historical benchmark.
The earlier three-sample installation audit is retained as
[historical validation evidence](validation/reviewer-2026-10-02.json).

The Conda recipe builds pinned MAGSplitter and Camelot sources into the package.
The Quay image uses that same validated Conda artifact. The bundled test data
are included; licensed Pathway Tools, its images and production reference data
are excluded from public packages.

## Source validation (2026-10-09)

- 170 runtime tests and 20 release-control tests passed.
- Strict documentation build, links/anchors in 40 rendered pages, and citation metadata passed.
- Version 4.0.0 source distribution and Linux/Python 3.11 wheel built successfully.
- Installed-wheel assets, version, bundled test preparation, MAGSplitter CLI and Camelot import passed outside the source checkout.
- Final Conda integration and Docker/Apptainer validation remain release-build checks.

## Release checks

1. Pass unit, release-control, documentation, source and wheel checks.
2. Build the final committed Conda package and run its installed workflow test.
3. Build/test Docker and Apptainer artifacts and review the vulnerability scan.
4. Retain checksums, dependency exports and the tested commit for each artifact.
5. Create the matching version tag after production validation. Manually select
   Anaconda and Quay publication only when ready; both default to off.
6. Publish a GitHub release separately when DOI archiving is intended.

See the [release process](releasing.md) for exact commands and recovery steps.
