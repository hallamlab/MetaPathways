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

## Deferred dependency security remediation (2026-10-09)

The [4.0.0 validation run](https://github.com/hallamlab/MetaPathways/actions/runs/38000115287)
passed Conda installation and workflow integration, Docker workflow integration,
and MP Apptainer artifact checks. Publication remains blocked by the repository's
container security gate. This is a HallamLab release policy, not an Anaconda
upload requirement.

Trivy 0.74.0 reported 649 fixable high/critical occurrences across 89 distinct
vulnerability IDs, primarily in Apptainer and its CNI networking dependencies,
plus Java libraries bundled with Nextflow. Repeated occurrences across binaries
are not independent vulnerabilities. Individual applicability to MP has not
been established. The retained `container-security-report` artifact contains the
full and filtered JSON reports. The licensed Pathway Tools SIF was not scanned.

The tested packages were Apptainer 1.5.4, CNI 1.0.1, CNI plugins 1.3.0 and
Nextflow 26.04.7. An upstream Nextflow 26.09.2-edge candidate had no fixable
high/critical Java findings in a separate scan, but remains an untested prerelease
for MP. Official Apptainer 1.5.4 Debian binaries still had findings. No candidate
replaced the validated dependencies.

Custom dependency rebuilds are deferred. Track patched upstream packages,
remaining findings, and required compatibility/security validation in
[issue #15](https://github.com/hallamlab/MetaPathways/issues/15). Deferral does not
disable the security gate, suppress findings, or authorize publication.

## Release checks

1. Pass unit, release-control, documentation, source and wheel checks.
2. Build the final committed Conda package and run its installed workflow test.
3. Build/test Docker and Apptainer artifacts and review the vulnerability scan.
4. Retain checksums, dependency exports and the tested commit for each artifact.
5. Create the matching version tag after production validation. Manually select
   Anaconda and Quay publication only when ready; both default to off.
6. Publish a GitHub release separately when DOI archiving is intended.

See the [release process](releasing.md) for exact commands and recovery steps.
