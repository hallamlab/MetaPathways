# Releasing MetaPathways

[User documentation](../README.md) · [Reproducibility](reproducibility.md)

This chapter is for maintainers publishing packages, not users installing or running MP. Its release examples describe the existing 3.5.x release process. Every release candidate must pass package and workflow validation before publication.

The source version in `metapathways/_version.py` is authoritative. An installed
`3.5.0.dev64` package is an older distribution; changing GitHub does not upgrade
that environment. The original stable release in this release series was **3.5.0**.
Use **3.5.1**, **3.5.2**, etc. for subsequent patch releases, or **3.5.1rc1**
for a preview. Tags have a leading `v`; package versions do not.

The supported release package targets **Linux x86-64 / Python 3.11**.
It includes the annotation pipeline, bundled FAST/metacount executables, reviewer
data, MAGSplitter and Camelot. Helper sources are pinned to exact commits with
SHA256 checksums. Licensed Pathway Tools is obtained separately by the user.

## Normal release using GitHub CI

After tester sign-off and merging the PR, from a clean checkout of the reviewed release branch (`dev` by default), using your existing Git SSH/HTTPS credentials:

```bash
python scripts/release.py prepare 3.5.2
git diff
git add README.md CITATION.cff metapathways/_version.py conda_recipe/meta_template.yaml
git commit -m "Prepare MetaPathways 3.5.2"
python scripts/release.py publish --branch dev
```

Commit the release script, workflow, tests, and other implementation files as well
when installing this workflow for the first time. If `prepare` changes nothing,
there is no version-only commit to make.

`publish` requires a clean checkout on the selected release branch (`dev` by default; `--branch master` selects `master`). It creates an annotated tag and
pushes the branch and tag atomically to `hallamlab/MetaPathways`. It never force
pushes, moves a published tag, or automatically commits your work. GitHub CLI
login is unnecessary on your machine.

The tag workflow:

1. Checks that the tag matches the source version and commit.
2. Creates a ZIP of the committed repository and a Python source tarball.
3. Renders a checksum-verified Linux Conda recipe and builds the package.
4. Installs that package from a temporary local channel into a new environment.
5. Runs database preparation and the included K12 integration test.
6. Requires all 16 named stages, nonempty key outputs, and no error markers.
7. Records logs, dependency exports, revision, and SHA-256 checksums.
8. Publishes a GitHub release with those downloadable assets.
9. Optionally uploads the same validated Conda package to Anaconda.org.

Only successful full builds are publishable. Release candidates are marked as
GitHub prereleases and use the `hallamlab/label/rc` Conda channel; stable versions
use `hallamlab` / label `main`. No wheel is advertised as platform-independent.

To test CI without publishing, select **Actions → Release → Run workflow** on
the feature branch (or the reviewed `dev`/`master` branch). Leave `release_tag` and `source_run_id` blank; optionally enable `test_containers` to build and scan Docker/SIF artifacts too. This does not publish packages, registry images, a GitHub release, or a DOI. Assets are retained as an Actions artifact. Dispatching the
workflow on a release tag also enables publishing.

## Enable Anaconda.org uploads

In GitHub repository **Settings → Secrets and variables → Actions**:

- Add the secret **ANACONDA_API_TOKEN**, with upload access to `hallamlab`.
- Set the repository variable **PUBLISH_CONDA** to `true`.

GitHub releases work without this configuration. GitHub authentication does not
authenticate Anaconda.org. The token is passed only to the upload step through
`BINSTAR_API_TOKEN`, never through a command argument or committed file.

## Local build and validation

Install packaging tools once; runtime dependencies are resolved in a separate
build/test environment:

```bash
mamba create -n metapathways-build -c conda-forge \
  python=3.11 conda-build conda-index anaconda-client pyyaml \
  pip wheel 'setuptools>=83,<85'
conda activate metapathways-build
python -m pip install build

python -m unittest discover -s tests/release -v
python scripts/release.py build
```

Run the build from a clean, committed checkout. Assets go into `dist/release/`.
The output directory must be empty; use `--output dist/release-attempt2` for a
retry. The build streams command names to the terminal and records tool output
in that directory, including failure diagnostics. Database preparation requires
internet access and downloads current taxonomy and enzyme data.

For a quick source-packaging check, use `build --source-only`; its artifacts
cannot pass the publication gate.

To upload the exact package after local validation:

```bash
anaconda login
python scripts/release.py upload-conda
```

Choose Anaconda.org if login asks for a destination. No force-overwrite option
is used. `--output` selects an alternative validated artifact directory.

## Rebuilds and failed releases

Use `prepare VERSION --build-number N` to change Conda's build number.
Without that option, preparing the same version preserves its build number;
preparing a new version resets the build number to zero. For a new
GitHub release, use a new version/tag; incrementing only the build number does
not permit replacing a published Git tag. The recipe template is authoritative;
generated `conda_recipe/meta.yaml` is not committed.

If a tag push fails, the local annotated tag remains. Fix the connectivity or
branch issue and retry `publish`; the script accepts that local tag only when
it still points to HEAD. If the tag reached GitHub, rerun the failed workflow
job instead of moving the tag. After a successful GitHub release but failed
Anaconda upload, rerun only the failed Anaconda job.

`make release-prepare VERSION=3.5.0`, `make release-build`, and
`make release-publish` wrap these commands. `make full-build` now builds,
validates, and uploads to Anaconda using this controller; it requires the
packaging environment above and does not publish a GitHub tag. Publishing
containers, PyPI wheels, and licensed Pathway Tools runs is outside this release
workflow. The K12 example validates installation and core operation; it does
not reproduce the manuscript's CAMI2 performance benchmark. Reference downloads
are moving resources, and dependency specifications are not a complete lock.
The recorded Conda export includes a temporary local URL for MetaPathways itself;
when recreating it, replace that URL with the downloaded release package.

## Recover the interrupted v3.5.0 GitHub upload

The failed upload attempted to attach an empty `global_errors_warnings.txt`.
Empty success logs now stay inside the complete release ZIP; only nonempty
files are standalone assets. Publishing resumes a partial draft and verifies
existing assets without overwriting them.

After committing and pushing the workflow fix to `dev`, open **Actions → Release
→ Run workflow**, select `dev`, and enter:

- `release_tag`: `v3.5.0`
- `source_run_id`: `36760751744`

This reuses the successful Conda build from that run and checks its source commit
against the original tag. Do not move the tag or run `publish` again for it.
An ordinary “Re-run jobs” on the old run still uses the old workflow.

## Docker, Apptainer, and Quay

Every tagged release (including manual recovery with `release_tag`) builds a
Linux amd64 Docker image from the exact validated Conda package and explicit
dependency export. It runs the core integration test in Docker, requires all 16
stages and six nonempty outputs, then converts that same image into an Apptainer
SIF. Apptainer is checked for the exact version and CLI startup; its read-only
filesystem is not used for the bundled database-writing integration test.
The SIF, container checksums, and validation logs are attached to GitHub.
They are also retained in the workflow's `container-assets` artifact.

In Quay, open the **hallamlab organization → Robot Accounts → Create Robot
Account**. Give the robot **Write** permission on `hallamlab/metapathways`.
In GitHub **Settings → Secrets and variables → Actions**, configure:

| Kind | Name | Value |
| --- | --- | --- |
| Secret | `QUAY_USERNAME` | Full robot username, e.g. `hallamlab+github_releases` |
| Secret | `QUAY_PASSWORD` | Robot token generated by Quay |
| Variable | `PUBLISH_QUAY` | `true` |
| Optional secret | `QUAY_API_TOKEN` | Quay OAuth token with `repo:write` access |

The robot credentials upload images. The OAuth token updates the repository
**description**, a separate Quay API permission. Create an OAuth application
under the Quay organization's **Applications**, then generate an access token
with permission to write repositories. Do not put credentials in source files,
Git commands, Docker build arguments, or Dockerfiles.

The overview is maintained in `docker/README.quay.md`. With `QUAY_API_TOKEN` set,
CI synchronizes it after an image push. Without that token, paste the file into
the repository description editor in Quay; image publishing still works.
You can update it locally with `python scripts/release.py quay-description` when
the token is available in your terminal environment.

Stable releases publish tags `3.5.0`, `v3.5.0`, and `latest` (substitute the actual
version). Release candidates publish their version tags without changing
`latest`. Version tags identify a release; record the Quay digest for an exact
image reference. Recovering an older stable release also updates `latest`, so
perform recovery before publishing a newer release.

For local builds after `python scripts/release.py build`:

```bash
python scripts/release.py container-build
# Local build verifies and tests both formats; these next commands publish:
docker login quay.io
python scripts/release.py container-push
python scripts/release.py container-attach  # requires GitHub CLI API authentication
```

Use `--ref v3.5.0` when using artifacts from an existing tag rather than HEAD.
Use `--output` for the package artifact directory and `--container-output` for
the container directory. Build requires an empty container output directory;
push/attach verify its checksums and provenance. Existing GitHub assets with
different checksums are never overwritten. Keep the validated container output
and local Docker image to retry a failed push/attachment without rebuilding.

Docker usage and Apptainer pull/exec examples are in `docker/README.quay.md`.
The older Makefile Docker targets and `docker/Dockerfile_ptools` are legacy
manual paths; release CI uses `scripts/containers.py` and `Dockerfile.release`.
Licensed Pathway Tools is not included in public release containers.

## Security rebuilds without changing the application version

The historical security rebuild of MetaPathways 3.5.1 used Python 3.11, Snakemake minimal
9.27.0, urllib3 >=2.8.0, and setuptools >=83. The minimal Snakemake distribution
provides the local CLI used for database preparation without the legacy stopit
runtime dependency. This records that release configuration; the current development workflow uses Nextflow, as described in the user guide. Container builds refresh the base image and apply Debian
package updates before installing the validated Conda environment.

A nonzero Conda build number produces a separate Git tag and release, for example
`v3.5.1-build1`; the application version remains `3.5.1`. This preserves the
original tag and published files. Prepare it with:

```bash
python scripts/release.py prepare 3.5.1 --build-number 1
# Commit the reviewed changes, then:
python scripts/release.py publish --branch dev
```

Quay receives `3.5.1-build1` and `v3.5.1-build1` tags as well as updated `3.5.1`,
`v3.5.1`, and stable `latest` aliases. The new GitHub release carries its own SIF,
archives, checksums, and dependency exports. Conda gets a new build rather than
an overwritten package. For reproducible use, choose a build tag or image digest.

Pip is a build-time tool, not a dependency of the packaged runtime. The validated
runtime is created with Conda's automatic pip insertion disabled; dependency
exports use Python package metadata. This avoids shipping vulnerable libraries
vendored inside pip while leaving normal package-building tools available to CI.

Before publishing containers, CI runs checksum-pinned Trivy 0.74.0 directly on
the runner and blocks publication when it finds a high/critical vulnerability
for which a fix is available. Both the full scan and filtered gate reports are retained as
`container-security-report`. Unfixed distribution advisories need separate
review; passing this gate does not mean the image has no vulnerabilities.

## Feature-branch review before publication

The intended sequence is **feature branch → single-server and HPC benchmarks plus tester review → PR merge to dev/master → release tag → publication**. Pushing a branch or opening a PR does not publish a release. Smoke checks run for PRs targeting `master`, `dev`, or a future `main`, and pushes to those branches or `feat/**`. They validate the citation file, unit tests, docs, source/wheel builds, installed CLI entry points, and packaged reviewer/report assets.

Use the [PR tester checklist](pr-testing.md) for independent installation, the single- and two-sample CAMI workflows, restart behavior, optional licensed Pathway Tools, and report exports. Record the tested Git commit. Record the source commit as well as the package version when testing a candidate.

The selected next release is **3.5.2**. Preparing this version does not publish it; branch and PR testing still happen before a release tag is pushed. Do not reuse 3.5.1 or its existing tags. Change the version with `prepare` only at the release step; it also updates the citation version and resets the Conda build number for a new version. Avoid editing runtime code/version files while a benchmark is executing from an editable checkout.

Before merging, manually dispatch the Release workflow on the feature branch with blank publication inputs. Select `test_containers` to test the exact Conda artifact inside Docker and its SIF conversion. Download `release-assets`, `container-assets`, and the security reports from Actions for review. Core integration uses explicit 2 GB task reservations, a 4 GB memory budget, and two CPUs so the small fixture fits CI runners. These limits are not production metagenome sizing recommendations.

The Conda package includes the Nextflow controller, reporting assets, reviewer data, MAGSplitter and Camelot. Helper sources and SHA256 checksums are pinned in `requirements-workflow.txt` and fetched during the package build. The Quay image installs that same artifact; users do not install helpers separately. Source installations resolve the same pinned helpers through MP’s Python package metadata. The public build gate exercises the core K12 workflow; **passing that gate does not certify `analysis_wf` with MAGs, Slurm, nested Apptainer, or licensed Pathway Tools**. Those require the tester checks. A container containing Apptainer does not guarantee the host permits nested container execution. Test licensed inference with the recommended host Conda installation and your own image first.

## Zenodo and software citation

`CITATION.cff` supplies software authors (from MP's existing author metadata), title, repository, and license. Confirm the author list before publication. No DOI or release date has been invented. `prepare` adds/updates the selected software version. We maintain one citation metadata file; Zenodo prioritizes `.zenodo.json` over CFF if both exist, so do not add a conflicting second file. See [Zenodo's supported metadata](https://help.zenodo.org/docs/github/describe-software/).

A repository administrator must link their GitHub account to Zenodo, grant the required organization access, and enable `hallamlab/MetaPathways` in Zenodo's GitHub settings. See [enable a repository](https://help.zenodo.org/docs/github/enable-repository/). This connection is account-side configuration and cannot be inferred from committed files or the existence of an earlier release. Our GitHub Actions workflow does not call the Zenodo API or store a Zenodo token.

Keep pre-merge tests as branches, PRs, and Actions artifacts. Once enabled, the integration ingests new releases; do not publish a release candidate casually if you are not ready for it to become an archive record. After the intended release, verify its Zenodo record, version, source commit/tag, author metadata, license, and DOI before announcing it. Record both the version-specific DOI and concept DOI where applicable, and use the version-specific DOI for an exact release citation.

Treat the software snapshot, built packages/container assets, and manuscript benchmark data as separate deliverables. Do not assume the GitHub integration copies every release attachment or later-added SIF, nor that it deposits your benchmark outputs. Inspect the archive contents and arrange a separate data deposit with its own provenance if required. Never deposit the licensed Pathway Tools installer, SIF, or MetaCyc databases as public MP assets.

## Account-side release checklist

| Service | Repository preparation | Maintainer verification before publishing |
| --- | --- | --- |
| GitHub | PR smoke checks, manual package/container validation, tag workflow | Protect the target branch, require smoke checks and tester approval; choose a new version |
| Anaconda.org | Recipe, installed-package integration, checksums, RC/main labels | `ANACONDA_API_TOKEN` has upload access to `hallamlab`; `PUBLISH_CONDA=true` |
| Quay | Build from validated package, Docker test, SIF conversion, vulnerability gate | Repository exists; robot has Write permission; `QUAY_USERNAME`, `QUAY_PASSWORD`, `PUBLISH_QUAY=true` |
| Zenodo | CFF metadata and documented release process | GitHub integration enabled, correct organization authorization, archived release/DOI verified |

Secret values must never be committed or printed. Repository files alone cannot verify the current account permissions or secret configuration. A local source build is not proof of a successful Conda solve, registry upload, or Zenodo deposit.
