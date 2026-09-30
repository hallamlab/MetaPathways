# Releasing MetaPathways

The source version in `metapathways/_version.py` is authoritative. An installed
`3.5.0.dev64` package is an older distribution; changing GitHub does not upgrade
that environment. The first stable release from this checkout is **3.5.0**.
Use **3.5.1**, **3.5.2**, etc. for subsequent patch releases, or **3.5.1rc1**
for a preview. Tags have a leading `v`; package versions do not.

The supported release package targets **Linux x86-64 / Python 3.10**.
It includes the core annotation pipeline and bundled FAST/metacount executables.
MAGSplitter, camelot-frs, and licensed Pathway Tools remain separately installed
optional dependencies. The release build does not download moving Git branches
and silently include them in the Conda package.

## Normal release using GitHub CI

From a normal checkout of `dev`, using your existing Git SSH/HTTPS credentials:

```bash
python scripts/release.py prepare 3.5.0
git diff
git add README.md metapathways/_version.py conda_recipe/meta_template.yaml
git commit -m "Prepare MetaPathways 3.5.0"
python scripts/release.py publish
```

Commit the release script, workflow, tests, and other implementation files as well
when installing this workflow for the first time. If `prepare` changes nothing,
there is no version-only commit to make.

`publish` requires a clean checkout on `dev`. It creates an annotated tag and
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
the `dev` branch. Assets are retained as an Actions artifact. Dispatching the
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
  python=3.10 conda-build conda-index anaconda-client pyyaml \
  pip wheel 'setuptools<81'
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
