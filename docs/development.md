# Development and maintenance

For contributors changing MP, validating releases, maintaining references or
publishing packages and documentation. Normal analysis does not require these
steps. Start with the [user installation guide](installation.md) to run MP.

## Install from source

```bash
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
mamba env create -f docker/conda_base.yml
mamba run -n metapathways python -m pip install .
conda activate metapathways
```

Check out the revision you intend to test before installing. Reinstall after
changing the checkout; use `python -m pip install -e .` only when intentionally
working on an editable installation. Record `git rev-parse HEAD`. Do not update
code or its environment while an analysis is using it. Then follow the normal
[three-sample test](test.md).

## Contributor guides

```{toctree}
:maxdepth: 1

architecture
workflow-internals
ptools-maintenance
container-validation
validation-development
developer-cli-reference
pr-testing
releasing
release-readiness
documentation
```
