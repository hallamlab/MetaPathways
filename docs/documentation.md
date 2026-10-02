# Maintaining the documentation

The Markdown files in `docs/` are the source for Read the Docs. Keep the README limited to the original abstract, installation/reviewer essentials, the conceptual diagram and links to the complete guide. Edit detailed instructions here; do not maintain separate copies on the website.

## Build and review locally

Use a separate documentation environment:

```bash
mamba create -n metapathways-docs -c conda-forge python=3.11 pip
conda activate metapathways-docs
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m http.server 8000 --directory docs/_build/html
```

Open `http://localhost:8000`. The documentation build does not install MP, reference databases, or Pathway Tools. It uses the committed CLI reference; after changing CLI arguments, regenerate that page from an MP development environment with `python scripts/generate_cli_docs.py`.

Mermaid diagrams use ordinary fenced `mermaid` blocks, which render on GitHub and through the Sphinx extension. Keep the conceptual diagram in the README consistent with `overview.md`. HTML diagrams load the pinned Mermaid JavaScript version from jsDelivr. The GitHub documentation job builds with warnings treated as errors.

## Connect Read the Docs

In the existing **metapathways** project, connect `hallamlab/MetaPathways`. The configuration file is `.readthedocs.yaml`, the Sphinx configuration is `docs/conf.py`, and dependencies are in `docs/requirements.txt`.

Activate the branch you want to preview under **Versions** and trigger a build. Keep the public default on the intended main/development branch; switch it to the updated branch after the documentation is merged. Enable pull-request previews if desired. Branch previews allow review before changing the site's default version.

The project address stays **https://metapathways.readthedocs.io/**. Repository changes supply the build configuration; linking the project, version activation and default-version settings are managed in Read the Docs. Confirm a successful build there before announcing updated hosted pages.

The previous `quick_start.html`, `install.html` and `usage.html` addresses redirect to the corresponding new pages.
