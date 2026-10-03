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

Mermaid sources live in `docs/diagrams/`; GitHub and Read the Docs display the committed SVG previews at a shared scale. Keep the conceptual diagram in the README consistent with `overview.md`. The GitHub documentation job builds with warnings treated as errors.

## Connect Read the Docs

In the existing **metapathways** project, connect `hallamlab/MetaPathways`. The configuration file is `.readthedocs.yaml`, the Sphinx configuration is `docs/conf.py`, and dependencies are in `docs/requirements.txt`.

Activate the branch you want to preview under **Versions** and trigger a build. Keep the public default on the intended main/development branch; switch it to the updated branch after the documentation is merged. Enable pull-request previews if desired. Branch previews allow review before changing the site's default version.

The project address stays **https://hallamlab-metapathways.readthedocs.io/**. Repository changes supply the build configuration; linking the project, version activation and default-version settings are managed in Read the Docs. Confirm a successful build there before announcing updated hosted pages.

The previous `quick_start.html`, `install.html` and `usage.html` addresses redirect to the corresponding new pages.

## Shared diagram scale

ASPIRE and MetaPathways use a shared 2,240-unit-wide white canvas for documentation previews. MP’s main workflow in the README and documentation home page uses its cropped SVG at full available width for readability. Smaller diagrams are centered without stretching; Mermaid diagrams use a common 1.5× scale to bring their 16-pixel labels close to the publication figures’ typography. Preview width is responsive, but relative scale stays consistent across pages. Click a diagram to open its SVG for closer inspection. Original publication SVG/PDF downloads stay tightly cropped.

Edit Mermaid sources under `docs/diagrams/`, or the original publication SVGs under `docs/assets/`. To rebuild the committed previews from the repository root, in the documentation environment:

```bash
python -m pip install -r docs/diagram-requirements.txt
python -m playwright install chromium
python scripts/render_workflow_diagrams.py
```

Chromium requires its usual Linux system libraries. Rendering downloads the pinned Mermaid bundle; ordinary Sphinx builds need neither Chromium nor network access for diagrams. `docs/diagrams/figures.json` records the source mapping and shared canvas size. If a future diagram needs a wider canvas, update both projects together. Do not hand-edit generated files in `docs/assets/diagrams/`.
