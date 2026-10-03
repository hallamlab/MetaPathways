"""Build the guides without importing MP or installing analysis tools."""
import ast
import os
from pathlib import Path

project = 'MetaPathways'
author = 'Hallam Lab and MetaPathways contributors'
copyright = '2026, MetaPathways contributors'
metadata = ast.parse((Path(__file__).parents[1] / 'metapathways/_version.py').read_text())
release = next(ast.literal_eval(n.value) for n in metadata.body
               if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == '__version__' for t in n.targets))
version = release
extensions = ['myst_parser', 'sphinxcontrib.mermaid']
source_suffix = {'.md': 'markdown', '.rst': 'restructuredtext'}
root_doc = 'index'
exclude_patterns = ['_build', 'build', 'src', 'validation', 'requirements.txt', 'includes']
myst_heading_anchors = 4
myst_fence_as_directive = ['mermaid']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {'collapse_navigation': False, 'navigation_depth': 2}
html_title = f'MetaPathways {release}'
html_baseurl = os.environ.get('READTHEDOCS_CANONICAL_URL', 'https://metapathways.readthedocs.io/en/latest/')
html_static_path = ['assets']
html_css_files = ['docs.css']
mermaid_version = '11.12.1'
mermaid_init_config = {'startOnLoad': False, 'theme': 'neutral', 'flowchart': {'htmlLabels': False}}
mermaid_fullscreen = True

# Preserve the entry points from the previous Sphinx site.
templates_path = ['templates']


def legacy_pages(app):
    for old, new in {'quick_start': 'installation', 'install': 'installation', 'usage': 'commands'}.items():
        yield old, {'redirect_target': new + '.html'}, 'redirect.html'


def setup(app):
    app.connect('html-collect-pages', legacy_pages)

# The extension runs Mermaid on window load; avoid a second automatic render.
mermaid_light_theme = 'neutral'
mermaid_dark_theme = 'neutral'
mermaid_height = 'auto'
