"""Native accessible disclosures with linkable subsection labels."""
from html import escape
from docutils import nodes
from sphinx.util.docutils import SphinxDirective


class AdvancedOptions(nodes.General, nodes.Element):
    pass


class AdvancedDirective(SphinxDirective):
    has_content = True
    required_arguments = 1
    final_argument_whitespace = True

    def run(self):
        node = AdvancedOptions(summary='Advanced user options: ' + self.arguments[0])
        self.set_source_info(node)
        node.extend(self.parse_content_to_nodes(allow_section_headings=False))
        return [node]


def visit_html(translator, node):
    translator.body.append('<details class="mp-advanced"><summary>' +
                           escape(node['summary']) + '</summary>\n')


def depart_html(translator, node):
    translator.body.append('</details>\n')


def plain(translator, node):
    # Non-HTML builds display the full content without a disclosure control.
    pass


def register_panel_anchors(app, doctree):
    # MyST records heading slugs for page.md#section links. Panel titles use
    # rubrics to avoid introducing invalid section nesting inside disclosures.
    # Register their explicit IDs with the same metadata contract (MyST 4).
    slugs = app.env.metadata[app.env.docname].setdefault('myst_slugs', {})
    for panel in doctree.findall(AdvancedOptions):
        for rubric in panel.findall(nodes.rubric):
            for identifier in rubric['ids']:
                slugs[identifier] = (rubric.line, identifier, rubric.astext())


def setup(app):
    app.add_node(AdvancedOptions, html=(visit_html, depart_html),
                 text=(plain, plain), latex=(plain, plain))
    app.add_directive('advanced', AdvancedDirective)
    app.connect('doctree-read', register_panel_anchors)
    return {'version': '1', 'parallel_read_safe': True, 'parallel_write_safe': True}
