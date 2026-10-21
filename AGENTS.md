# Repository instructions

## Workflow figure aesthetic — standing user preference

Read `docs/WORKFLOW_STYLE.md` before changing workflow diagrams, their generators,
exports or documentation embeds. Preserve the MetaPathways nodal publication
style: numbered square module spine, diamond compute nodes, circular input /
data / output nodes, Times serif labels, thin black connectors and the shared
blue / gray / green palette. Update scientific content within that layout.
Do not replace these figures with stacked text cards or generic box flowcharts.
A style change requires an explicit user request.

Rebuild source-linked SVG / PDF / preview artifacts together, run
`python3 scripts/check_workflow_style.py`, and visually inspect the result.
This governs workflow figures; preserve the separately configured study palettes
for scientific plots. Apply this convention to new workflow figures too.

MetaPathways-specific user override (2026-10-21): the primary brief and detailed
workflows use standard computational flowchart shapes with the existing palette
and Trebuchet MS typography. See docs/WORKFLOW_STYLE.md for the scope.

The approved nodal brief is now the README and RTD landing-page figure; the
standard flowcharts remain in the workflow section. Follow the nodal-specific
layout and styling rules in docs/WORKFLOW_STYLE.md when editing that figure.
