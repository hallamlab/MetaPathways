# Advanced user guide

Use this section when you need custom inputs, cluster tuning, SQL access or
controlled comparisons. **Advanced user does not mean administrator:** these are
optional analysis settings, not sudo privileges. The default workflow is in the
[installation](installation.md), [database](databases.md), [analysis](analysis.md)
and [Pathway Tools](pathway-tools.md) guides. Expand the matching **Advanced user
options** panels on those pages only when needed.

- [Custom sample manifests](analysis.md#custom-analysis-manifest), including long reads and sample subsets.
- [Alternative licensed database sources](databases.md#already-have-a-licensed-metacyc-source).
- [Taxonomic scope](pathway-tools.md#choose-a-taxonomic-scope) and [transport inference](pathway-tools.md#transport-inference-and-sequence-backed-inputs).
- [Stage reruns](annotation.md#stage-controls) and [compact outputs](benchmarking.md#compact-results-on-limited-storage).

```{toctree}
:maxdepth: 1

resources
results-schema
benchmarking
cli-reference
```

Code changes, fixture regeneration, compatibility-list maintenance and release
publication belong in [Development](development.md).
