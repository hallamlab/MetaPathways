# Reports and the EDA portal

Successful `run`, `mag_split` and `ptools` commands refresh the report automatically. When a sample belongs to an existing parent report, the parent report is refreshed so all samples remain accessible. Report indexing runs on the controller after Nextflow completes, uses SQLite on disk, and does not repeat annotation, mapping or pathway inference.

For existing or partial results, build a report independently:

```bash
metapathways report -o results
```

`-o` can be one sample output directory or a parent containing sample directories as immediate children. Reports are saved under:

```text
results/reports/
├── MP_run_report.html       # Accounting, import notes, files and execution links
├── EDA_portal.html          # Search and export interface
├── results.sqlite          # Relational data for all included samples
├── schema.json             # Columns, row counts, keys and relationships
├── output_inventory.tsv    # Paths, sizes and modification times
├── portal.js
└── portal.css
```

Open `MP_run_report.html` directly for navigation. To search, subset and export:

```bash
metapathways report -o results --serve --no-rebuild
```

This opens a browser and serves the existing snapshot on loopback only. Omit `--no-rebuild` to refresh changed results first. `--no-browser` prints the URL; `--port 8765` chooses a fixed local port. Ctrl-C stops the server. All assets are local, and the explorer makes no third-party requests. Copy the complete output directory to a workstation to retain file links; no cloud upload is involved.

Protein taxonomy supports SwissProt (including the small test database), UniRef and eggNOG. Each annotation row reports the taxon of its own reference hit in `taxonomy`; `lca_taxonomy` is computed from score-qualified hits in that same database, with support counted independently for each database. The primary ORF row follows its selected annotation's `reference_db` and target. No database priority or cross-database fallback is used. These describe reference-hit evidence, not a definitive organism assignment for the query. Unsupported databases report `Not computed`; missing or unknown taxon IDs in supported databases report `Unclassified`. An actual LCA of `root` remains `root`. SILVA rRNA taxonomy is separate. Existing outputs need annotation-table regeneration to gain these fields.

The workflow in the portal is:

1. Choose a table: ORFs/taxonomy, functional annotations, pathways, pathway genes, MAG membership, abundance or file inventory.
2. Search text, add column filters, or filter ORF-based tables by linked EC/reaction, reference database and pathway/entity identifiers.
3. Follow a row's related-result buttons. Sample, pathway and entity keys remain attached, avoiding cross-sample collisions.
4. Select export columns, apply filters, and export **all matching rows** as CSV. Pagination does not truncate an export.
5. Save the query definition or bookmark its URL to record your selection.

The portal is a table explorer: it adds no biological plots or interpretation. The [schema guide](results-schema.md) explains row units, relation keys, safe joins and the distinction between missing output and biological absence. Large joins remain in SQLite; CSV is streamed rather than assembled entirely in browser memory. Broad scans can still take time. Each query has a two-minute execution budget, and at most four queries execute concurrently.

The full inventory links remaining RNA, sequence, alignment, statistics and raw PGDB products. The relational views cover the recognized formats documented in the schema guide; an arbitrary legacy table is not silently assumed to fit that schema. Missing primary annotations leave explicit placeholder ORFs and import notes; malformed recognized tables stop the rebuild and preserve the previous database.
