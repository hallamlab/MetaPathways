# Explore results, follow links and export tables

[Home](index.md) · [Exact schema and SQL](results-schema.md) · [Benchmark statistics](benchmarking.md)

## Open the report or start the explorer

The report and explorer serve different purposes. `MP_run_report.html` provides accounting and navigation to results/logs. `EDA_portal.html` queries the SQLite index so you can search, filter, follow related records and export CSVs. The explorer requires its local Python server; double-clicking the explorer HTML alone does not start that server.

After a completed workflow:

```bash
metapathways report -o /path/to/analysis --serve --no-rebuild
```

Keep the terminal open. MP prints a local address and normally opens your browser. Ctrl-C stops the server. To refresh the index after results change, omit `--no-rebuild`. Report generation reads existing output; it does not rerun annotation or infer missing pathways.

Use the parent directory containing sample folders to combine samples in one report. Using a single sample directory limits that report to the sample. Report files live under the selected root's `reports/` directory.

## View a remote report through SSH

On the **remote server**, activate the MP environment and run:

```bash
metapathways report -o /remote/path/to/analysis \
  --serve --no-rebuild --no-browser --port 8765
```

On **your local computer**, in a second terminal, replacing `user@server` with the same SSH destination you normally use:

```bash
ssh -N -L 8765:127.0.0.1:8765 user@server
```

The tunnel command may appear idle; that is normal. Leave both terminals running. Open the token-bearing URL printed by MP in your local browser, preserving its query string. The report paths are:

- Explorer: <http://localhost:8765/reports/EDA_portal.html>
- Run report: <http://localhost:8765/reports/MP_run_report.html>

`-L` forwards the local port to the server's loopback port, and `-N` asks SSH to provide forwarding without a remote shell. Use your site's SSH host alias or jump-host configuration when required. MP listens only on loopback; you do not need to expose the server publicly or open a firewall port.

If local port 8765 is busy, forward another local port with `-L 8766:127.0.0.1:8765`, then browse to `localhost:8766`. If the remote port is busy, choose a different `--port` and update the tunnel's destination port too. A page that loads but cannot query may be an old static copy or an HTML file opened without the server.

## Start with accounting

Open the run report before interpreting tables. Review which samples and products were found, source/import notes, and task execution records. Optional MAG failures can coexist with overall workflow success. A file from an earlier successful attempt can coexist with a failed newer attempt; output availability and latest task status answer different questions.

Check a few known sample IDs and original contig IDs. If expected tables are empty, inspect the inventory and task logs before interpreting zeros. Missing output does not mean the organism lacks a function.

## Follow one annotation through the results

1. Choose **ORFs and taxonomy** and filter to a sample.
2. Search for an ORF ID or product term you recognize.
3. Follow its related annotation records to inspect database accessions and EC/reaction assignments.
4. Follow pathway associations, keeping both sample and entity IDs selected.
5. Export the filtered table with identifying columns included.

A primary ORF annotation is not identical to all database hits. **Functional annotations** retains reference records; a gene can have more than one. **Pathway genes** contains pathway/ORF associations, so a gene in several pathways appears several times. The [schema](results-schema.md) describes each row unit and join key.

## Select useful subsets

| Question | Starting table and filters |
| --- | --- |
| Which SwissProt annotations mention kinase in this sample? | Functional annotations: sample, reference database, product text |
| Which genes support a particular MAG pathway? | Pathways: sample, entity, pathway ID; follow Pathway genes |
| What are all reported genes on contigs assigned to a MAG? | MAG ORFs: sample and entity |
| Which genes were supplied to MAG pathway inference? | MAG input genes: sample and entity |
| What abundance measurements exist for an ORF? | Abundance: sample, feature type, feature ID, measurement |
| Where is the raw GFF, RNA output, BAM or PGDB archive? | File inventory / run-report output links |

**MAG ORFs** uses the full contig map. **MAG input genes** is the smaller set selected for Pathway Tools; it is not the total gene content of the bin. Similarly, input genome-bin counts need not equal successful PGDB counts.

The portal provides text searches and column/related-result filters, not biological plotting. Choose the columns needed for downstream work and export all matching rows as CSV. Pagination limits what is displayed, not the full filtered export. Retain the query definition or bookmarked URL and record the report snapshot used for an analysis.

## Review abundance measurements side by side

Select **Read abundance** for one row per sample and feature, with separate
columns for count, length, coverage, RPKM and TPM. Use **Choose table and export
columns** to turn measurements on or off; the table updates immediately and CSV
exports use the same selection. Missing measurements remain blank. Filter
`feature_type` to compare contigs or ORFs separately.

The raw measurement table is retained in SQLite for provenance and older saved
queries, but is omitted from the table menu. Source files remain available through
the file inventory. To update an existing report after installing the new explorer,
run `metapathways report -o /path/to/analysis --serve` without `--no-rebuild`, then
reload the browser page. This rebuilds the report without rerunning the analysis.

## Avoid accidental double counting

Always keep `sample_id`. Keep `entity_id` for pathways and genome results. The same ORF-looking identifier can occur in different samples, and the same pathway can occur in the community and many bins.

Do not sum abundance after joining each gene to all its annotations and all its pathways: many-to-many joins multiply rows. Decide whether you need unique genes, annotation records or pathway associations before aggregating. Pathway scores are reported inference values, not universal probabilities of biological truth. See [safe SQL examples](results-schema.md#sql-with-explicit-joins).

## Share or archive a result

Filtered CSVs can be shared independently. To retain the explorer and original-file navigation, preserve the complete output root, including its `reports/`, sample results and logs, and serve it with a compatible MP installation. `results.sqlite` and `schema.json` also support direct read-only analysis in Python/R/SQLite.

Not every raw product becomes a relational table: the inventory links additional sequences, alignments, RNA reports and PGDB products. Import notes identify missing or unsupported sources. The portal does not silently infer missing data, recalculate abundance or validate historical mapping provenance.

### Starting at the sample level

A fresh explorer URL opens **Samples**. Use a sample row's related-results buttons to inspect its ORFs, then follow annotations, pathways and abundance. Saved query URLs retain their selected table and filters. The footer shows the recorded run version, sample count, command, executor and status when available, plus the report update time and a GitHub link for bug reports and feature requests. For older outputs without a recorded run version, the footer labels the report-builder version instead. Diagnostic details remain available in the Import notes and Execution history tables.

## Follow RNA features to contigs and MAGs

Choose **RNA features** to browse 16S/23S rRNA and tRNA loci. Filter `subunit`
for a particular rRNA, or `feature_type` for tRNAs. tRNA rows retain their
reported isotype and anticodon; tRNAs do not receive a SILVA taxonomy assignment.

Use **RNA taxonomy** for an rRNA locus's reference hits. **Contig ORFs** opens
other genes on that contig; **Contig MAGs** shows its recorded bin membership.
From a gene, annotation or pathway-gene row, **Contig RNA** and **Contig rRNA
taxonomy** show the same context. From a MAG, use **MAG RNA**. **RNA abundance**
shows measurements when matching feature IDs are present in the read-count table.

RNA taxonomy is evidence on a locus and its contig, not an automatic taxonomic
assignment to every gene or the entire MAG. Different loci or references can
report different taxonomy; MP preserves these separately. Missing evidence is
not a negative result. Unmatched taxonomy records remain visible with an
explanation in **Import notes**.

After updating MP, rebuild an older explorer with `metapathways report -o OUTPUT`
(or add `--serve` to browse immediately). This reads retained output files without
rerunning analysis. Older compact outputs may lack the RNA source files; rebuilding
cannot recover files that were already removed.

## Download results for one MAG

Choose **Sample**, then **MAG** at the top of the explorer. The selection stays
active as you switch biological tables, including contigs, genes, all features,
RNA evidence, annotations, abundance and pathways. **All features** combines
annotated protein-coding genes with RNA loci. Run-wide metadata is clearly marked
and is not filtered by this selection.

Use **Export all matching rows as CSV** for the current table and filters.
Use **Download MAG results (ZIP)** for the complete indexed results for that MAG,
independent of the current table filters. The ZIP contains separate CSVs, source
checksums, a row-count manifest and a README. It contains tables, not FASTA
sequences or complete PGDB archives. Empty tables can indicate unavailable outputs.

MAG membership comes from the preserved contig-to-MAG map. A MAG without that map
cannot be exported as a complete subset. Pathways are the results inferred for
that MAG; community pathways are not copied into the MAG export. The optional
**Pathway Tools entity ID** linked filter still refers specifically to PGDB input
genes; use the top-level MAG selector for all mapped features.
