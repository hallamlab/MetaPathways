# Workflow and output guide

[Start with the README](index.md) · [Complete CLI flags](cli-reference.md)

For guided usage, start with the [command cookbook](commands.md) or [Pathway Tools chapter](pathway-tools.md). This page describes the stage and output contracts.

**Before running the complete workflow with PGDBs, follow the [Pathway Tools installation guide](pathway-tools.md)** to build and register your licensed SIF. Use `analysis_wf --skip_ptools` when you do not need pathway inference.

## Command boundaries

| Command | Work performed | Requires |
| --- | --- | --- |
| `analysis_wf` | One DAG for multiple assemblies, per-sample reads/MAG maps, PGDBs and combined reports | [Layout or manifest](analysis.md#analysis-wf-input-layout), MPDB, registered SIF |
| `build_db` | Download references, build FAST/BLAST indexes and supporting maps through Nextflow | Network, writable database directory, indexing tools |
| `run` | Assembly QC, prediction, search, annotation, optional mapping, Pathway Tools input preparation | Assembly and compatible MPDB; reads optional |
| `mag_split` | Reuse existing annotations to create MAG Pathway Tools inputs; preserve the original contig map | Completed community outputs, MAGSplitter, headerless two-column contig map |
| `build_pt` | Build, validate and register a private Pathway Tools SIF | Licensed installer, Apptainer, Nextflow, build prerequisites |
| `ptools` | Build/extract community and optional MAG PGDBs | Pathway Tools inputs, SIF/native installation, Camelot extraction dependency |
| `report` | Index existing outputs; optionally open a local search/export portal | Existing MP output; no biological tools needed for indexing |

Nextflow schedules tasks; existing MP scripts still perform the biological calculations and write the established sample paths. Public database construction uses the new Nextflow planner. The historical Snakefiles and old `bin/metapathways-data-install.sh` helper remain legacy artifacts; use `metapathways build_db` for supported construction. `--snakemake` on that command accepts only explicitly supported compatibility options (`cores`, `jobs`, `dryrun`, `dry-run`, `forceall`, `keep-going`, `rerun-incomplete`, `printshellcmds`, `latency-wait`); unsupported options fail rather than being interpreted as arbitrary Nextflow settings. Prefer the named resource flags.

## Annotation stages

The pipeline follows sample dependencies while permitting independent database searches/parses to run concurrently. Each sample maintains its own context and final output paths.

| Stage | Purpose | Representative output |
| --- | --- | --- |
| `PREPROCESS_INPUT` | Filter sequences; establish MP/original contig identifiers | `preprocessed/` |
| `ORF_PREDICTION` / `ORF_TO_AMINO` | Predict ORFs and derive sequences | `orf_prediction/` |
| `FILTER_AMINOS` | Apply amino-acid sequence filters | QCed protein sequences |
| `FUNC_SEARCH` | Search each selected functional database | `blast_results/` |
| `COMPUTE_REFSCORES` | Compute reference scores used by parsing | Refscore file |
| `PARSE_FUNC_SEARCH` | Apply search thresholds and parse database results | Parsed search tables |
| `SCAN_rRNA` | Barrnap prediction followed by reference search/summary | `results/rRNA/` |
| `SCAN_tRNA` | Predict tRNAs | `results/tRNA/` |
| `ANNOTATE_ORFS` | Assign functional annotations | Annotated GFF and annotation products |
| `CREATE_ANNOT_REPORTS` | Write primary functional/taxonomic and reference tables | `results/annotation_table/` |
| `GENBANK_FILE` | Export annotated sequence record | `genbank/*.gbk` |
| `PATHOLOGIC_INPUT` | Prepare selected Pathway Tools inputs and EC/reaction mappings | `ptools/`, `*.EC_RXN_map.tsv`, `*.ptinput.tsv` |
| `COMPUTE_TPM` | Map supplied reads and derive abundance/count tables | `results/rpkm/` |

FASTA-amino inputs take the compatible protein-input path. Stage choices and thresholds are in [CLI reference](cli-reference.md); internal derived stages do not necessarily have independent public flags. MP preserves the existing `yes`, `skip`, `redo` controls. Skipping a producer does not manufacture the inputs needed by its consumers.

QC defaults, ORF prediction parameters, alignment algorithm/mode, functional score/e-value/identity thresholds and rRNA database thresholds remain configurable through the CLI. Search flags must match the configured database indexes. The portal reports existing assignments; it does not change thresholding or infer taxonomic ranks from free-text lineage strings.

## Output layout

```text
OUTPUT/
├── SAMPLE/
│   ├── preprocessed/                 # Contigs and original-name map
│   ├── orf_prediction/               # Predicted and filtered sequences
│   ├── blast_results/                # Search and parsed hits
│   ├── genbank/                      # GFF/GenBank annotation
│   ├── ptools/                       # Community Pathway Tools inputs
│   ├── magsplitter/
│   │   ├── contig_to_mag.tsv          # Preserved original two-column map
│   │   └── results/MAG_ID/            # MAG-specific Pathway Tools inputs
│   ├── results/
│   │   ├── annotation_table/
│   │   ├── rpkm/
│   │   ├── rRNA/
│   │   ├── tRNA/
│   │   └── pgdb/
│   │       ├── community/
│   │       └── MAGs/MAG_ID/
│   ├── run_statistics/
│   ├── metapathways_steps_log.txt
│   └── errors_warnings_log.txt
├── reports/                          # Run report, EDA portal and SQLite
├── logs/                             # Durable CLI/Nextflow records
└── .metapathways/                     # Receipts, locks and temporary work
```

Commands operating on one existing sample (`mag_split`, `ptools`) also store their logs below that sample directory. A report over the parent includes these execution summaries. Running the report directly on a sample creates `SAMPLE/reports/`; running it on the parent creates `OUTPUT/reports/`.

## Dependency and cache behavior

Resource requests belong to each task. pProdigal, barrnap, ptRNAscan, FAST/BLAST, and read mapping/counting use the requested thread count capped by the total CPU budget. Serial stages, including Pathway Tools and Python parsers, reserve one CPU. Nextflow can mix threaded and serial tasks within the CPU and memory budgets. Independent database contexts share a preceding-stage dependency, and later stages wait for all contexts in that group. A failing required task stops dependent work. Optional MAG PGDB failures remain recorded while other entities proceed.

MP retains durable task receipts outside Nextflow work directories. Reuse compares command signatures and tracked file state (size and modification time, including tracked reference/index files). These are not full content hashes of all biological inputs. The report source manifest independently records SHA-256 for files it parses. Use external input/reference checksums for stronger provenance, and `redo` when changing tool implementations or environments that are not represented by command/file fingerprints.

Some legacy stages rewrite earlier outputs, such as the final GenBank record. The final producer owns that file's cached fingerprint. Old abundance output is not adopted without rerunning the corrected read mapping. Database index receipts track all index shards, so a missing shard invalidates reuse.

Temporary cleanup archives diagnostic files before removal. Explicit work/cache paths and interruptions are retained. Final results and logs are never considered disposable caches. Multiple controllers targeting the same command/output are blocked by a file lock; per-invocation resource budgets do not govern unrelated MP commands.

## Report limitations and error recovery

Report generation runs after successful workflow completion. It is disk-backed, single-controller postprocessing rather than a scheduled biological stage. Large multi-sample report builds still use headnode I/O/CPU; run standalone `report` where site policy permits metadata indexing. It does not alter biological outputs. For a failed workflow, generate a report explicitly over the partial output to inspect what exists.

A report requires recognized sample directories. Headers and identifiers matter: the importer rejects ambiguous duplicate source tables and malformed records. Do not merge different samples by concatenating tables containing unscoped `C1-G1` identifiers. Report schema versioning and source records make supported mappings explicit.

Missing source tables, unannotated placeholders and expected optional MAG failures are visible instead of silently becoming zeros. See [results schema](results-schema.md) for exact semantics and [reproducibility](reproducibility.md) for benchmark limitations.

For complete multi-sample execution, see [automatic input matching](analysis.md#analysis-wf-input-layout) and the [custom manifest](analysis.md#custom-analysis-manifest). `analysis_wf` assigns independent per-sample dependencies, saves the resolved input manifest, and builds the combined report after Nextflow completes so final statuses are included.

### RNA identifiers and abundance integrity

RNA annotations are emitted once per locus, including contigs with no CDS predictions. rRNA IDs include contig, subtype, coordinates and strand; distinct copies of the same rRNA must not share a gene ID. GFF-to-GTF conversion rejects duplicate IDs before mapping. ORF abundance uses featureCounts' `Length` column, requires a one-to-one gene-ID join with the GTF, and reports zero RPKM/TPM when all counts are zero. Tool output and errors are streamed into the task log, and failed sorting/conversion/counting/calculation commands stop the task immediately.

### Concurrent FAST searches

MP supplies FAST's `-X` option with a unique temporary directory for each search and reference-score invocation. FAST's default temporary names use a time-based seed and can collide when independent searches start in the same second. Older wrappers could therefore mix hits from different databases even when both processes returned success. A parser error showing Swiss-Prot `sp|...` targets in a MetaCyc result is one symptom; the databases themselves may be intact.

The corrected planner invalidates previous functional-search and reference-score receipts and does not adopt untracked outputs for these stages. Repeat the original run command to regenerate the searches and affected downstream results; unrelated unchanged preprocessing can be reused. Search volumes stop on the first tool error, and merged results are published only after all volumes succeed. This change does not require rebuilding MPDB or the Pathway Tools SIF.

### Pathway Tools failure diagnostics

Before each community or MAG PGDB invocation, MP expands its compact `0.pf` annotations into per-contig PathoLogic inputs in private staging. Each genetic element has an annotation file and a `SEQ-FILE` containing its real preprocessed contig sequence. Feature IDs are preserved. The sample's `ptinput.tsv` supplies each feature's original contig, coordinates and strand; these restore MAG member coordinates that the splitter copied from a representative annotation. Missing sequences, unmapped/duplicate features or invalid source coordinates stop preparation before Pathway Tools starts. Unannotated contigs without PGDB features are omitted. The sequence counts and coordinate-restoration counts are retained in `diagnostics/<invocation>/input/sequence-input.json`. This enables Pathway Tools to derive protein sequences for its PGDB BLAST databases. The source annotation files and MAG split outputs are not modified.

Each container PGDB attempt saves its internal `pathologic.log`, other available logs, reports and `execution.json` under the entity output's `diagnostics/ATTEMPT/` directory. The execution record identifies the image, exit code and phase (`build`, `export`, or `archive`). On failure MP also prints the last 80 lines of the internal Pathologic log. These diagnostics survive temporary container-state cleanup. An optional MAG task failing does not by itself establish that the failure is biological or expected; inspect its diagnostics. Community PGDB failures remain fatal.

Older container wrappers deleted the private Pathologic log along with temporary state, so the root cause of a historical exit 255 may be unavailable. Updating MP and repeating the same `ptools` invocation retains successful task receipts and retries unsuccessful entities, saving the internal error if it recurs. This logging fix works with an existing SIF. Newly built images also include NCBI BLAST+; the image recipe hash changes, so `build_pt` creates a separate image when rebuilt.

`build_pt` also snapshots official SRI patches for the selected release and validates a synthetic BLAST database/search before registration. Analysis tasks keep patch downloads disabled. Rebuilding creates a distinct image; changing the selected image can invalidate prior PGDB receipts. A patched image is not evidence that a particular vendor inference bug is resolved. Add `-d /path/to/MPDB -a fast` to prepare matching MetaCyc sequences, indexes and tables after the image build; see [MetaCyc from Pathway Tools](pgdb-workflow.md#metacyc-from-pathway-tools). This optional preparation operates on the reference, not on sample annotations or sample PGDBs.

Pathway Tools EC inputs are written as one `EC` line per identifier, including
when source annotations contain comma-, semicolon-, or pipe-separated lists.
PGDB staging applies the same normalization to existing PF inputs without
modifying the originals or rerunning annotation searches. The staging diagnostic
`sequence-input.json` records how many features had their EC entries normalized.
Provisional identifiers (for example `3.6.5.n1`) are preserved; Pathway Tools may
reject these, and its warnings remain in the saved logs. A successful Nextflow
wrapper for an optional MAG is not proof of a successful PGDB: consult MP's task
receipts and the per-entity `execution.json` status.

PGDB staging separates the amino-acid label and anticodon in recognized MP tRNA
names (for example, `C1094.tRNA2-GluTTC` becomes `C1094.tRNA2-Glu-TTC` in the
`NAME` field). This prevents the label's last letter from merging with the
three-base anticodon during name parsing. Feature `ID` fields, sequences,
coordinates and source files remain unchanged, so annotation and abundance
joins retain their original identifiers. The original and staged names are
recorded in `sequence-input.json` under `normalized_trna_names`. Unknown
anticodons such as `NNN` and unrecognized name formats are left unchanged;
MP does not guess their biological assignments. Updating this preparation
invalidates existing PGDB receipts, but does not require rebuilding the SIF
or rerunning annotations.

When the original `orf_prediction/<sample>.cds.gff` is available, PGDB staging
also preserves Prodigal's per-contig translation table through PathoLogic's
`CODON-TABLE` field. This avoids treating table-4 TGA codons as premature stops.
Older outputs without that metadata retain Pathway Tools' default genetic code;
no alternative code is guessed. The chosen codes are recorded in
`sequence-input.json`. This uses the input format documented in Pathway Tools'
installed `sample-genetic-elements.dat`; it does not patch Pathway Tools.

Concurrent SIF PGDB tasks use Xvfb with a private filesystem display socket.
Apptainer's `--containall` isolates `/tmp` but does not isolate Linux abstract
X sockets. MP disables Xvfb's abstract listener (`-nolisten local`) and creates
the private `/tmp/.X11-unix` directory. This avoids exhausting `xvfb-run`'s ten
display attempts when many containers start together, which can otherwise
return exit 1 during cleanup even after Pathway Tools saves its PGDB and prints
`Done`. Build and export X-server logs are retained as `build-xvfb.log` and
`export-xvfb.log` in the invocation diagnostics. This wrapper change requires
no SIF rebuild or vendor patch; nonzero build/export exits still fail the task.

For a controlled transport-inference comparison with an installed SIF, use
`metapathways ptools -o OUTPUT/SAMPLE --entity community --no_transport_inference`.
Omit `--entity` to process the community and MAGs. `analysis_wf` also accepts
`--no_transport_inference`. TIP remains enabled by default; changing the setting
invalidates the corresponding PGDB task cache. Failed SIF builds retain their
on-disk databases under `diagnostics/<run-id>/failed-pgdbs/`, together with the
staged inputs and original failure status. This preserves recovery material; it
does not automatically export or label partial databases as successful.

For an explicit PGDB taxon override, `ptools` and `analysis_wf` accept
`--taxon_id NCBI_ID`. This changes only private staged `organism-params.dat`
inputs, which are retained in diagnostics. It applies to every selected entity;
use `--entity community` with `ptools` for a community-only experiment.
`--taxprune --taxon_id 131567` uses cellular organisms with taxonomic pruning
enabled, avoiding the separate unpruned rescoring pass documented for
`-no-taxonomic-pruning` in the Pathway Tools User Guide (printed pp. 16–17).
This broad taxon includes multicellular eukaryotes too and is not a microbial-only
filter. Pruned and unpruned inference are different analysis settings.
Neither flag changes MP's gene taxonomic annotations.

`--taxonomic_scope all|bacteria|archaea|eukaryotes` provides readable aliases
for taxon IDs 131567, 2, 2157, and 2759, respectively; `euks` is accepted as an
alias for `eukaryotes`. It is mutually exclusive with `--taxon_id`. For example,
`--taxprune --taxonomic_scope all` is equivalent to
`--taxprune --taxon_id 131567`. No new default is imposed: omitting both retains
the input taxon and the existing pruning setting. A `prokaryotes` scope is not
yet supported because it needs a union of Bacteria and Archaea; MP rejects it
explicitly rather than silently substituting cellular life (which includes
eukaryotes). Scope changes guide PGDB inference and do not filter input contigs.

### If only `--threads 4` is specified

For the default local executor, capable tools receive up to four CPUs each (capped by the available CPU budget). MP detects available CPUs and memory for the total scheduling budget; four threads is not a four-CPU limit for the whole workflow. General tasks reserve the default 16 GB each. In `analysis_wf`, PGDB tasks inherit the same memory request and always reserve one CPU; `--ptools_memory` is an optional override. Nextflow schedules ready tasks together while their reservations fit. For example, a detected budget of 32 CPUs and 64 GB permits at most four simultaneous tasks that each request four CPUs and 16 GB, even though CPUs remain free. Dependencies can reduce concurrency further.

Detection is a startup snapshot, not exclusive ownership of the server. Set `--max_cpus` and `--max_memory` on shared machines. Local memory reservations govern scheduling, not hard memory enforcement. If the available memory budget is smaller than a task's reservation, planning fails; choose a suitable explicit `--memory` for a tiny test or provide more capacity. Slurm defaults to four submitted jobs and six submissions per minute; aggregate maxima apply only when explicitly supplied.

### Abundance checkpoints

`COMPUTE_TPM` checks the assembly, annotation GFF, supplied read files (including R2 for paired inputs), commands, and final abundance outputs when deciding whether to reuse results. Its `bwa/` directory contains generated intermediate files and is not an input dependency. Changes to those intermediates alone do not require recomputing valid abundance tables.

Checkpoints written before this correction included `bwa/` as an input. Their first invocation with the corrected code recomputes abundance once to replace that old checkpoint; subsequent unchanged runs reuse it. Other tasks retain their existing checkpoint rules.

### Large multi-sample workflows

MP writes `main.nf` plus small files under `modules/` in each invocation log directory. This keeps thousands of sample/bin tasks below the JVM compiled-method size limit. Modules expose individual task completion channels: they do not introduce a wait for the whole module, a separate resource budget, or serial sample execution. The shared executor limits and durable MP task checkpoints still apply. Keep the module files with `main.nf` when archiving a generated workflow.

If an older version failed at startup with `Method too large: Main.runScript`, update MP and repeat the same command with the same output directory. Compilation fails before any tasks are submitted; deleting the output is unnecessary. The next invocation generates a new modular workflow.
