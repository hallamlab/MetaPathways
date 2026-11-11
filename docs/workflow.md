# Workflow and output guide

For the tool-by-tool diagrams and citations, see the [detailed workflow](detailed-workflow.md).

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
| `ptools` | Build/extract community and optional MAG PGDBs | Pathway Tools inputs, registered or explicit SIF, Camelot extraction dependency |
| `report` | Index existing outputs; optionally open a local search/export portal | Existing MP output; no biological tools needed for indexing |

Nextflow schedules the analysis. Use the named resource flags (`--max_cpus`, `--max_tasks`, `--executor`) and `--dryrun`; use `build_db --update` to replace selected references in an existing MPDB.

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

Repeat the same command with the same output directory to reuse completed,
compatible work and retry unfinished tasks. Keep the original inputs, references,
results and logs. See [restarting](execution.md) before forcing a rerun.

## Report limitations and error recovery

Reports show available results and task outcomes. Missing outputs and optional
MAG failures are not zero-valued biological results. You can run `metapathways
report -o OUTPUT` over partial results to inspect what completed.

### Pathway Tools failure diagnostics

Keep `diagnostics/ATTEMPT/execution.json`, `input/pathologic.log` and the available
build/export logs under the failed entity's PGDB output. Compact runs retain them
in diagnostic archives. Community failures stop the workflow; optional MAG
failures remain recorded while other entities may finish.
See the [Pathway Tools troubleshooting table](pathway-tools.md#outputs-warnings-and-failures).

### Concurrent FAST searches

If a search fails, inspect its task log and rerun the original command after
fixing the cause. Do not copy standalone search outputs into a workflow to mark
it complete. See [FAST troubleshooting](troubleshooting.md#fast-searches-stall-on-a-shared-filesystem).

### RNA identifiers and abundance integrity

Keep sample and feature identifiers when joining RNA, annotation and abundance
tables. Use the [results schema](results-schema.md) for supported joins and units.

### If only `--threads 4` is specified

Each threaded task may use up to four CPUs; this is not a four-CPU limit for the
whole workflow. Set `--max_cpus` and `--max_memory` to control total local usage.
See [resources and Slurm](resources.md).

### Abundance checkpoints

Repeat the original workflow command to reuse compatible abundance outputs.
Missing required outputs cause that work to rerun. Do not delete completion
markers or edit receipts to bypass validation.

### Large multi-sample workflows

Follow the [controller and scratch guidance](execution.md#large-slurm-runs).
Keep the run logs when reporting controller failures.
