# Benchmarking and supplementary statistics

[Home](../README.md) · [Reproducibility records](reproducibility.md) · [Result schema](results-schema.md)

## Define the experiment before starting

Specify which assemblies, reads, genome assignments, reference releases, tool revisions and resource budgets are being evaluated. Use a complete input manifest and a new output directory for a fresh measurement. A reused task is a cache check, not a newly measured biological stage.

For CAMI, distinguish ground-truth source-genome bins from recovered MAGs. Ground-truth bins demonstrate annotation reuse, splitting, per-genome inference and report integration; they do not evaluate binning accuracy or robustness to contamination. Report the number of sample-specific bins separately from the number of distinct source genomes. Contig filtering and absence of selected annotated genes can reduce the number of actual PGDB attempts.

Freeze the MPDB and SIF during the run. Do not change analysis code, dependencies or reference files while tasks are running. Editing documentation does not change the biological workflow, but keep a source record identifying the code actually used. A Git commit alone is insufficient when the working tree has uncommitted changes.

## Example full workflow

```bash
metapathways analysis_wf --manifest /project/benchmark/samples.tsv \
  -o /project/benchmark/fresh-results -d /project/MPDB \
  --annotation_dbs swissprot metacyc \
  --taxprune --taxonomic_scope all \
  --threads 8 --max_cpus 32 --max_memory '64 GB'
```

The command uses the registered image. Add `--image /path/to/frozen.sif` to pin it explicitly. Include a maximum CPU and memory budget in the methods. Omitting them uses detected local availability, which can vary by host and cgroup. This is a throughput benchmark with concurrency, not necessarily an isolated per-sample latency measurement.

## What is recorded automatically

| Record | Location / meaning |
| --- | --- |
| Resolved inputs and sample IDs | Output-root `inputs.resolved.tsv` |
| Planned community/genome entities | Output-root `inputs.entities.json` |
| User command and console | `logs/cli/` and `logs/analysis_wf/RUN_ID/console.log` |
| Planned commands, CPU/memory reservations and dependencies | `tasks.json` in the invocation directory |
| Task status and measured wrapper elapsed time | `summary.json` and `tasks/*.json` |
| Nextflow timing and sampled resource fields | `trace.tsv`, `report.html`, `timeline.html` |
| Tool-specific counts and warnings | Per-task logs and sample `run_statistics/`, annotation/RNA/abundance outputs |
| PGDB build/export status and internal diagnostics | Per-entity `results/pgdb/.../diagnostics/` |
| Indexed biological tables and source accounting | `reports/results.sqlite`, `schema.json`, `output_inventory.tsv` |

Default Nextflow traces include status/exit, submit time, duration, realtime, CPU utilization, peak RSS/virtual memory and read/write character counters when available. Availability and units depend on the executor and trace format; retain the original headers and values. A missing measurement is unknown, not zero. `rchar`/`wchar` are process I/O counters, not a direct measurement of physical disk traffic.

Peak RSS is sampled task memory, not its requested reservation and not the peak of the whole concurrent workflow. Do not sum per-task peaks to claim a simultaneous machine peak. CPU utilization can exceed 100% for threaded work; it is not automatically normalized by the requested CPU count. Approximate CPU-seconds derived from utilization and elapsed time must be labeled as derived, not directly measured.

The standard trace does **not** provide a continuous whole-server CPU/RAM/disk time series, energy usage or complete network-traffic accounting. Those require a separate host monitor started with the benchmark. For the planned stage-runtime/resource figure, task-level trace measurements are appropriate; they do not support claims about an unmeasured whole-host peak.

## Supplementary table structure

Keep linked tables rather than forcing differently scoped quantities into one repeated wide table:

| Table | Recommended fields |
| --- | --- |
| Sample/input inventory | Sample ID, body site/group, assembly/read/map paths and hashes, read layout, original contig count/bases, input genome-bin count |
| Sample result statistics | Retained contigs/bases, predicted/annotated features with feature types, reference-hit counts, read mapping/counting statistics |
| Entity/pathway statistics | Sample/entity, ground-truth or recovered-bin method, selected-input genes, outcome, base pathways, reactions and explicit unique pathway–gene associations |
| Task performance | Invocation, sample, stage, database/entity where applicable, status, requested CPUs/memory, elapsed time, trace CPU and memory/I/O fields |
| Run provenance | Host/OS/CPU/RAM, MP/environment/reference/SIF versions and hashes, complete command, start/end, concurrency budget and notes |

Obtain each biological statistic from its authoritative source and define its unit. For example, an annotation-table row count is not necessarily the number of unique genes; a pathway-to-ORF export may contain multiple annotation rows per association. Count unique keys when claiming unique genes or links. The schema importer collapses some repeated associations deliberately.

Separate input bin count, generated PF input count, attempted PGDB count, successful PGDB count, failed count and skipped count. Do not convert a missing/failed inference into zero pathways. Distinguish base pathway counts from superpathways or the total records in `pathways.dat`.

## Runtime and resource figures

For stage comparisons, group tasks by stage and, where relevant, reference database or entity. Show distributions across samples/bins rather than hiding all variation in one mean. Report whether a panel includes fresh successes only, failed attempts, or both.

Parallel task times overlap. Summing their durations measures accumulated task elapsed time, not workflow makespan. Measure end-to-end wall time from controller start through final report completion. Nextflow's task trace does not include all validation/planning and post-workflow report-indexing time. MP's task elapsed time also includes wrapper overhead; it is not necessarily only the biological executable's compute time.

Keep separate measurements for MPDB construction and Pathway Tools image construction if they are discussed. Sample workflow traces do not retroactively capture setup performed before the run. Copying/downloading inputs is likewise outside sample-stage measurements unless explicitly included in the experiment.

## Save environment and machine information

Run these on the analysis host before the benchmark, adapting the checkout path. Save them beside the manifest, not in a temporary Nextflow work directory:

```bash
mkdir -p /project/benchmark/provenance
conda list --explicit > /project/benchmark/provenance/conda-explicit.txt
python -m pip freeze > /project/benchmark/provenance/pip-freeze.txt
metapathways version > /project/benchmark/provenance/mp-version.txt
nextflow -version > /project/benchmark/provenance/nextflow-version.txt
apptainer --version > /project/benchmark/provenance/apptainer-version.txt
lscpu > /project/benchmark/provenance/lscpu.txt
free -b > /project/benchmark/provenance/memory.txt
git -C /path/to/MetaPathways rev-parse HEAD > /project/benchmark/provenance/mp-commit.txt
git -C /path/to/MetaPathways status --short > /project/benchmark/provenance/mp-working-tree.txt
```

Also retain the source snapshot or checksums for uncommitted/untracked implementation files, the SIF and `.sif.json`, MetaCyc provenance, and reference release/checksum files. For very large raw datasets, compute checksums as a separate acquisition/validation step so that extra I/O does not distort measured analysis performance.

## Finish, inspect, then summarize

Read `summary.json` and per-entity `execution.json`, not just the top-level completion message. Inspect required task success, optional failures/skips, input coverage and report import notes. Preserve failed-attempt logs when retrying; a final successful attempt does not erase their computational cost.

Default cleanup removes disposable work/cache state after archiving diagnostics. It retains final results, receipts, logs, traces and Nextflow report/timeline. Explicit `--work_dir`, `--conda_cache` or `--keep_work` retain additional scratch; retaining it is not necessary merely to make the supplementary resource table. Interrupted runs may also retain scratch for diagnosis.

The final supplementary table and figure still require extraction, joining and quality checks after completion. The explorer provides biological subsets and an inventory; it does not automatically generate a publication-ready benchmark figure or certify manuscript statements. Audit every claim against the frozen final records.

## Compare one server with a Slurm cluster

Test the three small reviewer samples on the cluster before submitting the full benchmark. Install the same feature-branch revision on shared storage, and activate its environment before starting MP. Nextflow submits with the logged-in user's Slurm identity; there are no MP password flags. Confirm that `sbatch`, `squeue`, and `scancel` are available and that compute nodes can use the same software, database and input paths.

For a first three-sample cluster test, adapt the shared paths and allocation names:

```bash
metapathways analysis_wf \
  --manifest /shared/project/cami-reviewer/all.tsv \
  -o /shared/project/mp-reviewer-slurm \
  -d /shared/project/MPDB \
  --annotation_dbs swissprot \
  --skip_ptools \
  --threads 4 --max_cpus 8 --memory '4 GB' --max_memory '16 GB' \
  --executor slurm --account my_project --partition compute \
  --max_tasks 2 --submit_rate 6 --time_limit 2h
```

The example uses the full public SwissProt/SILVA references in your MPDB. If using the bundled test references, select `swissprot_test` and `SILVA_SSU_test SILVA_LSU_test` as in the reviewer walkthrough. The two-hour limit is a small-input example, not a full CAMI task limit. Follow your site's policy for keeping the Nextflow controller running: some sites permit a persistent headnode session, others require a controller allocation. MP preparation and report generation run in that controller, while scheduled biological stages run on compute nodes.

For the complete benchmark, copy the original assemblies, reads, CAMI genome maps, reference database and your licensed Pathway Tools SIF onto storage accessible to all selected nodes. Write an HPC-specific manifest pointing there; local server paths and symlinks are not portable. Use a new output directory. Supply the same analysis options, reference release, pruning scope, and image content as on the single server. Pass `--image /shared/project/pathway-tools.sif` explicitly if its user registration differs on the cluster.

Record a run/scenario ID (`local` or `slurm`), Git commit, environment export, input/reference/SIF checksums, task resource requests, node hardware, scheduler settings, start/end times, and log directory. Compare elapsed run time separately from summed task CPU time; cluster queue delays are part of operational elapsed time but not biological compute time. Do not mix a resumed run with a fresh timing run. If total budgets differ, report that as a throughput/scaling scenario rather than attributing the entire speed difference to Nextflow or Slurm. Keep all per-task traces and both scenario manifests for the supplementary tables.

For a larger concurrency test, the same per-job settings work on either executor: `--threads 8 --memory '64 GB' --max_tasks 100`. Omit aggregate maxima to avoid manual arithmetic. The local executor uses detected host capacity; Slurm permits up to 100 submitted jobs with no implicit aggregate caps. Add `--submit_rate 60` to allow up to one Slurm submission per second. Partition may be omitted to use the cluster default.

## Compact results on limited storage

Add `--compact_results` to `analysis_wf` on either local or Slurm execution. The default keeps all sample outputs. Compact mode is intended for the complete workflow; individual `run`, `mag_split`, and `ptools` commands do not perform this cleanup.

Each sample gets a final cleanup task that waits for **all** its annotation, read-abundance, splitting and requested PGDB tasks. Expected optional MAG failures count as finished attempts; their diagnostics and status remain available. Required failures prevent that sample's cleanup. Other successfully completed samples can already be compacted while the remaining samples run.

Retained files include the source tables and identifier maps used by the explorer, final supporting result tables, pathway TSVs, MAG input gene membership, ORF groups, run statistics and diagnostic logs. MP validates report relationships before deleting anything. Reports and CSV exports can be rebuilt normally with `metapathways report -o OUTPUT`.

Removed files include BAMs, sequence intermediates, raw alignment results, intermediate GenBank/GFF files, Pathway Tools working inputs except report-required membership records, and archived PGDBs (`*cyc.tar.bz2`). Compact output therefore cannot be used to reopen a complete PGDB in Pathway Tools. Original input assemblies/reads, the MPDB, and the SIF are untouched and must reside outside sample output directories.

Nextflow work and its run-local Conda cache are removed at normal controller exit. A successfully completed compact workflow also removes staging links and task-reuse receipts. Small controller locks, sample completion markers, input manifests, task plans, resource traces and logs remain as provenance. Active tasks still need temporary space: this flag reduces accumulated storage after sample completion, not the peak storage needed by concurrent unfinished samples. An interrupted controller retains temporary work until a successful resume, to avoid removing files potentially still in use by scheduler jobs.

`--compact_results` cannot be combined with `--keep_work`, `--work_dir` or `--conda_cache`; shared external caches are never deleted. It does not remove your installed Mamba environment or user-wide Nextflow installation.

Repeat the same compact workflow command to resume an interrupted run. Completed compact samples are retained and skipped; incomplete samples use normal checkpoints. An interrupted cleanup resumes from its marker. Changed settings, input/reference metadata, implementation, missing retained files, or `--force_redo` require a new output directory for compacted samples. Do not delete `compact-results.json` to try to restore checkpoint behavior: the intermediates have been intentionally removed.

For example, add this flag to the existing benchmark command:

```bash
metapathways analysis_wf \
  --manifest /project/benchmark/samples.tsv \
  -o /project/benchmark/compact-results -d /project/MPDB \
  --image /project/containers/pathway-tools.sif \
  --annotation_dbs swissprot metacyc \
  --taxprune --taxonomic_scope all \
  --threads 8 --compact_results
```

The cleanup task appears separately in the execution/resource records; distinguish its time from biological stages when preparing benchmark tables.
