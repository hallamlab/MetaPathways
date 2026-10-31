# Troubleshooting and reproducibility

- **A command fails before execution:** check the CLI/error logs, activate the environment and use `command -v nextflow coverm fastal` to check dependencies.
- **A task exceeds the memory budget:** lower concurrency or adjust `--memory`/`--max_memory` to measured requirements; Slurm can enforce its allocation.
- **Some MAG pathways are missing:** inspect the entity's last task status and logs. Expected MAG inference failure is distinct from missing annotation inputs.
- **The portal says local explorer required:** launch `metapathways report -o OUTPUT --serve --no-rebuild` and open the URL it prints.
- **Results were changed after report generation:** rebuild the report. It is an indexed snapshot, not a live view of files.
- **A report rebuild fails:** inspect the named source and expected schema. Fix/correct the source or report importer; do not invent missing identifiers or replace absent values with zero.

Record commands, Git commit, software versions, reference versions/checksums and input checksums for reproducible work. See [reproducibility notes](reproducibility.md) for the historical benchmark and the distinction between installation tests, synthetic software tests and full biological validation.

## Slurm controller and file quota failures

`unable to create native thread`, `pthread_create ... EAGAIN` and JVM allocation
failures can reflect PID/thread or memory limits on the controller node. Check
the printed preflight limits and `logs/analysis_wf/RUN/nextflow.log`. Actor thread
numbers in logs are sequence numbers, not simultaneous thread counts. Use the
[controller job guidance](execution.md#large-slurm-runs) for large datasets;
reducing `--max_tasks` does not directly cap JVM threads.

`Disk quota exceeded` can mean a file-count quota even when disk space remains.
Check both bytes and files with your site's quota tools. Unpacked community
PGDBs may create one or more files per contig. Choose `--compact_results` before
starting and provide sufficient worker-local scratch. `--file_quota` accepts a
remaining file allowance for an advisory estimate; MP does not change compaction
policy or delete partial PGDB directories automatically.

After a failed run, work and loose task logs remain in place. Stop or confirm
completion of that run's remaining Slurm jobs before deleting any retained paths;
MP does not cancel unrelated user jobs. Inspect `tasks.json`, `summary.json` and
the failed task's log to locate partial results. When failed PGDB output directories exist,
`partial_pgdb_candidates.json` lists candidates and observed task statuses for
inspection; it does not certify them as safe to delete. Keep validated final outputs and
durable receipts when resuming. Do not remove compact completion markers or edit
the resolved manifest to bypass provenance checks. Existing compacted runs require
fresh output after an implementation upgrade; missing featureCounts summaries
from older compaction cannot be restored automatically.


## FAST searches stall on a shared filesystem

FAST can spend most of its walltime waiting for random index reads on Lustre.
A running process with very low CPU use is not necessarily making useful progress.
Use `--fast_local_db` with `analysis_wf` to copy indexes sequentially to worker-local
storage before searching. In the default `--annotation_run_mode pervol`, each
volume is copied, searched and removed before the next volume. In `default` mode,
the complete database must fit. Each concurrent search has its own private copy.
Scratch must accommodate the largest volume, temporary hit files and merged output;
concurrent tasks multiply this requirement.

The scratch directory is resolved on the compute worker from `SLURM_TMPDIR`.
On clusters without it, export `METAPATHWAYS_FAST_SCRATCH` to a node-local path
available on every worker. This setting is independent of compact-result scratch.
Staging fails explicitly when indexes or scratch space are unavailable; it does
not silently fall back to searching the shared filesystem. Normal success receipts
are written, so subsequent workflow invocations can reuse the results.

For example, add `--fast_local_db --func_search_time_limit 2h --time_limit 48h`
to the existing Slurm command. The shorter limit applies to FUNC_SEARCH tasks;
other tasks retain the general limit. Select a search limit appropriate to your
inputs and cluster. Automatic retries and adaptive stall detection are not enabled.
A timeout remains a failed task with retained diagnostics.

The standalone `python -m metapathways.MetaPathways_func_search` command accepts
`--stage_db` and returns a nonzero exit status on failure. A standalone output is
not a workflow success receipt: copying it into `blast_results` does not guarantee
`ALREADY_COMPUTED`. Failed or mismatched receipts cause the workflow to rerun the
search. Prefer staging within the managed workflow; do not edit receipts to mark
unverified outputs successful. Existing compacted results remain subject to the
upgrade restrictions described in the release notes.


### Test only affected samples after upgrading from 4.0.1

`--samples` selects rows from the original full manifest; keep that manifest
unchanged. The remaining samples are validated but no computational tasks are
scheduled for them. Use sample IDs exactly as written in the `sample_id` column.

Previously compacted 4.0.1 samples cannot resume under 4.0.3 because completion
markers are implementation-specific and intermediates have been removed. Use
new output directories to test just the failed sample/assembler combinations.
This reruns the selected samples' workflows, including other annotation databases
and downstream steps; it does not rerun the other samples or merge the test
results back into their original reports. Do not delete old outputs or receipts.
Noncompacted results can be reused where success receipts still match, but an
upgrade may invalidate some tasks; reuse is not guaranteed.

For the reported SPAdes `S245_S55` and MEGAHIT `S567_S65` failures, stop old
controllers/jobs before updating the software environment, then use the commands
below. Replace paths, retain any custom reference/Pathway Tools settings from
the original invocation, and ensure the Pathway Tools SIF is registered in the
updated environment (or supply the original `--image /absolute/path/to/ptools.sif`).

```bash
# Stop the old controllers/jobs before upgrading their environment.
mamba install -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.3
conda activate metapathways
metapathways version

# Replace these paths. Keep the ORIGINAL full manifests unchanged.
DB=/absolute/path/to/MPDB
SPADES_MANIFEST=/absolute/path/to/manifest_spades.tsv
MEGAHIT_MANIFEST=/absolute/path/to/manifest_megahit.tsv
TEST_OUT=/absolute/path/to/new/fast-local-test-4.0.3

# Run one controller at a time, from a host permitted to submit Slurm jobs.
metapathways analysis_wf \
  --manifest "$SPADES_MANIFEST" --samples S245_S55 \
  -o "$TEST_OUT/spades" -d "$DB" \
  --annotation_dbs cazy swissprot uniref50 uniref90 metacyc \
  --threads 32 --executor slurm --account def-shallam \
  --memory 128GB --ptools_memory 16GB --max_tasks 4 \
  --fast_local_db --func_search_time_limit 2h --time_limit 48h \
  --compact_results --keep_work

# After the SPAdes invocation finishes, run the MEGAHIT case.
metapathways analysis_wf \
  --manifest "$MEGAHIT_MANIFEST" --samples S567_S65 \
  -o "$TEST_OUT/megahit" -d "$DB" \
  --annotation_dbs cazy swissprot uniref50 uniref90 metacyc \
  --threads 32 --executor slurm --account def-shallam \
  --memory 128GB --ptools_memory 16GB --max_tasks 4 \
  --fast_local_db --func_search_time_limit 2h --time_limit 48h \
  --compact_results --keep_work
```

Fir provides `SLURM_TMPDIR` to each compute job; do not substitute the controller's
scratch path. `--max_tasks 4` conservatively limits concurrent scratch consumers
for this test. The two-hour search limit is a starting point, not a guarantee for
all databases. `--keep_work` retains Nextflow diagnostics after success and uses
additional storage; `--compact_results` still removes selected sample
intermediates at completion. Full worker logs remain retained or archived.

If the run fails, share the printed error, invocation `summary.json`, `trace.tsv`
and the affected task log under `OUTPUT/logs/analysis_wf/RUN_ID/`. Those records
can distinguish staging, search, scheduler and downstream failures.
