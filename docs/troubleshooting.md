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
