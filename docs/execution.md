# Logs, temporary files and restarting

Every command with an output location saves its CLI transcript under `OUTPUT/logs/cli/`. Workflow invocations also retain:

```text
OUTPUT/logs/COMMAND/RUN_ID/
├── console.log             # Controller and labeled tool output
├── tasks.json              # Commands, dependencies, resources and paths
├── tasks/                  # Individual output logs and status records
├── summary.json            # Invocation outcome and per-task results
├── nextflow.log
├── nextflow_tasks/          # Archived .command.* diagnostics
├── trace.tsv
├── report.html
├── timeline.html
├── main.nf
└── nextflow.config
```

Console output is displayed and saved. Reports link the retained Nextflow records. Original sample logs remain available as well. Historical workflows without trace files cannot retrospectively gain CPU/RAM measurements by generating a report.

Pathway Tools containers place their private `/tmp` and `/var/tmp` on the disk
backing the task state directory, rather than Apptainer's limited in-memory
session filesystem. In compact mode this follows the selected task scratch
location (normally `SLURM_TMPDIR` on Slurm). Each invocation remains isolated;
this controller setting does not require rebuilding an existing SIF. Shared
JSON receipts and image-digest caches use unique temporary files and atomic
replacement, so concurrent publishers do not share a temporary filename.

Default Nextflow work, session cache and Conda package/environment caches live under `OUTPUT/.metapathways/COMMAND/tmp/RUN_ID`. After a completed invocation, including an ordinary task failure, diagnostics are archived and these temporary directories are removed. Interrupted invocations retain them because cancellation may still be in progress.

Use `--keep_work` to retain the defaults. Explicit `--work_dir` and `--conda_cache` paths are retained automatically; MP does not delete a user-provided shared cache. Nextflow's installed runtime and MP's Conda environment are not per-run caches and are retained.

Durable receipts remain under `OUTPUT/.metapathways/COMMAND/receipts`. Repeating a command checks tracked inputs, outputs and command signatures, permitting reuse even after temporary cleanup. Never remove final outputs expecting receipts alone to recover them. To force one annotation stage, use its `redo` flag; `--force_redo` reruns all annotation stages. Old read-mapping output without a current receipt is deliberately recomputed once.
