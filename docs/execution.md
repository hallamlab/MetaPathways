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

Controller output and task completion summaries are displayed and saved. Add `--stream_tool_output` to also mirror full worker output; `--verbose` alone does not enable it. Full tool output remains in individual task logs. Reports link the retained Nextflow records. Original sample logs remain available as well. Historical workflows without trace files cannot retrospectively gain CPU/RAM measurements by generating a report.

Pathway Tools containers place their private `/tmp` and `/var/tmp` on the disk
backing the task state directory, rather than Apptainer's limited in-memory
session filesystem. In compact mode this follows the selected task scratch
location (normally `SLURM_TMPDIR` on Slurm). Each invocation remains isolated;
this controller setting does not require rebuilding an existing SIF. Shared
JSON receipts and image-digest caches use unique temporary files and atomic
replacement, so concurrent publishers do not share a temporary filename.

Default Nextflow work, session cache and Conda package/environment caches live under `OUTPUT/.metapathways/COMMAND/tmp/RUN_ID`. After a successful invocation, diagnostics are archived and these temporary directories are removed. Failed and interrupted invocations retain work, caches and task logs because scheduler jobs may still be writing. The controller prints its failure and the first available Nextflow ERROR line before finalization. A later successful run also leaves earlier failed-run work in place. Confirm that the run’s jobs have stopped before manually removing retained directories.

Use `--keep_work` to retain the defaults. Explicit `--work_dir` and `--conda_cache` paths are retained automatically; MP does not delete a user-provided shared cache. Nextflow's installed runtime and MP's Conda environment are not per-run caches and are retained.

Durable receipts remain under `OUTPUT/.metapathways/COMMAND/receipts`. Repeating a command checks tracked inputs, outputs and command signatures, permitting reuse even after temporary cleanup. Never remove final outputs expecting receipts alone to recover them. To force one annotation stage, use its `redo` flag; `--force_redo` reruns all annotation stages. Old read-mapping output without a current receipt is deliberately recomputed once. FAST searches also require a matching success receipt: a standalone output file does not establish workflow completion. Use `--fast_local_db` within the workflow to preserve normal receipt-based reuse.

## Large Slurm runs

Run large controllers in a compute allocation where cluster policy permits job
submission from compute nodes. `--max_tasks` bounds submitted jobs; it does not
bound controller threads or memory. Use the {download}`controller sbatch template <../examples/slurm/analysis-controller.sbatch>`
with a script containing your environment activation and full `analysis_wf`
command. Submit with `sbatch --account=YOUR_ACCOUNT analysis-controller.sbatch run-analysis.sh`.
The template requests 4 CPUs and 32 GB for seven days; adjust to cluster policy
and observed usage. Its default `NXF_OPTS=-Xmx4g` limits Java heap only.

For file-count quotas, choose `--compact_results` before starting. It retains
featureCounts assignment summaries in `bwa/` along with final tables and PGDB
archives. Allow sufficient node-local scratch for concurrent builds. Existing
compact outputs require the same flag on restart; deleted intermediates and
previously deleted assignment summaries are not restored.

See [cluster troubleshooting](troubleshooting.md#slurm-controller-and-file-quota-failures)
for failure recovery.
