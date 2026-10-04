# Resources and Slurm

## Local execution is the default

| Setting | Default / meaning |
| --- | --- |
| `run -t / --threads` | Eight CPUs per threaded tool, capped by total CPU budget |
| Threaded tools | pProdigal, barrnap, ptRNAscan, FAST, BLAST, CoverM, featureCounts, and samtools sort |
| Serial stages | One CPU each, including Pathway Tools and Python parsers; native math pools are capped too |
| `--max_cpus` | CPUs available to the process, including affinity/cgroup limits |
| `--memory` | Per-task reservation: normally 16 GB; standalone `ptools` defaults to 4 GB |
| `analysis_wf --ptools_memory` | Optional PGDB-only override; otherwise inherits `--memory` |
| `--max_memory` | Currently available host memory, bounded by cgroups |
| `--max_tasks` | Additional concurrency cap; CPU/memory reservations still apply |
| `build_db -t` | Legacy total CPU budget; omitted means available CPUs |
| `build_pt -t` | Image-compression CPUs |

For example, permit a total of 16 CPUs while allowing eight threads per capable tool:

```bash
metapathways run -i sample.fasta -o results -d /data/MPDB \
  -t 8 --max_cpus 16 --max_memory '64 GB'
```

Independent reference searches may run together, and Nextflow schedules as many eligible tasks as fit their reservations. Memory reservations are scheduling requests, not measured peaks; local execution does not impose a container memory ceiling. Each MP invocation has its own budget. `--threads` controls each threaded task, while `--max_cpus` controls the total: for example, a 16-CPU budget can fit two eight-CPU jobs, one eight-CPU job plus eight single-CPU jobs, or sixteen single-CPU jobs when dependencies and memory permit. ptRNAscan uses that many single-threaded tRNAscan-SE workers; samtools receives one fewer additional thread so its main thread fits the reservation. Thread counts are upper bounds: small inputs and serial portions of a tool may use fewer CPUs. These settings apply when planning a new invocation; they do not resize running jobs.

## Submit from a Slurm headnode

MP explicitly requests one node (`--nodes=1`) for every Slurm task. Threaded tools use their allocated CPUs on that node; concurrency across nodes comes from separate jobs. No node-count flag is needed in the MP command.

Run from an authenticated login/head node where `sbatch`, `squeue` and `scancel` are available:

```bash
metapathways run -i /shared/sample.fasta -o /shared/results -d /shared/MPDB \
  --executor slurm --account my_project --partition compute \
  -t 8 --memory '64 GB' --max_tasks 100 --time_limit 24h
```

MP uses your existing Slurm identity; it does not take a password or SSH private key. `--partition` is optional: omitted means the cluster default. Run `sinfo` to list partitions; these are named node groups/queues with access and time limits. `--qos` and `--reservation` are optional. Inputs, outputs, references, work/cache paths, the MP installation and its environment must have the same absolute paths on compute nodes. Pathway Tools on Slurm requires a SIF. The report's HTML explorer runs locally; it is not a Slurm service.

Slurm defaults to at most four submitted jobs, six submissions per minute and 24 hours per task. `--max_tasks 100` directly permits up to 100 queued/running jobs; Slurm decides actual running concurrency. There are no implicit aggregate CPU or memory caps on Slurm. If you explicitly provide `--max_cpus` or `--max_memory`, MP conservatively reduces the job cap using the largest task request. `--submit_rate` independently throttles submissions; for example, `--submit_rate 60` permits one per second. There are no automatic task retries. MP uses an explicit Nextflow configuration so ambient profiles do not silently override these limits.

The basic interface is the same locally and on HPC: `--threads 8 --memory '64 GB' --max_tasks 100`. Serial tasks still request one CPU. In `analysis_wf`, the memory request also applies to PGDBs unless overridden with `--ptools_memory`. Locally, detected available CPUs and memory automatically limit execution; on Slurm, the scheduler determines capacity. Aggregate maxima are optional additional controls, not numbers the user must calculate.

Live Slurm execution still needs validation at your site. Configuration and synthetic local scheduling tests do not establish compatibility with every cluster's authentication, filesystem or resource policies.
