# Software architecture: local and HPC

For the tool-by-tool diagrams and citations, see the [detailed workflow](detailed-workflow.md).

The MP command validates inputs and creates a dependency graph. Nextflow schedules ready tasks using the selected executor; the biological tools perform the same work in either mode. The MP controller and Nextflow continue running until the workflow finishes.

[![Architecture](assets/diagrams/architecture.svg)](assets/diagrams/architecture.svg)

[Zoom diagram](assets/diagrams/architecture.svg) · [Mermaid source](diagrams/architecture.mmd)

The two executor branches are alternatives selected by `--executor`; local execution is the default. Workers need access to the software environment, inputs, reference databases, image and output paths. On Slurm these must be accessible from the compute nodes, usually through shared storage. This diagram shows logical components; shared boxes do not mean jobs share a Pathway Tools instance.

Thread-capable tools use the requested `--threads`; serial tools, including Pathway Tools, request one CPU. Nextflow can run independent tasks simultaneously. Local concurrency fits detected or specified CPU/memory budgets. Slurm requests resources per job, with `--max_tasks` and `--submit_rate` bounding submissions. See [resource settings](resources.md).

Each containerized PGDB attempt has private Pathway Tools state, allowing separate entities to run concurrently. A community PGDB failure fails the required workflow task; MAG failures are retained as optional outcomes. Reports include task statuses and diagnostics. See [execution and restart behavior](execution.md) and [Pathway Tools](pathway-tools.md).

## Bounded controller graph

MP serializes task descriptors and dependency routing in `schedule.json` next to
`tasks.json`. Each dependency depth calls an alias of the same generic Nextflow
worker process with dynamic per-task CPU and memory requests. Completion tuples
are routed to consumers and grouped by a key carrying the expected prerequisite
count. Each consumer is emitted immediately when its own prerequisites arrive;
other tasks and samples at the same depth need not finish first. Process/operator
count grows with dependency depth rather than the number of MAG entities.
Biological commands and durable receipt checks remain in the existing worker.
