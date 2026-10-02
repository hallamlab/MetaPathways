# Troubleshooting and reproducibility

- **A command fails before execution:** check the CLI/error logs, activate the environment and use `command -v nextflow coverm fastal` to check dependencies.
- **A task exceeds the memory budget:** lower concurrency or adjust `--memory`/`--max_memory` to measured requirements; Slurm can enforce its allocation.
- **Some MAG pathways are missing:** inspect the entity's last task status and logs. Expected MAG inference failure is distinct from missing annotation inputs.
- **The portal says local explorer required:** launch `metapathways report -o OUTPUT --serve --no-rebuild` and open the URL it prints.
- **Results were changed after report generation:** rebuild the report. It is an indexed snapshot, not a live view of files.
- **A report rebuild fails:** inspect the named source and expected schema. Fix/correct the source or report importer; do not invent missing identifiers or replace absent values with zero.

Record commands, Git commit, software versions, reference versions/checksums and input checksums for reproducible work. See [reproducibility notes](reproducibility.md) for the historical benchmark and the distinction between installation tests, synthetic software tests and full biological validation.
