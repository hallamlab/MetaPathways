# Complete multi-sample analysis

**To build PGDBs, first complete the [Pathway Tools installation guide](pathway-tools.md)** and build/register your licensed SIF with `metapathways build_pt`. Once that setup is complete, run the workflow below. If you do not want pathway inference, add `--skip_ptools`; no Pathway Tools installation is needed.

`analysis_wf` runs annotation, optional read mapping, optional MAG splitting, community/MAG pathway inference, and a combined report. It accepts one or many metagenomes. All computational tasks share **one Nextflow DAG and one resource budget**; samples do not wait for other samples to finish. Community PGDB construction, read mapping, and MAG splitting can proceed together once their annotation inputs are ready. MAG PGDB jobs follow splitting.

## Analysis wf input layout

Use this exact layout for automatic discovery (directory and sample names are case-sensitive):

```text
inputs/
  assemblies/
    SampleA.fasta
    SampleB.fasta
  reads/
    SampleA_R1.fastq.gz
    SampleA_R2.fastq.gz
    SampleB_interleaved.fastq.gz
  mag_maps/
    SampleA.tsv
    SampleB.tsv
```

Assembly suffixes are `.fa`, `.fna`, or `.fasta`, optionally `.gz`. Read suffixes are `.fq` or `.fastq`, optionally `.gz`. Read names must end with `_R1` and `_R2` for paired files, `_interleaved` for interleaved pairs, or `_single` for single-end reads. Sample IDs use letters, digits and underscores, beginning with a letter. Names `logs`, `reports`, `inputs`, `assemblies`, `reads`, and `mag_maps` are reserved.

Every assembly needs exactly one read layout and one map by default. Each map is a **headerless, two-column TSV**: original assembly contig ID, then MAG ID. A contig may appear only once, and its ID must match the assembly FASTA header's first whitespace-delimited token. MAG IDs use letters, digits, underscores and periods, beginning with a letter. Periods become underscores in MAG output directory names; IDs that collide after that conversion are rejected. `community` and IDs containing `non_binned` are reserved.

```bash
metapathways analysis_wf \
  -i /path/to/inputs \
  -o /path/to/analysis \
  -d /path/to/MPDB \
  --threads 8 --max_cpus 32
```

The registered SIF from `metapathways build_pt` is used automatically; `--image /path/to/ptools.sif` overrides it. Container isolation allows multiple single-CPU Pathway Tools jobs. All workflow tasks inherit `--memory` (16 GB by default); `--ptools_memory` optionally overrides only PGDB jobs; memory availability can limit concurrency before CPUs do. Slurm uses the same [resource flags](resources.md) and requires all inputs, outputs, software and the SIF to be accessible on compute nodes.

For existing flat directories, point `-i` at the assemblies directory and supply `--reads_dir /path/to/reads` and `--mag_maps_dir /path/to/maps`. If only one of those flags is supplied, the other branch defaults to `reads/` or `mag_maps/` under `-i`. Use `--no_reads` or `--no_mags` to explicitly omit those branches for every sample during discovery. Use `--skip_ptools` to omit PGDB construction while retaining annotation, read mapping, MAG splitting and reporting. The manifest below supports a different combination for each sample.

Discovery does not recurse, merge sequencing lanes, or guess unmarked FASTQ layouts. Unexpected files/subdirectories, duplicate assembly names, orphan reads/maps, missing mates, reused input files, duplicate contig IDs and map IDs absent from their assembly stop the command **before any jobs are submitted**, with a link to this section. Hidden directory entries are ignored. Explicitly disabled branches are not scanned. Assembly headers and maps are checked fully; FASTQ files are checked for readability and nonzero size, not full sequencing integrity.

Append `--dryrun` to validate inputs and write the combined task plan without running biological tools. This still creates output directories, logs, staged symlinks and `inputs.resolved.tsv`. Remove `--dryrun` from the same command to execute. Input paths may contain spaces because MP stages symlinks; the output and MPDB paths must currently use letters, digits, underscores, hyphens, periods and slashes for compatibility with legacy tool commands.

## Custom analysis manifest

Use a tab-separated file with this **exact header and column order**:

```text
sample_id	assembly	read_layout	reads_1	reads_2	mag_map
SampleA	assemblies/assembly_a.fasta	paired	reads/a_1.fq.gz	reads/a_2.fq.gz	maps/a.tsv
SampleB	assemblies/assembly_b.fasta	interleaved	reads/b.fq.gz		maps/b.tsv
SampleC	assemblies/assembly_c.fasta	none	""	""	""
```

The example uses actual tabs, with `""` denoting empty fields in the final row. Paths are relative to the manifest's directory or absolute. `read_layout` is `paired`, `interleaved`, `single`, or `none`. Paired reads require both paths; single/interleaved reads require only `reads_1`; `none` requires both paths empty. A blank `mag_map` omits MAG splitting and MAG PGDBs for that sample. Community PGDB construction still runs unless `--skip_ptools` is supplied. Sample IDs come from the manifest, so original filenames do not need to match them. Assemblies must still have a supported FASTA suffix.

```bash
metapathways analysis_wf \
  --manifest /path/to/samples.tsv \
  -o /path/to/analysis -d /path/to/MPDB \
  --threads 8 --max_cpus 32
```

Do not combine `--manifest` with `-i`, discovery-directory flags, `--no_reads` or `--no_mags`. The single-sample `-1`, `-2`, `--interleaved`, `--samples` and `--test` options are not used by `analysis_wf`; each row declares its own reads and sample ID. Required annotation stages cannot be skipped; successful tasks are reused through receipts.

## Outputs, restarting and exploration

By default MP keeps all sample outputs. Add `--compact_results` to `analysis_wf` to remove large intermediates after each sample finishes while retaining rebuildable reports, final tables, logs and benchmark traces. PGDB archives are also removed. See [compact results and restart behavior](benchmarking.md#compact-results-on-limited-storage).

Each sample writes to `analysis/SAMPLE/`. The dataset root contains `inputs.resolved.tsv` with absolute original paths, combined `reports/`, and `logs/analysis_wf/RUN/` with the plan, tool logs, resource trace and task statuses. Staged symlinks and receipts stay under `.metapathways/analysis_wf/`; temporary Nextflow work/cache cleanup follows the ordinary MP resource options.

Repeat the same command/output to resume. The first invocation requires a new output directory. A changed sample list, MAG ID set, or input path mapping requires a new output directory; this prevents silently mixing datasets. Input contents, commands and resource requests are checked through task receipts, and completed outputs without a matching receipt are not adopted by this workflow. `--force_redo` reruns all selected annotation, splitting and PGDB tasks. Do not run another MP command into the same output while an analysis is active.

A community PGDB failure is a required-task failure. MAG PGDB failures are expected for some MAGs: they are logged as optional failures and do not stop the remaining MAGs. A MAG with no retained `0.pf` input after successful splitting is recorded as `SKIPPED`. The combined report retains per-sample/entity task statuses, including failures and skips.

After completion, open `analysis/reports/MP_run_report.html`, or start the searchable explorer:

```bash
metapathways report -o /path/to/analysis --serve --no-rebuild
```

The portal combines samples and keeps sample IDs on all related records so identically named MAGs and ORFs in different samples remain distinct.
