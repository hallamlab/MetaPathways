# Command cookbook: run only the modules you need

[Home](../README.md) · [Complete CLI reference](cli-reference.md) · [Stage/output reference](workflow.md)

Use `metapathways COMMAND --help` for your installed revision. `metapathways version` reports the package version; also record the source Git revision. These examples assume an activated environment and writable output paths. Replace `/path/to/MPDB` and input filenames before running them.

## Choose a command

| Goal | Command | Meaning of `-o` |
| --- | --- | --- |
| Annotate one or more assemblies | `run` | Parent of sample directories |
| Annotate, map reads, split genomes, infer pathways and report | `analysis_wf` | Parent of sample directories |
| Prepare public/licensed reference indexes and tables | `build_db` | Use `-d` for the MPDB root |
| Build the licensed Pathway Tools image | `build_pt` | Directory holding generated SIF images |
| Split existing annotations by genome | `mag_split` | One existing sample directory |
| Infer community/genome pathways from existing annotations | `ptools` | One existing sample directory |
| Index/view existing results | `report` | One sample or a parent of samples |

## build_db: references before analysis

```bash
metapathways build_db -d /path/to/MPDB --func swissprot -a fast --dryrun
metapathways build_db -d /path/to/MPDB --func swissprot -a fast
```

The first command plans; the second downloads/formats. `--func` selects functional references; RNA, enzyme and taxonomy support are also prepared. Public options include SwissProt, CAZy and UniRef50/90. eggNOG requires the documented pre-supplied FASTA rather than an automatic acquisition rule. MetaCyc requires a licensed source; see [Pathway Tools](pathway-tools.md).

`-a fast` and `-a blast` choose the annotation index format. Match them with `run/analysis_wf --annotation_algorithm FAST` or `BLAST`. Build into a new directory when changing reference versions during an ongoing experiment. Database download/build time is separate from sample-analysis runtime.

On this command `-t` is the legacy aggregate CPU limit, unlike per-tool `run --threads`. Prefer explicit `--max_cpus` when communicating a shared resource budget. `--dryrun` performs no downloading/indexing. Historical Snakefiles are not the supported execution path.

## run: annotations and optional abundance

```bash
metapathways run -i SampleA.fasta -o results -d /path/to/MPDB \
  -1 SampleA_R1.fastq.gz -2 SampleA_R2.fastq.gz \
  --annotation_dbs swissprot metacyc \
  --threads 8 --max_cpus 16 --max_memory '64 GB'
```

This writes `results/SampleA/`. For interleaved pairs replace the two read arguments with `-1 SampleA_interleaved.fastq.gz --interleaved`. Single-end uses `-1` alone. No reads means no read abundance. MP passes the actual R2 to CoverM and uses the exact expected CoverM BAM for counting; a leftover sorted BAM is not a replacement for missing mapping output.

The command prepares Pathway Tools inputs but does not run Pathway Tools. For protein FASTA, standalone `run --input_format fasta-amino` uses compatible annotation stages; the complete nucleotide workflow and DNA/read/genome-splitting examples do not apply unchanged.

To inspect the plan add `--dryrun`. To rerun only read mapping/counting with the original other arguments, add `--COMPUTE_TPM redo`. Stage defaults, dependencies and outputs are documented in [workflow stages](workflow.md#annotation-stages). Do not use `--force_redo` merely to retry a failed downstream task: it forces annotation-stage recomputation.

## analysis_wf: one complete run for N samples

```bash
metapathways analysis_wf --manifest /project/samples.tsv \
  -o analysis -d /path/to/MPDB \
  --annotation_dbs swissprot metacyc \
  --taxprune --taxonomic_scope all \
  --threads 8 --max_cpus 32 --max_memory '64 GB'
```

Use [automatic directories or a manifest](inputs.md). All samples share one scheduling budget and dependency graph. Independent tasks can overlap, including mapping and PGDB work. `--skip_ptools` omits inference; discovery's `--no_reads` and `--no_mags` explicitly omit those inputs. In a manifest, blank optional input columns express per-sample omissions.

The workflow validates inputs and writes a resolved manifest before submitting analyses. Required annotation stages cannot be skipped through `skip`; valid completed tasks are reused automatically. `--ptools_memory '8 GB'` changes only each PGDB reservation. `--memory` changes all workflow task reservations, including PGDBs unless `--ptools_memory` is explicitly set. Increasing the total memory budget does not automatically increase either per-task setting.

## mag_split: reuse community annotation

```bash
metapathways mag_split -o results/SampleA -m /project/SampleA.tsv \
  --max_cpus 8 --max_memory '32 GB'
```

This needs the sample's completed annotation/PF files, original-contig mapping, authoritative feature table, and your headerless contig-to-genome map. It does not perform binning or redo reference searches. It produces MAG-specific Pathway Tools inputs and preserves the original map for full membership reporting.

Changing assignments can change per-genome pathway inference even if the community annotation is unchanged. Run `ptools` afterward, or let `analysis_wf` manage the dependency.

## build_pt and ptools: licensed setup, then inference

```bash
metapathways build_pt -i /path/to/pathway-tools-29.5-linux-64-tier1-install \
  -o ~/mp-containers -d /path/to/MPDB -a fast
metapathways ptools -o results/SampleA \
  --taxprune --taxonomic_scope all --max_cpus 8 --max_memory '32 GB'
```

The first builds/registers the image and optionally prepares MetaCyc annotation references. The second uses existing sample annotations to infer community and available MAG pathways. Use `--entity community` or `--entity MAG_001` to select a single entity. The [Pathway Tools chapter](pathway-tools.md) covers installer acquisition, image selection, patches, licensing, taxonomy, TIP, errors and outputs in detail.

## report: inspect and export without reanalysis

```bash
metapathways report -o results
metapathways report -o results --serve --no-rebuild
```

The first indexes current results. The second serves that snapshot for searching and CSV exports. Omit `--no-rebuild` when you want a refreshed index. Reports can be built from partial outputs, but cannot reconstruct missing resource measurements or repair old biological results. [Explorer instructions](reports-tutorial.md) include SSH access and safe table joins.

## Resources, failure and restart

Use the [resource and Slurm guide](../README.md#resources-and-slurm) for CPU/memory budgets and cluster flags. With `--threads 8 --max_cpus 32`, up to four ready eight-CPU tasks can fit, or a mixture of threaded and single-CPU tasks, subject to memory. Local execution is default; Slurm uses your logged-in identity, not credentials passed to MP.

Repeat the same command and output location after a recoverable failure. MP checks durable receipts and tracked files, reuses successful compatible tasks, and retries failed work. No explicit Nextflow `-resume` is required on the MP CLI. A new output directory is appropriate for a fresh benchmark, not for a simple restart.

Changing parameters, references, image, inputs or software can change reuse behavior. Output existence alone is not proof of provenance. Preserve `logs/` and `.metapathways/` receipts if you want to resume. See [restart details](../README.md#logs-temporary-files-and-restarting).
