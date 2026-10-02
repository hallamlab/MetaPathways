# MetaPathways

## Abstract

The development of high-throughput sequencing technologies over the past decade has generated a tidal wave of environmental sequence information from a variety of natural and human engineered ecosystems. The resulting flood of information into public databases and archived sequencing projects has exponentially expanded computational resource requirements rendering most local homology-based search methods inefficient. MetaPathways v1.0 is a modular annotation and analysis pipeline for constructing environmental Pathway/Genome Databases (ePGDBs) from environmental sequence information capable of using the Sun Grid engine for external resource partitioning. However, a command-line interface and facile task management introduced user activation barriers with concomitant decrease in fault tolerance.

MetaPathways has since advanced as a modular tool, deepening our understanding of microbial metabolism at various biological levels. With this release, we have addressed previous challenges in modularity and database management. v3.5 enhances user accessibility through streamlined installation via package indexes or containers, refined modules, and interface upgrades. It boasts updated algorithm support for sequence feature prediction, annotation, metabolic inference, and coverage metrics. Tested on mock community data, Metapathways v3.5 demonstrates improved performance and usability. With automated installation and database management, this open-source tool makes advanced metagenomic analysis more accessible. Metapathways v3.5 represents a significant step forward in automated, comprehensive metagenomic analysis, facilitating a deeper exploration of microbial interactions and metabolic functions in environmental genomics.

## Quick start

Install MetaPathways, try the included three-sample dataset, then use your own data. Linux x86-64 is supported. Choose one installation method below; **Mamba is recommended**.

### 1. Conda package with Mamba (preferred)

```bash
mamba create -n metapathways --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=3.5.2
conda activate metapathways
```

This installs MP and its workflow dependencies, including MAGSplitter and Camelot. Prepare the included data, build the small reference database, and run all three samples:

```bash
metapathways prepare_test -o ~/mp-reviewer
cd ~/mp-reviewer
metapathways build_db --test -d MPDB
metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools --threads 4 --memory '4 GB' --max_tasks 2
metapathways report -o all --serve --no-browser --port 8765
```

Open the URL printed by the report server. On a remote server, use an [SSH tunnel](docs/reports-tutorial.md#view-a-remote-report-through-ssh). The **2.4 MiB input dataset is included** in the package; database preparation downloads enzyme and taxonomy support records. The test covers annotation, paired-read abundance, genome splitting, reports and exploration. Pathway inference is skipped because it requires your own Pathway Tools license.

### 2. Quay: Docker or Apptainer

```bash
# Docker
docker pull quay.io/hallamlab/metapathways:3.5.2

# Or Apptainer
apptainer pull metapathways.sif docker://quay.io/hallamlab/metapathways:3.5.2
```

The image includes the same workflow dependencies and reviewer data. Follow the [Docker three-sample test](docker/README.quay.md#docker-three-sample-test) or [Apptainer three-sample test](docker/README.quay.md#apptainer-three-sample-test) to run the commands with your working directory mounted for persistent results. Licensed Pathway Tools is a separate image.

### 3. Local installation from GitHub

```bash
git clone https://github.com/hallamlab/MetaPathways.git
cd MetaPathways
mamba env create -f docker/conda_base.yml
conda activate metapathways
mamba install --yes -c conda-forge pip
python -m pip install .
```

Then run the **same three-sample commands under option 1**, starting with `metapathways prepare_test -o ~/mp-reviewer`. MP installs its Python workflow helpers automatically. The data comes from the installed package; the test does not depend on your checkout location.

### Try your own data

Build a production reference database, then annotate an assembly:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways run -i /path/to/assembly.fasta -o results -d ~/MPDB --threads 8
```

For assemblies with reads and genome maps, follow the [complete workflow](docs/inputs.md). **If you want PGDBs, complete the [Pathway Tools installation guide](docs/pathway-tools.md) before starting that workflow.** The small reviewer references are for testing only.

## Start here

New to the terminal? Follow the chapters in order. Already installed MP? Start with the reviewer walkthrough or choose a command below. All documentation is maintained here on GitHub; a separate Read the Docs site is not required.

1. **[Getting started](docs/getting-started.md):** terminal basics, installation, activation and your first result.
2. **[Reviewer walkthrough](docs/reviewer-test.md):** three tiny CAMI samples for single-sample and two-sample CLI tests, including reads, genome splitting and reports; Pathway Tools is optional.
3. **[Organize your inputs](docs/inputs.md):** automatic matching, symbolic links, custom sample IDs/manifests and genome maps.
4. **[Command cookbook](docs/commands.md):** each module's purpose, required inputs, commands and outputs.
5. **[Pathway Tools](docs/pathway-tools.md):** where to obtain the licensed installer, one-command image builds, MetaCyc reference preparation, taxonomy, transport inference and troubleshooting.
6. **[Explore and export](docs/reports-tutorial.md):** reports, SSH access, related tables and filtered CSVs.
7. **[Resources and Slurm](#resources-and-slurm):** per-tool threads, total budgets, cluster submission and restart behavior.
8. **[Benchmarking](docs/benchmarking.md):** resource measurements, supplementary tables, figures and interpretation limits.
9. **[Maintainer release guide](docs/releasing.md):** PR testing, Conda/Quay publication and Zenodo setup; [current readiness audit](docs/release-readiness.md).

### Which command do I need?

| Goal | Command |
| --- | --- |
| Complete workflow for one or many metagenomes | `metapathways analysis_wf` |
| Annotate assemblies, optionally map reads | `metapathways run` |
| Prepare the included reviewer data | `metapathways prepare_test` |
| Prepare reference databases | `metapathways build_db` |
| Build/register a licensed Pathway Tools SIF | `metapathways build_pt` |
| Split existing community annotations into genome bins | `metapathways mag_split` |
| Infer pathways from prepared community/MAG inputs | `metapathways ptools` |
| Browse existing outputs and export subsets | `metapathways report` |

### Reference chapters

| Guide | Contents |
| --- | --- |
| [CLI reference](docs/cli-reference.md) | Every command, flag and default, generated from current CLI help |
| [Result schema](docs/results-schema.md) | Table units, keys, joins, missing data and SQLite examples |
| [Workflow and output reference](docs/workflow.md) | Stages, directories, caching and detailed diagnostics |
| [Reproducibility](docs/reproducibility.md) | Software/reference records, historical provenance and validation limits |
| [Containers](docker/README.quay.md) | Released MP Docker/Apptainer images; distinct from your licensed Pathway Tools SIF |
| [Maintainer releases](docs/releasing.md) | Packaging, release checks and publication |

The remainder of this README is a compact working reference. The linked chapters provide the step-by-step explanations.

## Installation details and reviewer test

The [getting-started guide](docs/getting-started.md) explains installation choices, activation, troubleshooting and basic terminal use. The [reviewer walkthrough](docs/reviewer-test.md) explains expected outputs, single- and two-sample variants, automatic discovery and optional licensed pathway inference.

The standard installation check is the three-sample workflow above. The legacy `run --test` K12 example remains available for compatibility, but is not the complete reviewer test.

## Prepare reference databases

For a minimal reference database:

```bash
metapathways build_db -d /data/MPDB --func swissprot -a fast
```

SILVA and supporting enzyme/taxonomy files are built alongside functional references. Select multiple functional databases with `--func swissprot cazy uniref50`; `uniref90` is also supported and substantially larger. `-a blast` builds BLAST indexes instead of FAST indexes. Use `--dryrun` to inspect the planned downloads and indexing jobs first.

The builder downloads public references from their configured upstream locations. These are not all version-pinned. eggNOG requires a local FASTA at `MPDB/functional/eggnog`; the old builder had no working eggNOG acquisition rule. Licensed MetaCyc reference acquisition is separate from this public builder. Existing compatible MPDB installations can be used directly with `-d`.

Choose annotation databases actually present in your MPDB. The FAST/BLAST choice in `run` must match its indexes. Keep reference release records and checksums; changing a database can change biological results.

## Annotate your own data

### One assembly

```bash
metapathways run \
  -i sample.fasta \
  -o results \
  -d /data/MPDB
```

The sample output is `results/sample/`; the run-wide report is `results/reports/`. Input can be nucleotide FASTA, compressed FASTA, or amino-acid FASTA with `--input_format fasta-amino`. GFF and GenBank are output formats, not accepted primary inputs of this CLI.

Choose references and inspect the plan:

```bash
metapathways run -i sample.fasta -o results -d /data/MPDB \
  --annotation_dbs swissprot cazy --annotation_algorithm FAST --dryrun
```

A dry run checks dependencies and writes the execution plan and logs without executing annotation tasks. It can create sample directory scaffolding; it does not produce a completed biological report.

### Read mapping

Paired reads in separate files:

```bash
metapathways run -i sample.fasta -o results -d /data/MPDB \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz
```

Interleaved paired reads:

```bash
metapathways run -i sample.fasta -o results -d /data/MPDB \
  -1 sample.interleaved.fastq.gz --interleaved
```

Single-end reads use `-1` alone. With no reads, abundance calculation is skipped. `-1` and `-2` must identify different files; MP now rejects identical paired input paths, including aliases to the same file. The paired command passes the actual reverse file to CoverM.

Historical output generated by the incorrect mate selection must not be treated as corrected merely because it appears in a report. Repeating the annotation command with the corrected reads can reuse completed annotation stages and recompute mapping. The portal copies existing abundance values; it does not certify their provenance or rerun mapping.

### Multiple assemblies

Pass a directory of FASTA files to `-i`; each file becomes a sample with its own output directory. Use unique, simple filenames and `--samples` to select a subset. Read mapping arguments apply to one sample only: run separate commands for samples with different FASTQs. CPU/memory budgets are per invocation, so independent MP commands do not share a global resource limit.

### Stage controls

Stage flags accept `yes`, `skip` or `redo`. For example:

```bash
metapathways run -i sample.fasta -o results -d /data/MPDB \
  --COMPUTE_TPM redo -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz
```

`yes` reuses valid outputs, `redo` executes the stage again, and `skip` requires any needed downstream inputs to already exist. `--force_redo` forces all annotation stages. Changes to tracked inputs or outputs invalidate task receipts. See the [workflow guide](docs/workflow.md) before skipping dependencies.

## MAGs and pathway inference

### Build Pathway Tools once

Follow the [Pathway Tools licensing and installer guide](docs/pathway-tools.md#get-the-installer) to obtain your licensed Linux x86-64 installer, then provide it to MP:

```bash
metapathways build_pt \
  -i ./pathway-tools-29.5-linux-64-tier1-install \
  -o ./containers
```

Nextflow runs the Apptainer build, validates Pathway Tools startup, and registers the finished SIF in `~/.config/metapathways/ptools.json`. Subsequent `ptools` commands use that image automatically. Without `-o`, images go under `~/.local/share/metapathways/containers`; the XDG equivalents are honored. `--image` on `ptools` or `METAPATHWAYS_PTOOLS_IMAGE` overrides registration.

Apptainer must support unprivileged builds on your host. Where subordinate UID/GID mappings are configured, its `newuidmap` and `newgidmap` helpers must be installed by the system administrator (Ubuntu provides them in `uidmap`). Building needs internet access for the base image and operating-system packages. Docker is not required. The installer itself is supplied by you; MP does not obtain a license.

`build_pt -t` controls compression CPUs, default two. If shared storage makes package installation slow, set `APPTAINER_TMPDIR` to sufficiently large local scratch. Startup validation has been exercised with 29.5; do not assume every installer release has the same unattended interface. If the installer was renamed, supply `--ptools_version 29.5` explicitly.

Every explicit `build_pt` invocation fetches the installer release's official Linux patches from SRI over HTTPS, following the [vendor's patch installation instructions](https://bioinformatics.ai.sri.com/ptools/faq.html). There is no prompt and no custom Pathway Tools code or MetaCyc data patch. Failed patch downloads stop the build. The build loads these patches before freezing the SIF; analysis invocations disable further patch downloads. Each invocation creates a separate image, preserving existing images. The adjacent `.sif.json` records installer/image/recipe hashes, patch source URLs and individual checksums, snapshot hash, and startup validation output; the patch manifest is also embedded in the image. Keep this private image to reproduce its exact patch set.

New images include NCBI BLAST+ and its configuration. Before registration, validation creates and searches a tiny synthetic protein database, then checks Pathway Tools startup. This BLAST installation supports Pathway Tools' own sequence databases and hole filling; it is separate from the FAST/BLAST annotation indexes in MPDB. Installing official patches or BLAST does not establish that any particular Pathway Tools inference failure has been fixed.

### MetaCyc from Pathway Tools

To build the SIF and prepare its bundled MetaCyc reference in an existing MPDB with one command:

```bash
metapathways build_pt \
  -i ./pathway-tools-29.5-linux-64-tier1-install \
  -o ./containers -d /path/to/MPDB -a fast
```

Use `-a blast` for BLAST annotation indexes instead. Nextflow first builds and validates the SIF, then exports its bundled MetaCyc flat files alongside `protseq.fsa` in private staging. It runs the repository's `metacyc_mapping_build.py` and `metacyc_build_ont.py`, checks the reference files and tables, indexes the proteins, and installs:

| Location relative to MPDB | Contents |
| --- | --- |
| `functional/metacyc` | Protein FASTA from the selected MetaCyc release |
| `functional/formatted/metacyc.*` | Selected FAST or BLAST indexes |
| `functional/formatted/metacyc-names.txt` | Protein identifiers and descriptions |
| `functional_categories/MetaCyc-monomer-rxn-pairs.tsv` | Protein-to-reaction mappings, including containing complexes |
| `functional_categories/MetaCyc-PWY-RXN-CMP-map.tsv` | Pathway, reaction, enzyme and primary compound mappings |
| `functional_categories/MetaCyc_PWY_Ontology.tsv` | Pathway ontology |
| `functional_categories/MetaCyc_reldate.txt`, `MetaCyc_provenance.json` | Release, source/image hashes, input checksums and mapping counts |

The MetaCyc preparation task uses one CPU. Other MPDB references are not refreshed. Preparation replaces the existing MetaCyc reference and removes obsolete MetaCyc index files, including indexes for the other aligner, so run it when no analyses are reading that MPDB. Annotation runs must use the selected index format. Proteins without reaction mappings are counted in provenance; absence of a reaction assignment is permitted.

To add MetaCyc later using the registered SIF:

```bash
metapathways build_db -d /path/to/MPDB --func metacyc -a fast
```

Select a particular SIF with `--metacyc_source /path/to/pathway-tools.sif`, or provide a complete licensed flat-file `data/` directory. That directory must contain `protseq.fsa`, `proteins.dat`, `enzrxns.dat`, `reactions.dat`, `pathways.dat`, `compounds.dat`, and `classes.dat`. A FASTA alone cannot supply the reaction and pathway mappings. Native installations may contain only sequence files in `data/` because their remaining reference data are built into the executable; the SIF export handles this. Files supplied over SFTP must first be copied locally.

MetaCyc is opt-in and is not downloaded by the ordinary default `build_db` command. The user supplies the licensed installation/data and remains responsible for its permitted use; MP does not grant a license or publish the installer, reference data, patches, or SIF. `--dryrun` plans either route without exporting references, fetching patches, indexing, or registering an image.

### Community pathways

```bash
metapathways ptools -o results/sample --taxprune --taxonomic_scope all
```

These examples use broad cellular-life taxonomic pruning. Choose the scope to match your analysis and record it in your methods; see [scope choices](docs/pathway-tools.md#choose-a-taxonomic-scope).

Each SIF task gets private Pathway Tools data, home and temporary state, allowing concurrent isolated instances. Pathway Tools uses one CPU per PGDB. A community PGDB failure fails the command. Native Pathway Tools is still supported when no image is selected, but serialized to protect shared state. The legacy `--container` flag retains its original meaning and bypasses automatic SIF selection.

### Add MAG pathways without reannotating

Create a tab-separated map with **no header**, containing original assembly contig identifiers and MAG identifiers:

```text
original_contig_1	MAG_001
original_contig_2	MAG_001
original_contig_3	MAG_002
```

Then run:

```bash
metapathways mag_split -o results/sample -m contig_to_mag.tsv
metapathways ptools -o results/sample --taxprune --taxonomic_scope all
```

MAG splitting reuses the community annotation and contig mapping. It also preserves the supplied map at `magsplitter/contig_to_mag.tsv`, allowing the report to connect all retained MAG contigs and their ORFs. MAG Pathway Tools inputs are a smaller, selected set and are shown separately.

Pathway Tools can fail for individual MAGs. These failures are recorded and allowed to remain optional; they do not fail an otherwise successful community workflow. Missing or failed inference is not equivalent to “zero pathways.” Inspect `entities` and execution history in the portal. `--taxprune` enables taxonomic pruning; the default leaves it off.

## Reports and the EDA portal

Successful `run`, `mag_split` and `ptools` commands refresh the report automatically. When a sample belongs to an existing parent report, the parent report is refreshed so all samples remain accessible. Report indexing runs on the controller after Nextflow completes, uses SQLite on disk, and does not repeat annotation, mapping or pathway inference.

For existing or partial results, build a report independently:

```bash
metapathways report -o results
```

`-o` can be one sample output directory or a parent containing sample directories as immediate children. Reports are saved under:

```text
results/reports/
├── MP_run_report.html       # Accounting, import notes, files and execution links
├── EDA_portal.html          # Search and export interface
├── results.sqlite          # Relational data for all included samples
├── schema.json             # Columns, row counts, keys and relationships
├── output_inventory.tsv    # Paths, sizes and modification times
├── portal.js
└── portal.css
```

Open `MP_run_report.html` directly for navigation. To search, subset and export:

```bash
metapathways report -o results --serve --no-rebuild
```

This opens a browser and serves the existing snapshot on loopback only. Omit `--no-rebuild` to refresh changed results first. `--no-browser` prints the URL; `--port 8765` chooses a fixed local port. Ctrl-C stops the server. All assets are local, and the explorer makes no third-party requests. Copy the complete output directory to a workstation to retain file links; no cloud upload is involved.

Protein taxonomy supports SwissProt (including the small reviewer database), UniRef and eggNOG. Each annotation row reports the taxon of its own reference hit in `taxonomy`; `lca_taxonomy` is computed from score-qualified hits in that same database, with support counted independently for each database. The primary ORF row follows its selected annotation's `reference_db` and target. No database priority or cross-database fallback is used. These describe reference-hit evidence, not a definitive organism assignment for the query. Unsupported databases report `Not computed`; missing or unknown taxon IDs in supported databases report `Unclassified`. An actual LCA of `root` remains `root`. SILVA rRNA taxonomy is separate. Existing outputs need annotation-table regeneration to gain these fields.

The workflow in the portal is:

1. Choose a table: ORFs/taxonomy, functional annotations, pathways, pathway genes, MAG membership, abundance or file inventory.
2. Search text, add column filters, or filter ORF-based tables by linked EC/reaction, reference database and pathway/entity identifiers.
3. Follow a row's related-result buttons. Sample, pathway and entity keys remain attached, avoiding cross-sample collisions.
4. Select export columns, apply filters, and export **all matching rows** as CSV. Pagination does not truncate an export.
5. Save the query definition or bookmark its URL to record your selection.

The portal is a table explorer: it adds no biological plots or interpretation. The [schema guide](docs/results-schema.md) explains row units, relation keys, safe joins and the distinction between missing output and biological absence. Large joins remain in SQLite; CSV is streamed rather than assembled entirely in browser memory. Broad scans can still take time. Each query has a two-minute execution budget, and at most four queries execute concurrently.

The full inventory links remaining RNA, sequence, alignment, statistics and raw PGDB products. The relational views cover the recognized formats documented in the schema guide; an arbitrary legacy table is not silently assumed to fit that schema. Missing primary annotations leave explicit placeholder ORFs and import notes; malformed recognized tables stop the rebuild and preserve the previous database.

## Complete multi-sample analysis

**To build PGDBs, first complete the [Pathway Tools installation guide](docs/pathway-tools.md)** and build/register your licensed SIF with `metapathways build_pt`. Once that setup is complete, run the workflow below. If you do not want pathway inference, add `--skip_ptools`; no Pathway Tools installation is needed.

`analysis_wf` runs annotation, optional read mapping, optional MAG splitting, community/MAG pathway inference, and a combined report. It accepts one or many metagenomes. All computational tasks share **one Nextflow DAG and one resource budget**; samples do not wait for other samples to finish. Community PGDB construction, read mapping, and MAG splitting can proceed together once their annotation inputs are ready. MAG PGDB jobs follow splitting.

### Analysis wf input layout

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

The registered SIF from `metapathways build_pt` is used automatically; `--image /path/to/ptools.sif` overrides it. Container isolation allows multiple single-CPU Pathway Tools jobs. All workflow tasks inherit `--memory` (16 GB by default); `--ptools_memory` optionally overrides only PGDB jobs; memory availability can limit concurrency before CPUs do. Slurm uses the same resource flags described below and requires all inputs, outputs, software and the SIF to be accessible on compute nodes.

For existing flat directories, point `-i` at the assemblies directory and supply `--reads_dir /path/to/reads` and `--mag_maps_dir /path/to/maps`. If only one of those flags is supplied, the other branch defaults to `reads/` or `mag_maps/` under `-i`. Use `--no_reads` or `--no_mags` to explicitly omit those branches for every sample during discovery. Use `--skip_ptools` to omit PGDB construction while retaining annotation, read mapping, MAG splitting and reporting. The manifest below supports a different combination for each sample.

Discovery does not recurse, merge sequencing lanes, or guess unmarked FASTQ layouts. Unexpected files/subdirectories, duplicate assembly names, orphan reads/maps, missing mates, reused input files, duplicate contig IDs and map IDs absent from their assembly stop the command **before any jobs are submitted**, with a link to this section. Hidden directory entries are ignored. Explicitly disabled branches are not scanned. Assembly headers and maps are checked fully; FASTQ files are checked for readability and nonzero size, not full sequencing integrity.

Append `--dryrun` to validate inputs and write the combined task plan without running biological tools. This still creates output directories, logs, staged symlinks and `inputs.resolved.tsv`. Remove `--dryrun` from the same command to execute. Input paths may contain spaces because MP stages symlinks; the output and MPDB paths must currently use letters, digits, underscores, hyphens, periods and slashes for compatibility with legacy tool commands.

### Custom analysis manifest

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

### Outputs, restarting and exploration

By default MP keeps all sample outputs. Add `--compact_results` to `analysis_wf` to remove large intermediates after each sample finishes while retaining rebuildable reports, final tables, logs and benchmark traces. PGDB archives are also removed. See [compact results and restart behavior](docs/benchmarking.md#compact-results-on-limited-storage).

Each sample writes to `analysis/SAMPLE/`. The dataset root contains `inputs.resolved.tsv` with absolute original paths, combined `reports/`, and `logs/analysis_wf/RUN/` with the plan, tool logs, resource trace and task statuses. Staged symlinks and receipts stay under `.metapathways/analysis_wf/`; temporary Nextflow work/cache cleanup follows the ordinary MP resource options.

Repeat the same command/output to resume. The first invocation requires a new output directory. A changed sample list, MAG ID set, or input path mapping requires a new output directory; this prevents silently mixing datasets. Input contents, commands and resource requests are checked through task receipts, and completed outputs without a matching receipt are not adopted by this workflow. `--force_redo` reruns all selected annotation, splitting and PGDB tasks. Do not run another MP command into the same output while an analysis is active.

A community PGDB failure is a required-task failure. MAG PGDB failures are expected for some MAGs: they are logged as optional failures and do not stop the remaining MAGs. A MAG with no retained `0.pf` input after successful splitting is recorded as `SKIPPED`. The combined report retains per-sample/entity task statuses, including failures and skips.

After completion, open `analysis/reports/MP_run_report.html`, or start the searchable explorer:

```bash
metapathways report -o /path/to/analysis --serve --no-rebuild
```

The portal combines samples and keeps sample IDs on all related records so identically named MAGs and ORFs in different samples remain distinct.

## Resources and Slurm

### Local execution is the default

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

### Submit from a Slurm headnode

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

## Logs, temporary files and restarting

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

Default Nextflow work, session cache and Conda package/environment caches live under `OUTPUT/.metapathways/COMMAND/tmp/RUN_ID`. After a completed invocation, including an ordinary task failure, diagnostics are archived and these temporary directories are removed. Interrupted invocations retain them because cancellation may still be in progress.

Use `--keep_work` to retain the defaults. Explicit `--work_dir` and `--conda_cache` paths are retained automatically; MP does not delete a user-provided shared cache. Nextflow's installed runtime and MP's Conda environment are not per-run caches and are retained.

Durable receipts remain under `OUTPUT/.metapathways/COMMAND/receipts`. Repeating a command checks tracked inputs, outputs and command signatures, permitting reuse even after temporary cleanup. Never remove final outputs expecting receipts alone to recover them. To force one annotation stage, use its `redo` flag; `--force_redo` reruns all annotation stages. Old read-mapping output without a current receipt is deliberately recomputed once.

## Troubleshooting and reproducibility

- **A command fails before execution:** check the CLI/error logs, activate the environment and use `command -v nextflow coverm fastal` to check dependencies.
- **A task exceeds the memory budget:** lower concurrency or adjust `--memory`/`--max_memory` to measured requirements; Slurm can enforce its allocation.
- **Some MAG pathways are missing:** inspect the entity's last task status and logs. Expected MAG inference failure is distinct from missing annotation inputs.
- **The portal says local explorer required:** launch `metapathways report -o OUTPUT --serve --no-rebuild` and open the URL it prints.
- **Results were changed after report generation:** rebuild the report. It is an indexed snapshot, not a live view of files.
- **An old output has no complete MAG map:** preserve its original two-column map as `SAMPLE/magsplitter/contig_to_mag.tsv` before rebuilding. The identifiers must match that sample's original contig names.
- **A report rebuild fails:** inspect the named source and expected schema. Fix/correct the source or report importer; do not invent missing identifiers or replace absent values with zero.

Record commands, Git commit, software versions, reference versions/checksums and input checksums for reproducible work. See [reproducibility notes](docs/reproducibility.md) for the historical benchmark and the distinction between installation tests, synthetic software tests and full biological validation.

## Team, support and citation

Current team: Ryan J. McLaughlin, Tony X. Liu, Tomer Altman, Aditi N. Nallan, Aria S. Hahn, Julia Anstett, Connor Morgan-Lang, Kishori M. Konwar and Steven J. Hallam. Previous contributors include Niels W. Hanson and Shang-Ju Wu.

Source: [hallamlab/MetaPathways](https://github.com/hallamlab/MetaPathways). Historical code: [MetaPathways-legacy](https://github.com/hallamlab/MetaPathways-legacy). Questions and bug reports: [GitHub issues](https://github.com/hallamlab/MetaPathways/issues). License: [MIT](LICENSE), with bundled third-party license notices retained.

Please cite:

> McLaughlin RJ, Liu TX, Altman T, Nallan AN, Hahn AS, Anstett J, Morgan-Lang C, Konwar KM, Hallam SJ. *MetaPathways v3.5: Modularity and Scalability Improvements for Pathway Inference from Environmental Genomes*. bioRxiv (2024). [doi:10.1101/2024.06.04.597460](https://doi.org/10.1101/2024.06.04.597460).

Maintainers: use the [release controller and CI workflow](docs/releasing.md). The source version is declared in `metapathways/_version.py`. No release or benchmark is launched by editing this documentation.
