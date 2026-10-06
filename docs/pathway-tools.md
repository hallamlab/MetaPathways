# Pathway Tools: licensing, image builds, databases and inference

[Home](index.md) · [Test walkthrough](test.md) · [Complete flags](cli-reference.md)

## Defaults and explicit choices

| Setting | MP default | Your alternative |
| --- | --- | --- |
| Reaction compatibility screening when building MetaCyc | Enabled; successful full screens automatically publish an MPDB compatibility list | `--skip_pt_screen` on `build_pt` or `build_db` |
| Taxonomic pruning for PGDB inference | Enabled; avoids the unpruned rescoring pass | `--no_taxprune` |
| Organism taxonomic scope | `all`: NCBI 131567, cellular life; broad enough for mixed-domain inputs | `--taxonomic_scope bacteria`, `archaea`, `eukaryotes`, or `--taxon_id ID` |
| Transport inference in SIF runs | Enabled | `--no_transport_inference` |
| CPUs per PGDB | One; concurrency uses independent PGDBs | Control concurrency through the workflow resource flags |

These defaults apply to `ptools` and `analysis_wf`. Taxonomic scope guides PGDB
inference without changing the gene-level taxonomy in MP annotation tables.
The compatibility list removes only unsafe **explicit reaction assignments**,
not genes, annotations or sequences; the same reaction may still be inferred
later. The sections below explain the evidence, limits and diagnostic records.

## Understand the three different databases

| Name | What it contains | How you use it |
| --- | --- | --- |
| MPDB | Reference protein/RNA sequences, search indexes, taxonomy and mapping tables | `-d /path/to/MPDB` on annotation/workflow commands |
| MetaCyc in Pathway Tools | Reference biochemical knowledge used for pathway inference | Bundled in an appropriate licensed Pathway Tools distribution |
| Your PGDB | Inferred pathways/reactions and gene associations for a community or genome bin | Created by `ptools` or `analysis_wf` under the sample output |

Adding MetaCyc to **MPDB** enables another protein annotation search. Pathway Tools can infer pathways from SwissProt-derived annotations without that search: it still consults its own bundled MetaCyc knowledge base. Therefore “without MetaCyc annotation” does not mean “without MetaCyc pathway knowledge.”

EcoCyc describes E. coli; MetaCyc is a multi-organism reference; BioCyc is the collection of organism PGDBs. Compare sample pathway counts with the appropriate kind of database and counting unit. Base pathways, superpathways, and signaling pathways are not interchangeable totals.

## Get the installer

Start at the official [Pathway Tools site](https://www.pathwaytools.org/) and its licensing/download link. SRI provides separate academic and commercial licensing routes. The [academic license page](https://bioinformatics.ai.sri.com/ptools/licensing/ptools-academic-license.shtml) describes eligibility and the request process; SRI sends download instructions to the approved technical contact. Use the terms supplied by SRI for your organization.

Download the **Linux x86-64 installer**, even if your browser is on a Mac or Windows laptop: MP builds a Linux container on the analysis server. Choose a licensed distribution containing MetaCyc; the vendor describes its editions in the [installation guide](https://www.pathwaytools.org/installation-guide/released/index.html). Retain the installer filename, version, download records and checksum. Download addresses provided in license correspondence need not be public.

If it is on your laptop, copy it to the server using your normal file-transfer method. For example, in a local terminal, replacing the username and host:

```bash
scp ~/Downloads/pathway-tools-29.5-linux-64-tier1-install user@server:~/Downloads/
```

Create the destination directory on the server first if needed. Do not unpack or manually run the installer before giving it to `build_pt`. MP performs the unattended installation within the build. Version 29.5 is a tested example, not a claim that it is the newest version. Consult the vendor's [release notes](https://bioinformatics.ai.sri.com/ptools/release-notes.html) for newer releases and verify their installer compatibility separately.

## Host prerequisites

Activate the MP environment and check `apptainer --version` and `nextflow -version`. Building requires network access for the container base, operating-system packages and official patches, plus writable space for the installer, temporary build and finished SIF. A build can temporarily use substantially more disk than the final image. Use `df -h` to inspect free space.

The host must permit Apptainer's unprivileged/fakeroot build path. Administrators may need to configure subordinate UID/GID ranges and the `newuidmap`/`newgidmap` helpers (`uidmap` on Ubuntu). Installing a Conda executable cannot override a host security policy. Ask your administrator if a user-namespace or fakeroot error appears. Docker is not required.

Use shared storage for a finished image that compute nodes need to read. Local scratch can accelerate building on slow network filesystems. `APPTAINER_TMPDIR` is an optional Apptainer troubleshooting setting, not required to run MP. See your site's Apptainer setup when choosing scratch; it must have enough space.

## Build and register once

Run on the server, with the environment active:

```bash
metapathways build_pt \
  -i ~/Downloads/pathway-tools-29.5-linux-64-tier1-install \
  -o ~/mp-containers
```

MP uses Nextflow to build the image, download official release-specific SRI patches, validate the installation, and register the resulting SIF. Validation checks Pathway Tools startup and a small BLAST database/search. Wait for the build command to complete successfully before launching dependent analyses.

`-o` is the image directory. Without it, the default is `~/.local/share/metapathways/containers` (or its XDG equivalent). Registration is stored at `~/.config/metapathways/ptools.json` (or its XDG equivalent). The adjacent `.sif.json` records provenance, including hashes and patch records. `--ptools_version 29.5` supplies a version when MP cannot recognize a renamed installer. `-t 2` controls compression CPUs, not the number of CPUs used by later PGDB tasks.

An explicit new build creates a new image and preserves older images. Automatic patch downloads are disabled during analyses, so runs use the patch snapshot built into that image. There are no MP-authored edits to licensed Pathway Tools code or MetaCyc frames. A successful image validation does not prove that every possible genome will infer successfully.

Subsequent `ptools` and `analysis_wf` commands use the registered SIF automatically. No image environment variable is necessary. To pin a particular image for a comparison, add:

```text
--image /absolute/path/to/pathway-tools-VERSION-HASH.sif
```

Explicit `--image` takes precedence over the environment override and registration. Do not delete a SIF referenced by a running task. Rebuilding an image or changing the selected path can invalidate PGDB task reuse. On Slurm the selected path must be readable at the same location on all compute nodes.

## Build the MPDB MetaCyc annotation reference at the same time

Prepare public references first, then build the licensed image and add its MetaCyc reference:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
metapathways build_pt \
  -i ~/Downloads/pathway-tools-29.5-linux-64-tier1-install \
  -o ~/mp-containers -d ~/MPDB -a fast
```

The second command exports the bundled MetaCyc protein sequences and biochemical flat files from the SIF, constructs protein-to-reaction mappings and pathway/compound/ontology tables, and indexes the proteins for FAST. It does not refresh SwissProt or the other existing references. The [installed-file table](pgdb-workflow.md#metacyc-from-pathway-tools) lists the exact MPDB destinations.

To add or rebuild MetaCyc later using the registered image:

```bash
metapathways build_db -d ~/MPDB --func metacyc -a fast
```

To choose a different licensed source:

```bash
metapathways build_db -d ~/MPDB --func metacyc -a fast \
  --metacyc_source /path/to/pathway-tools.sif
```

`--metacyc_source` also accepts a local complete exported MetaCyc `data/` directory. It needs `protseq.fsa`, `proteins.dat`, `enzrxns.dat`, `reactions.dat`, `pathways.dat`, `compounds.dat`, and `classes.dat`. A `protseq.fsa` alone is insufficient even if its sequence statistics match another FASTA. Some native installs keep most knowledge-base data in the executable; the SIF export obtains the missing flat files.

Use `-a blast` if subsequent annotation will use `--annotation_algorithm BLAST`. MetaCyc rebuilding replaces its existing reference and removes obsolete index files, including the other aligner's indexes. **Do not rebuild a reference directory while analyses are reading it.** Use a separate MPDB for a concurrent reference-version comparison.

The ordinary public-reference build does not acquire licensed MetaCyc. You supply the authorized installer or data. MP does not grant permission to redistribute the installer, SIF, patches or reference files; preserve them according to the applicable vendor terms.

## Run pathway inference

For an annotated sample at `results/SampleA`:

```bash
metapathways ptools -o results/SampleA \
  --entity community --taxprune --taxonomic_scope all \
  --max_cpus 8 --max_memory '32 GB'
```

For all available entities, omit `--entity`:

```bash
metapathways mag_split -o results/SampleA -m /path/to/SampleA.tsv
metapathways ptools -o results/SampleA \
  --taxprune --taxonomic_scope all \
  --max_cpus 8 --max_memory '32 GB'
```

`--entity MAG_001` selects one actual normalized MAG ID. There is no `--entity MAGS` keyword. Omitting the flag includes the community and available genome bins; successful unchanged entities can be reused.

`-o` here is the **sample directory**, unlike `run` and `analysis_wf`, whose `-o` is the parent output directory. A wrong level can make MP look for missing `ptools/0.pf` inputs.

Every PGDB task uses one CPU. Concurrency comes from multiple independent PGDB tasks; assigning eight threads does not make one Pathway Tools process eight times faster. `ptools --memory` controls each standalone PGDB reservation; `analysis_wf --ptools_memory` optionally overrides PGDB reservations within the complete workflow; without it PGDBs inherit `--memory`. An explicit `--max_memory` additionally bounds scheduled reservations. Local reservations are not hard memory ceilings.

SIF tasks receive private home, data and temporary state. MP also isolates X-display sockets; independent PGDBs can run together without the native installation's shared state. Native Pathway Tools remains a legacy serialized option for standalone `ptools`; the complete `analysis_wf` requires the SIF route. The old `--container` flag is not a synonym for `--image`; use the image route shown here.

## Choose a taxonomic scope

`--taxprune` enables taxonomic pruning. `--taxonomic_scope` selects the organism taxon assigned to the staged PGDB inputs:

| Scope | NCBI taxon | Use |
| --- | ---: | --- |
| `all` | 131567 | Broad cellular-life scope for a mixed community |
| `bacteria` | 2 | A bacterial community or bacterial genome |
| `archaea` | 2157 | An archaeal community or archaeal genome |
| `eukaryotes` (`euks`) | 2759 | Eukaryotic inputs |

`all` includes multicellular eukaryotes; it is not a microbial-only or unicellular filter. There is no supported `prokaryotes` union scope. `--taxon_id NCBI_ID` allows a more specific positive numeric taxon and is mutually exclusive with the convenience scope flag.

The scope applies to all selected entities in that invocation. It guides pathway inference; it does not remove contigs or rewrite MP's gene-level taxonomic annotations. To use different scopes for different MAGs, run separate entity-specific `ptools` commands.

**Defaults: taxonomic pruning enabled, scope `all` (cellular life).** You can omit both flags for that behavior. Use `--no_taxprune` to disable pruning, `--taxonomic_scope bacteria`, `archaea` or `eukaryotes` to narrow the scope, or `--taxon_id` for a specific taxon. Pruning constrains inference using the chosen taxon; it does not disable other pathway-selection rules. Record these settings in your methods. Pathway Tools 29.5 can fail during its unpruned rescore pass; MP retains diagnostics rather than treating partial output as successful.

See the vendor User Guide's batch PathoLogic discussion and [MP's diagnostic detail](workflow.md#pathway-tools-failure-diagnostics). The installed PDF is authoritative for the installed version.

## Explicit reaction blacklist

For SIF runs, MP applies a matching MPDB compatibility list to private staged community and MAG PGDB inputs. Its bundled fallback, `metapathways/resources/ptools_reaction_blacklist.json`, is restricted to the exact SIF tested for those entries. Legacy native runs use the bundled known-trigger list because they have no SIF fingerprint. It removes only listed `METACYC` reaction assignments; original annotation tables, feature IDs, sequences, EC assignments and function names are retained. Each attempt records removed assignments and reasons in `ptools-reaction-filter.json` alongside its input diagnostics (native runs save the audit in the entity output directory). Compact SIF runs include this audit in the diagnostic archive.

The bundled known triggers are `TRANS-RXN8J2-121` and `RXN8J2-204`. For `TRANS-RXN8J2-121` in Pathway Tools 29.5, supplying it explicitly imports sublancin precursor proteins before name matching; an imported protein has no input raw-gene record, causing a NIL structure error. A one-ORF reproducer confirmed the failure. Omitting the explicit assignment allowed the build to finish and name matching recovered the same reaction later. This filter prevents the known early-import failure; it does not forbid later inference of the reaction or certify all Pathway Tools inputs. Changes to the blacklist invalidate affected PGDB checkpoints.

## Transport inference and sequence-backed inputs

Transport inference (TIP) is enabled by default for SIF runs. It uses annotations to identify probable transport proteins and associate or create transport reactions. High- and low-confidence predictions are reported separately. Low-confidence or ambiguous annotations are not equivalent to experimentally demonstrated transport.

`--no_transport_inference` disables TIP for a controlled comparison. It is accepted by `ptools` and `analysis_wf`. Record it when comparing pathway/reaction totals; TIP can change reactions without changing the number of base pathways.

MP stages real per-contig DNA sequences, feature coordinates, strands and available Prodigal genetic codes before invoking Pathway Tools. This allows Pathway Tools to derive its protein BLAST database. The BLAST tools **inside the SIF** serve that function; they are separate from FAST/BLAST annotation indexes **in MPDB**.

Annotations with multiple EC assignments are written as separate EC entries. Provisional EC identifiers are retained and may be rejected by Pathway Tools. Recognized tRNA display names separate amino-acid labels from anticodons, while stable feature IDs are preserved. Some ID-parsing warnings can remain. Neither formatting step invents biological assignments.

## Outputs, warnings and failures

Look under `results/SampleA/results/pgdb/community/` and `.../MAGs/MAG_ID/` for the archive, pathway table and pathway-to-ORF table. Each attempt retains `diagnostics/ATTEMPT/execution.json` and available logs. `input/pathologic.log` contains the internal inference messages; `build-xvfb.log` and `export-xvfb.log` contain display-wrapper diagnostics in the corrected version.

| Symptom | What to check or do |
| --- | --- |
| No registered image | Complete `build_pt`, or specify an existing readable `--image` |
| Installer version not recognized | Keep its original filename or supply `--ptools_version` |
| `debconf: delaying package configuration` | This message alone is normal in a minimal image; inspect the final build status |
| Missing compound / inverse-link failure | Check `pathologic.log`, image/patch provenance and scope flags; do not delete reference compounds to suppress it |
| Protein BLAST DB cannot be created | Check sequence staging and the internal log; installing BLAST alone cannot supply missing sequence |
| `Done` followed by failure | Check the recorded phase and wrapper logs; `Done` is not a substitute for successful export/archive |
| Provisional EC or ambiguous transport warning | Keep the evidence; a formatting change cannot establish a missing biological assignment |
| MAG lacks selected genes | The complete workflow can skip a bin with no generated PF input; distinguish this from zero inferred pathways |
| Nonzero exit during export | Preserve the saved PGDB and diagnostics; do not label the result complete from build messages alone |

Community failures fail the workflow. MAG failures may be optional, allowing other entities to finish, but remain failures in task records. `SUCCESS` at the overall scheduler level can coexist with optional failures. Never use it alone to claim all MAGs succeeded. The report distinguishes output availability from latest task outcome.

Failed on-disk PGDBs are retained under `diagnostics/ATTEMPT/failed-pgdbs/`. MP does not automatically certify or publish partial recovery. Retry with the same parameters after fixing the cause, and inspect the new execution record. [Restart guidance](execution.md#logs-temporary-files-and-restarting) explains reuse and targeted reruns.

## Screen reaction compatibility (maintainers)

After building a licensed Pathway Tools image, maintainers can screen the explicit
reaction IDs in an MPDB against that image:

```bash
metapathways screen_pt -d MPDB -o reaction-screen --max_tasks 1
```

The command uses the image registered by `build_pt`; use `--image /path/to/ptools.sif`
to select another image. Each container has a private home and temporary PGDB.
Publication downloads are disabled. `--scratch_dir /path/to/local/scratch` places
these temporary builds on local storage; only inputs, diagnostic logs and receipts
are retained in the output. This is a local maintainer command, not a Slurm workflow.

The screen first builds a no-reaction baseline, then tests batches of 100 reaction
IDs. Failed batches are split until individual triggers are isolated. A candidate
requires two isolated failures and a successful no-reaction control. Timeouts and
uncertain results are recorded as inconclusive. Batch failures whose two halves
pass are recorded separately as possible interactions.

For a targeted check:

```bash
metapathways screen_pt -d MPDB -o reaction-check \
  --reactions TRANS-RXN8J2-121
```

Repeat the same command to reuse completed attempts and retry interrupted or
inconclusive attempts; earlier diagnostics are preserved. The image, mapping table and
screen settings must match the saved checkpoint; use a new output directory when
changing them. `--max_tasks` can be changed on resume. Each build defaults to a
30-minute timeout (`--timeout`, in seconds).

Review `summary.json`, `blacklist-candidates.json` and each attempt's logs before
adding any entry to the shipped blacklist. Standalone candidates are not activated unless you explicitly use `--publish`.
Publication requires a full screen with no inconclusive or interaction failures. Synthetic explicit-ID screening does not certify all enzyme names,
taxonomic contexts or combinations of reactions.

### Automatic compatibility screening during MetaCyc builds

`build_pt -i INSTALLER -d MPDB` and `build_db -d MPDB --func metacyc`
now screen reactions **by default**, after preparing MetaCyc. No additional user
command is needed. The screen tests the licensed SIF and database together, bisects
failed batches, and requires repeated isolated failures plus successful controls.
It adds time to the build and retains screen evidence under
`MPDB/.metapathways/ptools-screens/`. Screening honors `--max_tasks`, capped by the CPU and memory budgets. Each
PTools container uses one CPU; `--memory` is its memory reservation. For example,
`--max_tasks 32 --memory '8 GB'` can screen 32 batches concurrently when 32 CPUs
and 256 GB of memory are available. Local runs default to available resources;
Slurm defaults to four screening containers if `--max_tasks` is omitted.
The outer Nextflow screening task reserves the aggregate CPU and memory for
those containers. On Slurm they run together inside one job allocation, not
as separate cluster jobs. Planning prints the effective concurrency and
reservation. Failed-batch splitting proceeds sequentially within each batch.

Use `--skip_pt_screen` on either build command to opt out. When importing a
MetaCyc data directory instead of a SIF, provide `--screen_image /path/to/ptools.sif`
or first register an image with `build_pt`; directory data alone cannot run the
compatibility checks. An unresolved screen fails the screening stage while
retaining the prepared database and diagnostic checkpoints.

A completed default screen publishes
`MPDB/functional_categories/ptools_reaction_compatibility.json`. MP finds it using
the reference database recorded in the sample run log and checks its mapping and
SIF fingerprints before using it for PGDB builds. It records the selected list in
`ptools-reaction-filter.json`. A mismatch uses only fallback entries confirmed for the selected SIF, rather
than applying an unrelated database-specific list. If neither matches, MP does
not assume other versions share the same defects; rebuild MetaCyc with screening
for the selected image. Compatibility-list changes
invalidate PGDB checkpoints.

To publish an already completed standalone full screen, repeat its command with
`--publish`; completed attempts are reused.

## Intermittent container startup failures

PGDB builds automatically retry an explicit Apptainer container-creation mount
failure up to two times (three attempts total), waiting 5 seconds and then 15
seconds. This applies only before the container shell starts. Input-validation
errors, missing images, and failures during container setup, pathway inference,
export or archiving are not retried by this mechanism. A persistent mount error
still fails the task and needs investigation of the cluster/container runtime.

Each attempt streams to the task log and is saved as `container-attempt-N.log`.
`execution.json` records each attempt's exit code, stage, elapsed seconds and any
retry delay. Compact mode preserves these diagnostics in its existing archives.
Task resource/runtime measurements include all attempts and delays; retain these
records when separating successful execution from infrastructure retry overhead.
No SIF rebuild or additional command-line flag is required.
