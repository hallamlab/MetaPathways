# Pathway Tools: licensing, image builds, databases and inference

[Home](index.md) · [Test walkthrough](test.md) · [Complete flags](cli-reference.md)

## Choose your setup before running commands

`build_pt` builds and registers the licensed SIF. `build_db` prepares MPDB
references. Use them separately, in that order when starting with a licensed
installer.

| What you want | Setup |
| --- | --- |
| MetaCyc protein annotation | Build/register the SIF below, then [build an MPDB with MetaCyc](databases.md) |
| Pathway inference from existing annotations | Build/register the SIF; adding MetaCyc protein references is optional |
| Add MetaCyc to an existing MPDB | Prepare a licensed source, then use `build_db --add --func metacyc` |
| Public annotation without pathway inference | Skip this page; use public references and `--skip_ptools` |

If you already have a working SIF or complete MetaCyc export, see
[existing sources](databases.md#already-have-a-licensed-metacyc-source).
Do not rerun the image build merely to add a database.

````{advanced} inference settings and defaults

(defaults-and-explicit-choices)=

```{rubric} Defaults and explicit choices
```

| Setting | MP default | Your alternative |
| --- | --- | --- |
| Reaction compatibility screening when building MetaCyc | Enabled; successful full screens automatically publish an MPDB compatibility list | `--skip_pt_screen` on `build_db` |
| Taxonomic pruning for PGDB inference | Enabled; avoids the unpruned rescoring pass | `--no_taxprune` |
| Organism taxonomic scope | `all`: NCBI 131567, cellular life; broad enough for mixed-domain inputs | `--taxonomic_scope bacteria`, `archaea`, `eukaryotes`, or `--taxon_id ID` |
| Transport inference in SIF runs | Enabled | `--no_transport_inference` |
| CPUs per PGDB | One; concurrency uses independent PGDBs | Control concurrency through the workflow resource flags |

These defaults apply to `ptools` and `analysis_wf`. Taxonomic scope guides PGDB
inference without changing the gene-level taxonomy in MP annotation tables.
The compatibility list removes only unsafe **explicit reaction assignments**,
not genes, annotations or sequences; the same reaction may still be inferred
later. The sections below explain the evidence, limits and diagnostic records.

````

````{advanced} building inside the MP container

(building-from-the-mp-apptainer-container)=

```{rubric} Building from the MP Apptainer container
```

The release container includes an unlicensed Ubuntu dependency base. `build_pt`
selects it automatically, avoiding package installation inside the restricted
nested build. Your licensed installer and official patches are still required
and remain private. See [nested containers](containers.md#nested-pathway-tools-containers)
for tested behavior and host limitations. On hosts affected by Java/FUSE stalls,
launch the outer MP SIF with `--unsquash`; this needs temporary extraction space.

````

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

With the Mamba environment active, run this once on the machine where you will
build the image. No installer, network download or sudo is needed for this check:

```bash
metapathways build_pt --check -o ~/mp-containers
```

MP first tries Apptainer from the active Mamba environment, then other installations
on `PATH` and standard system paths. It tests a small local container build and
execution before selecting a runtime. A working choice is saved per host and MP
environment and checked again in each new process. No PATH edits are needed.
Subsequent checks try that saved choice first and fall back if it no longer works.

The normal `build_pt` command also checks automatically before downloading licensed
patches. Slurm workers select and test their own runtime on the allocated node;
a login-node check cannot certify other nodes. Load your site's Apptainer module
before submitting if it is not available in the standard paths. MP does not install
global software or change host permissions.

Diagnostics are under the requested output's `logs/readiness/`: `runtime-selection.json`
records every attempted executable, and each candidate has its own log. If all
candidates fail, send that report and the candidate logs to your administrator.
For an explicit site requirement, set `METAPATHWAYS_APPTAINER` to its executable;
MP then tests only that choice, without silently switching to another runtime.

Older Apptainer installations without `--mksquashfs-args` use their site-configured
compression CPU settings instead of the `-t` compression limit.

**Already have a licensed SIF, or your HPC permits running but not building?**
Test execution alone, without an installer or build permissions:

```bash
metapathways build_pt --check --image /path/to/pathway-tools.sif -o ~/mp-containers
```

A successful check selects a runtime; it does not register the image. Supply
`--image /path/to/pathway-tools.sif` to `ptools` or `analysis_wf`, or
`--metacyc_source /path/to/pathway-tools.sif` to `build_db`.
Analyses also select a runtime automatically and test execution of the supplied SIF,
without requiring that the host can build images.

For the actual build, MP automatically downloads and caches a verified public
Ubuntu dependency image. It contains no Pathway Tools or MetaCyc. You do not need
to run APT, Docker, or install OS packages yourself. Network access is needed for
that first download and the official patches; allow disk space for the installer,
cached base, temporary build and finished SIF. The cache is under
`~/.cache/metapathways/ptools-base` (or `$XDG_CACHE_HOME`). An existing copy can be
supplied with `--dependency_base /path/to/base.sif`; official patch downloads
still require network access.

**If the check reports host permissions:** Mamba installed the runtime, but the
host must allow it to create containers. Ubuntu 24.04 and later can restrict
unprivileged user namespaces through AppArmor. An administrator must resolve
that restriction; MP does not change host security settings. Send them the printed
diagnostic log and the [Apptainer host requirements](https://apptainer.org/docs/admin/latest/user_namespace.html).
Fakeroot UID mappings/helpers may also need site configuration, depending on the
mode available. No-sudo users can build on hosts where this is already permitted.

If your site permits running containers but not building them, build the licensed
SIF on a compatible machine and transfer it, then select it with `--image`.
A transferred SIF still needs a host that permits Apptainer execution.
Public annotation and the bundled installation test do not require Pathway Tools.

Use shared storage for a finished image that compute nodes need to read. Local scratch can accelerate building on slow network filesystems. `APPTAINER_TMPDIR` is an optional Apptainer troubleshooting setting, not required to run MP. See your site's Apptainer setup when choosing scratch; it must have enough space.

## Build and register once

**Skip this step if your working SIF is already registered.** Obtain the installer
and check the host prerequisites above before running this command.

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

## Prepare MetaCyc annotation references

After the image build succeeds, follow [database setup](databases.md). That guide
covers MetaCyc-only, mixed references, and additions to an existing MPDB.
`build_pt` does not create or modify MPDB. The retired `build_pt -d` option stops
with a migration message; use `build_db` for all reference selections.

Building MetaCyc is optional for inference from existing SwissProt annotations.
To search MetaCyc proteins, select `--annotation_dbs metacyc` (or
`--annotation_dbs swissprot metacyc` for both) when annotating your samples.

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

SIF tasks receive private home, data and temporary state. MP also isolates X-display sockets; independent PGDBs can run together without the native installation's shared state. Both `ptools` and `analysis_wf` use the registered SIF or an explicit `--image`. A missing SIF stops planning with setup instructions; MP does not silently fall back to a native installation.

````{advanced} taxonomic scope

(choose-a-taxonomic-scope)=

```{rubric} Choose a taxonomic scope
```

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

````

## Explicit reaction blacklist

Before pathway inference, MP automatically filters explicit reaction assignments
known to cause problems with the selected Pathway Tools image. It applies the
filter to private copies of community and MAG inputs, preserving your original
annotations, sequences, EC assignments and function names. Pathway Tools can still
infer the same reactions later through other evidence.

Removed assignments and their reasons are recorded in `ptools-reaction-filter.json`
with the run diagnostics, including the diagnostic archive for compact runs.

````{advanced} transport inference

(transport-inference-and-sequence-backed-inputs)=

```{rubric} Transport inference and sequence-backed inputs
```

Transport inference (TIP) is enabled by default for SIF runs. It uses annotations to identify probable transport proteins and associate or create transport reactions. High- and low-confidence predictions are reported separately. Low-confidence or ambiguous annotations are not equivalent to experimentally demonstrated transport.

`--no_transport_inference` disables TIP for a controlled comparison. It is accepted by `ptools` and `analysis_wf`. Record it when comparing pathway/reaction totals; TIP can change reactions without changing the number of base pathways.



````

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
| Zero inferred pathways | A valid empty result produces header-only pathway and pathway-to-ORF tables. MP verifies the empty final pathway report and inference evidence list before handling the Pathway Tools empty-export message; raw PGDB files are preserved |
| MAG lacks selected genes | The complete workflow can skip a bin with no generated PF input; distinguish this from zero inferred pathways |
| Nonzero exit during export | Preserve the saved PGDB and diagnostics; do not label the result complete from build messages alone |

Community failures fail the workflow. MAG failures may be optional, allowing other entities to finish, but remain failures in task records. `SUCCESS` at the overall scheduler level can coexist with optional failures. Never use it alone to claim all MAGs succeeded. The report distinguishes output availability from latest task outcome.

Failed on-disk PGDBs are retained under `diagnostics/ATTEMPT/failed-pgdbs/`. MP does not automatically certify or publish partial recovery. Retry with the same parameters after fixing the cause, and inspect the new execution record. [Restart guidance](execution.md#logs-temporary-files-and-restarting) explains reuse and targeted reruns.

## Automatic compatibility screening during MetaCyc builds

`build_db --func metacyc` automatically checks reaction compatibility with the
selected licensed SIF. This adds time to database setup and helps prevent known
inference failures. No separate screening command is needed. If screening fails,
keep the diagnostic logs and repeat the same database command after resolving
the cause; prepared references are retained.

````{advanced} screening resources and alternative sources

Screening honors `--max_tasks`, `--max_cpus` and `--max_memory`; each container
uses one CPU and the requested `--memory`. On Slurm, screening containers run
inside one job allocation. MP prints the effective reservation before execution.
When importing a MetaCyc data directory, provide `--screen_image` or register a
SIF first. `--skip_pt_screen` on `build_db` explicitly skips compatibility testing.

````

````{advanced} container retry diagnostics

(intermittent-container-startup-failures)=

```{rubric} Intermittent container startup failures
```

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

````
