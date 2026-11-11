# Pathway Tools maintenance and diagnostics

[Development](development.md) · [Maintainer CLI reference](developer-cli-reference.md)

For contributors validating images, investigating failures or maintaining reaction compatibility lists. Normal users get automatic screening through `build_db`; see the [user guide](pathway-tools.md).

## Screen reaction compatibility

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

`build_db -d MPDB --func metacyc`
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

Use `--skip_pt_screen` on `build_db` to opt out. When importing a
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

## Taxonomic comments in PathoLogic input

MP writes per-feature taxonomic evidence as a `GENE-COMMENT`, preserving existing
comments. It no longer emits the `TAXONOMIC-ANNOT` input slot rejected by the
alpha tester's PathoLogic loader. Older inputs are converted in private staging;
original files and MP's taxonomic annotation tables remain intact. These comments
are evidence text, not a substitute for the organism taxonomic scope above.
The comment field is documented in [SRI's PathoLogic input example](https://brg.ai.sri.com/ptools-image-07/Curating.pdf).

## PGDB input preparation

MP stages real per-contig DNA sequences, feature coordinates, strands and available Prodigal genetic codes before invoking Pathway Tools. This allows Pathway Tools to derive its protein BLAST database. The BLAST tools **inside the SIF** serve that function; they are separate from FAST/BLAST annotation indexes **in MPDB**.

Annotations with multiple EC assignments are written as separate EC entries. Provisional EC identifiers are retained and may be rejected by Pathway Tools. Recognized tRNA display names separate amino-acid labels from anticodons, while stable feature IDs are preserved. Some ID-parsing warnings can remain. Neither formatting step invents biological assignments.
