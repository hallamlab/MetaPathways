# Prepare reference databases

**`build_pt` builds a licensed image. `build_db` prepares reference databases.**
You only need Pathway Tools setup when you select MetaCyc or want pathway inference.
First [install and activate MP](installation.md). Use a writable MPDB directory
and an internet connection. The tiny test MPDB is for the bundled test only.

## 1. Choose your references

| What you want | What to do first | Functional selection |
| --- | --- | --- |
| MetaCyc only | [Build/register your licensed SIF](pathway-tools.md#build-and-register-once), or supply an existing source below | `--func metacyc` |
| MetaCyc plus public references | Prepare the licensed source first | `--func swissprot metacyc` |
| Public references only | No Pathway Tools setup needed | `--func swissprot`, `cazy`, `uniref50`, or `uniref90` |

A new MPDB includes shared SILVA RNA references, taxonomy and enzyme support,
even when MetaCyc is its only functional database. SwissProt is downloaded only
when selected (or when `--func` is omitted, since SwissProt is the default).

## 2. Build a new MPDB

Choose **one** command below. The destination must be new or empty. For MetaCyc,
wait for the licensed image build to finish first; MP uses the registered image.

MetaCyc only:

```bash
metapathways build_db -d ~/MPDB --func metacyc -a fast
```

SwissProt and MetaCyc together:

```bash
metapathways build_db -d ~/MPDB --func swissprot metacyc -a fast
```

SwissProt without MetaCyc:

```bash
metapathways build_db -d ~/MPDB --func swissprot -a fast
```

Append `--dryrun` to inspect the plan first. It creates logs/plans but does not
download or modify references. Wait for the real build, including MetaCyc reaction
screening, to succeed before analysis. Use `-a blast` instead if annotation will
use `--annotation_algorithm BLAST`.

## 3. Select the references when analyzing

Building a reference does not automatically select it. Use
`--annotation_dbs metacyc` for MetaCyc-only annotation, or
`--annotation_dbs swissprot metacyc` for both, on `run` or `analysis_wf`.
Continue with [input organization](inputs.md) and [analysis](analysis.md).
Adding a database later does not change already-generated annotations.

## Add to an existing MPDB

Use `--add` for references that are not already present. Existing references and
shared support files are retained. If newly licensed, build/register your SIF
first, then add MetaCyc:

```bash
metapathways build_db -d ~/MPDB --add --func metacyc -a fast
```

To add a public database, no Pathway Tools installation is required:

```bash
metapathways build_db -d ~/MPDB --add --func cazy -a fast
```

A normal build pointed at an existing MPDB stops and suggests an explicit
operation. `--add` refuses to replace an existing reference.

## Update or resume

**Stop analyses reading this MPDB before updating it.** Use `--update` to rebuild
selected existing references, including a new MetaCyc release or a change of
index format. Other functional references and shared support are retained.
For example, after registering a new licensed image:

```bash
metapathways build_db -d ~/MPDB --update --func metacyc -a fast
```

To recover an interrupted create/add/update, repeat its **original command**.
Keep the same destination, selected references, licensed source and index format.
Successful compatible tasks are reused. An unfinished operation must be completed
before starting a different one. MP permits only one database builder per MPDB.
An incomplete older MPDB without a restart record should be rebuilt in a new
directory; arbitrary nonempty directories are rejected.

## Already have a licensed MetaCyc source?

You do not need to rebuild an image. Supply an existing SIF explicitly:

```bash
metapathways build_db -d ~/MPDB --func metacyc -a fast \
  --metacyc_source /path/to/pathway-tools.sif
```

Use `--add` or `--update` if appropriate for the destination. Alternatively,
`--metacyc_source /path/to/data` accepts a complete local export containing
`protseq.fsa`, `proteins.dat`, `enzrxns.dat`, `reactions.dat`, `pathways.dat`,
`compounds.dat` and `classes.dat`. A protein FASTA alone is insufficient.
Default reaction screening also needs a SIF: supply `--screen_image` or register
one. `--skip_pt_screen` explicitly omits that screening; see its
[scope and limitations](pathway-tools.md#defaults-and-explicit-choices).
Missing sources or required files stop planning before reference downloads begin.

## Other reference choices

UniRef databases are substantially larger than the test references. eggNOG needs
a supplied FASTA at `MPDB/functional/eggnog`; MP does not download it. Place that
file in an otherwise new MPDB before selecting `--func eggnog`. Public downloads
are not all version-pinned. Preserve release records and checksums with your
analysis, and keep licensed sources according to their vendor terms.
