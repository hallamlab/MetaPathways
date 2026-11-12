# Vitamin B5: from setup to results

This worked example starts with a function rather than an organism name. You will
search for a candidate enzyme, find its contig and genome bin, inspect ribosomal
identity evidence, check inferred pathways, and export tables for follow-up work.
You can explore completed results immediately or run the analysis yourself.

```{contents} In this vignette
:local:
:depth: 2
```

## 1. Choose your route

**Explore completed results:** install MP on Linux x86-64, then
[download the example](#3-download-the-completed-example). You need a browser,
but no Pathway Tools installation or MPDB.

**Run everything yourself:** install MP, then [build your own results](#4-build-your-own-results).
You also need your own licensed Pathway Tools installer, a compatible Apptainer
host, internet access, and storage for references.

Both routes meet at [Open the report](#5-open-the-report). **Choose one route.**
You do not need to repeat the analysis to browse the completed example.

The question is: **which fragment carries a plausible part of the machinery for
making pantothenate (vitamin B5)?** We begin with pantothenate synthetase, PanC,
also identified by EC 6.3.2.1. Its annotation is a lead to investigate, not a
measurement of vitamin production.

The inputs are three small, simulated CAMI II samples: `Gastrointestinal_5`,
`Skin_28`, and `Urogenital_22`. Each contains three 50 kb assembly fragments,
paired reads, and supplied contig-to-genome assignments. These bins are cropped
source-genome fragments, **not complete recovered MAGs**. Missing genes may lie
outside the crop. The [input provenance](test-bundle.md) explains selection and
limitations; this example cannot establish pathway completeness, nutritional
independence, or differences between human body sites.

The downloadable snapshot was analyzed with **MP 4.0.8, SwissProt 2026_04,
MetaCyc/Pathway Tools 29.5, and SILVA 138.1**. Its explorer was rebuilt with
**MP 4.0.13**. Screenshots below show that snapshot. Reference updates can change
a fresh run's results; exact row counts below are checks for the supplied snapshot.

## 2. Install and activate MP

Use a Linux x86-64 workstation or server. If Mamba is not installed, first follow
the [installation guide](installation.md). On an HPC system, use the site's
supported environment activation method. Run these commands in the terminal on
the machine that will hold the results:

```bash
mamba create -n mp-vignette --override-channels --strict-channel-priority \
  -c hallamlab -c conda-forge -c bioconda metapathways=4.0.13
conda activate mp-vignette
metapathways version
```

Confirm that MP reports version 4.0.13. Activate `mp-vignette` again in each new
terminal. If you already have this version in an environment, activate that
instead of creating another one.

```{figure} assets/vignette/terminal-install.png
:alt: Terminal command reference for installing MP 4.0.13 and activating the environment.

Installation command reference. Terminal command images accompany copyable code;
they do not depict a fresh installation performed for this vignette.
```

## 3. Download the completed example

This route needs no licensed installer, container runtime check, or reference
build. The compressed download is approximately **4.1 MiB**. Download [the example archive](https://github.com/hallamlab/MetaPathways/releases/download/v4.0.13/mp-vignette-v1.tar.gz)
and [its checksum](https://github.com/hallamlab/MetaPathways/releases/download/v4.0.13/mp-vignette-v1.tar.gz.sha256),
or use the commands below. Choose a new folder to avoid mixing runs:

```bash
mkdir -p ~/mp-vignette
cd ~/mp-vignette
curl -fL -O https://github.com/hallamlab/MetaPathways/releases/download/v4.0.13/mp-vignette-v1.tar.gz
curl -fL -O https://github.com/hallamlab/MetaPathways/releases/download/v4.0.13/mp-vignette-v1.tar.gz.sha256
sha256sum -c mp-vignette-v1.tar.gz.sha256
tar -xzf mp-vignette-v1.tar.gz
cd mp-vignette-v1
```

Continue only if the checksum reports `OK`. The archive contains `cami-test/`
(inputs and their provenance), `all/` (three samples and the report),
`provenance/` (reference release records), and `manifest.json` (file checksums).
It also contains worked CSV exports for comparison.

```{figure} assets/vignette/terminal-example.png
:alt: Commands to download, check, extract and open the completed example.

Download command reference. Once extracted, skip the licensed setup below and
continue at section 5.
```

The bundle contains analysis-derived annotation and pathway tables, sequences,
GFF/GenBank files, and RNA/abundance evidence. It excludes Pathway Tools software,
the SIF, MPDB references, PGDB databases and archives, alignment files, and
intermediate diagnostic workspaces. “Complete dataset” in Downloads means the
files retained in this bundle; it cannot restore excluded products. Public MP
packages do not provide a licensed Pathway Tools installation.

## 4. Build your own results

Skip this section if you downloaded the completed example. These steps produce
your own results from the included inputs with full SwissProt and MetaCyc
references. The inputs are tiny, but the reference downloads and initial
MetaCyc compatibility screen are substantially larger tasks. Allow enough disk
space for the MP environment, reference databases, image build/cache and results.
We do not provide a clean-install timing estimate: the recorded analysis reused
some completed tasks and would give a misleading estimate.

### 4.1 Check the machine before obtaining the installer

```bash
mkdir -p ~/mp-vignette-run
cd ~/mp-vignette-run
metapathways build_pt --check -o containers
```

Wait for the readiness check to pass. MP tries the environment's Apptainer and
available site installations automatically. A Mamba installation cannot grant
permissions denied by the host. If the check fails, use the printed diagnostic
and [host prerequisites](pathway-tools.md#host-prerequisites) before continuing.
You can still use the completed-example route while resolving host setup.

```{figure} assets/vignette/terminal-check.png
:alt: Commands to create a separate workspace and check container readiness.

Readiness command reference. A successful check precedes licensed setup.
```

### 4.2 Obtain your license and build the image

Request Pathway Tools through the official [licensing page](https://bioinformatics.ai.sri.com/ptools/licensing/ptools-academic-license.shtml),
using the route appropriate to your organization. Follow SRI's download
instructions and obtain a **Linux x86-64 distribution containing MetaCyc**.
The recorded example used 29.5. Do not substitute a private download link into
public scripts or documentation.

Place your installer in `~/Downloads/` on the analysis machine. If you downloaded
it to a laptop, [copy it to the server](pathway-tools.md#get-the-installer) first.
Replace the filename below if yours differs. Do not run or unpack it manually:

```bash
metapathways build_pt \
  -i ~/Downloads/pathway-tools-29.5-linux-64-tier1-install \
  -o containers
```

MP builds the licensed SIF and registers its location after a successful build.
Keep that image at its registered location. The next command discovers it
automatically. A readiness check alone does not register an image.

```{figure} assets/vignette/terminal-pt.png
:alt: Command to build a licensed Pathway Tools SIF from the user's installer.

Image-build command reference. Each user supplies their own licensed installer.
```

Already have a working SIF? Follow [existing-image setup](pathway-tools.md#host-prerequisites)
instead of rebuilding. HPC sites may allow image execution but require image
construction elsewhere. That alternative does not change the biological steps.

### 4.3 Build the full reference database

After the image build succeeds, run this from `~/mp-vignette-run`:

```bash
metapathways build_db -d MPDB --func swissprot metacyc -a fast \
  --max_cpus 8 --max_memory '32 GB'
```

This creates a new MPDB containing both functional references plus shared SILVA,
enzyme, and taxonomy support. It also screens explicit MetaCyc reaction
assignments for compatibility. Screening can take time without frequent terminal
updates; consult the printed task logs. Do not start a second builder in the same
MPDB. To resume an interrupted build, repeat its original command.

Use resource limits available on your machine. These limits are scheduling
budgets, not a guarantee that every reference task fits on every host. If you
already have a completed MPDB, use [database add/update guidance](databases.md)
rather than issuing a new build into it.

```{figure} assets/vignette/terminal-db.png
:alt: Command to build SwissProt and MetaCyc references with CPU and memory limits.

Reference-build command reference. Wait for all required tasks to succeed.
```

### 4.4 Prepare and inspect the inputs

```bash
metapathways prepare_test -o inputs
cat inputs/cami-test/all.tsv
```

The manifest has one row per sample, naming an assembly, two read mates, and a
MAG map. Paths are relative to the manifest. Keep the `cami-test` directory
intact. The preparation command also copies a tiny test MPDB under `inputs/`;
**we do not use that database here**. The analysis below selects the full `MPDB`
you just built.

```{figure} assets/vignette/terminal-inputs.png
:alt: Commands to prepare the included inputs and inspect their manifest.

Input-preparation command reference. Expect three distinct sample rows.
```

### 4.5 Run annotation, abundance, MAG splitting and pathways

```bash
metapathways analysis_wf \
  --manifest inputs/cami-test/all.tsv -o all -d MPDB \
  --annotation_dbs swissprot metacyc \
  --threads 4 --memory '4 GB' --max_tasks 2 \
  --max_cpus 8 --max_memory '32 GB'
```

This command includes pathway inference using the registered SIF. Default RNA
reference selection uses the full SILVA references prepared with MPDB. Do not
add the tiny test database names or `--skip_ptools`. Four GB per task was used
for this small example; it is not a recommendation for full metagenomes.

```{figure} assets/vignette/terminal-analysis.png
:alt: Complete workflow command using the three-sample manifest and both functional references.

Analysis command reference. Keep this terminal session alive until completion.
```

Review the final task status and the printed `logs/analysis_wf/RUN_ID/summary.json`.
Required failures must be resolved before interpreting results; a generated
report alone does not prove that every MAG pathway task succeeded. Repeating the
same command reuses compatible completed work.

```{figure} assets/vignette/terminal-verified.png
:alt: Verified historical summary: 40 successful tasks, 29 reused tasks, and the expected vitamin B5 evidence.

Verified summary of the supplied MP 4.0.8 run, transcribed from its saved summary
and result tables. This is historical evidence, not output from a new 4.0.13 run.
```

## 5. Open the report

From the directory containing `all/`, run:

```bash
metapathways report -o all --serve
```

This rebuilds the index from retained results and opens the explorer. It does not
repeat annotation or pathway inference. Keep the terminal running. On a remote
server instead use:

```bash
metapathways report -o all --serve --no-browser --port 8765
```

In a second terminal **on your laptop**, replace `user@server` with your SSH host:

```bash
ssh -N -L 8765:127.0.0.1:8765 user@server
```

Leave both terminals open and use the complete local URL printed by MP, including
its access token. See [remote report access](reports-tutorial.md#view-a-remote-report-through-ssh)
if your site uses a jump host or the port is busy. Stop the server with Ctrl-C
when finished.

```{figure} assets/vignette/terminal-report.png
:alt: Report launch command and optional SSH forwarding commands.

Report command reference. Use the local route or the remote route, not both.
```

Open **Run report** to check accounting, then return to **Explore and export**.
The Samples table should include all three sample IDs. The footer records the
historical analysis version and status; rebuilding the explorer does not rewrite
that version.

```{figure} assets/vignette/01-samples.png
:alt: Actual example explorer showing the three CAMI samples.

The completed example opens at Samples. Check the dataset before searching.
```

## 6. Start with the function

In **Explore a table**, choose **Functional annotations** from the annotation
group. Start with **All samples**, no MAG selection, and no column filters.
Enter `pantothenate synthetase` in **Search all fields** and press Enter (or click
**Apply filters**).

Open **Choose table and export columns** and retain `sample_id`, `orf_id`,
`contig_id`, `reference_db`, `product`, `ec`, and `reaction`. This also controls
the CSV columns. The supplied example returns **two rows for one gene**:
`Gastrointestinal_5 / C2-G10`. SwissProt supplies EC `6.3.2.1`; MetaCyc supplies
`PANTOATE-BETA-ALANINE-LIG-RXN`. Two reference annotations do not mean two genes.

```{figure} assets/vignette/02-function.png
:alt: Two functional annotation rows for C2-G10, one from MetaCyc and one from SwissProt.

A function search supplies the sample and contig address for the investigation.
```

Export these rows using **Export all matching rows as CSV**. Rename the downloaded
file `panc-annotations.csv`. Save the query definition alongside it so you can
reconstruct the filters. An exact `ec equals 6.3.2.1` column filter is an
alternative search, but it returns only the SwissProt row in this snapshot.

## 7. Find the contig's MAG and identity evidence

### 7.1 Discover membership

On either PanC row, click **Contig MAGs**. The related view supplies exact sample
and contig filters. It returns one membership:

| Sample | Contig | MAG |
| --- | --- | --- |
| Gastrointestinal_5 | Gastrointestinal_5-C2 | CAMI_OTU_97_68_0 |

```{figure} assets/vignette/03-membership.png
:alt: Contig-to-MAG membership linking Gastrointestinal_5-C2 to CAMI_OTU_97_68_0.

The MAG identifier is a bin label. It does not itself identify the organism.
```

Keep the sample with the MAG identifier: identifiers can recur between samples.
The **Sample** and **MAG** selectors can restrict subsequent biological tables;
this example also uses explicit column filters so the exact query is visible.

### 7.2 Inspect both ribosomal markers

From the membership row, click **Contig rRNA taxonomy**. Retain `sample_id`,
`contig_id`, `subunit`, `similarity`, and `taxonomy` as visible columns. Both
16S and 23S evidence rows refer to this same contig and report 100% similarity to
*Roseburia hominis* A2-183 in their SILVA matches.

```{figure} assets/vignette/04-rna.png
:alt: Linked 16S and 23S evidence on the PanC contig, both matching Roseburia hominis A2-183.

The ribosomal evidence places an identity clue on the fragment carrying PanC.
```

This is stronger context than treating the species name of a protein reference
hit as the identity of the query organism. The PanC SwissProt reference is from
*Ruminiclostridium cellulolyticum*; that does not rename this contig. RNA hits
remain locus-level evidence, not an automatic taxonomic assignment to an entire MAG.

For other RNA features, **Contig RNA** includes tRNAs with their reported isotype
and anticodon. tRNAs do not receive a SILVA taxonomic assignment.

## 8. Inspect pathways and the gene neighborhood

### 8.1 Check the inferred pathway

Return to **Contig-to-MAG membership**, restore the sample and contig filters if
needed, and click **Entity pathways** on the MAG row. Search for `PANTO-PWY`.
Show `sample_id`, `entity_id`, `pathway_id`, `name`, `reactions`,
`covered_reactions`, and `linked_orf_count`.

The MAG has **phosphopantothenate biosynthesis I**, with **2 of 4 reactions
covered** and two linked genes. This is a MAG-specific inference; community
pathway rows would describe a different entity.

```{figure} assets/vignette/05-pathway.png
:alt: PANTO-PWY for the selected MAG with two covered reactions out of four.

The pathway inference shows supported pieces and gaps. Its score is not a
percentage completeness or a probability of biological truth.
```

### 8.2 Follow supporting genes and export the second table

Click **Pathway genes** on that row. Retain `sample_id`, `entity_id`, `pathway_id`,
`orf_id`, `contig_id`, and `product`. The example links **C2-G9 (PanB)** and
**C2-G10 (PanC)**. Export all matching rows and rename the file
`panto-pathway-genes.csv`; save its query definition too.

```{figure} assets/vignette/06-genes.png
:alt: Two pathway-gene associations linking C2-G9 and C2-G10 to PANTO-PWY.

This export has one row per pathway–gene association, not per database annotation.
```

### 8.3 Inspect the neighbors

Click **Contig ORFs**. Keep the exact sample and contig filters. Add `start ≥
10576` and `start ≤ 12296`, then apply. Show `sample_id`, `orf_id`, `start`, `end`,
and `product`. The window contains G9, G10 and G11; G11 is an aspartate
decarboxylase (PanD) candidate. Its **Pathway links** include `PWY-5155`,
β-alanine biosynthesis III, with one covered reaction out of one in this MAG.

```{figure} assets/vignette/07-neighborhood.png
:alt: Neighboring G9, G10 and G11 annotations on the selected contig.

Neighboring precursor-related functions make a coherent lead for follow-up;
proximity alone does not establish an operon or expression.
```

Click **FASTA** on a feature, then choose an available nucleotide or amino-acid
record in the pop-out. **Download full FASTA** saves that individual sequence.
The viewer loads one selected record at a time; bulk sequences belong in Downloads.

```{figure} assets/vignette/08-fasta.png
:alt: A retained individual sequence opened in the FASTA viewer.

Inspect a single sequence and its source without loading the whole metagenome.
```

## 9. Use the exports downstream

You now have two small CSV files. A useful next step is an evidence table listing
the pathway's genes and which have the PanC annotations you exported. Save the
following as `summarize_evidence.py` beside both CSVs (or
{download}`download the script <assets/vignette/summarize_evidence.py>`). It uses Python's standard
library, so no additional analysis package is needed:

```python
import csv
from collections import defaultdict

annotations = defaultdict(set)
with open("panc-annotations.csv", newline="") as stream:
    for row in csv.DictReader(stream):
        annotations[(row["sample_id"], row["orf_id"])].add(row["reference_db"])

with open("panto-pathway-genes.csv", newline="") as stream:
    genes = list(csv.DictReader(stream))

with open("panto-evidence.csv", "w", newline="") as stream:
    fields = ["sample_id", "entity_id", "pathway_id", "orf_id", "product",
              "databases_in_panc_export"]
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    for gene in genes:
        evidence = annotations[(gene["sample_id"], gene["orf_id"])]
        writer.writerow({**{key: gene[key] for key in fields[:-1]},
                         "databases_in_panc_export": ";".join(sorted(evidence))})
```

Run it from that folder:

```bash
python summarize_evidence.py
```

```{figure} assets/vignette/terminal-export.png
:alt: Actual downstream summary output showing two pathway genes and the PanC database evidence.

Captured output from running the script on the browser-exported CSVs.
```

The resulting `panto-evidence.csv` has two rows. G10 has `metacyc;swissprot`;
G9's last field is blank because the first export was restricted to PanC—not
because G9 lacks annotations. This preserves the pathway's gene count while
summarizing multiple annotation records. Always join with `sample_id` as well as
`orf_id`, and preserve `entity_id` for MAG-specific interpretations.

For an optional abundance follow-up, use **ORF abundance** from G10. In **Read
abundance**, keep the sample, feature ID and feature type, and select the desired
measurements. Export that table to annotate the evidence summary with coverage
or TPM. Do not sum values after expanding genes into multiple annotation or
pathway rows. The cropped reads are for demonstrating measurement navigation,
not testing habitat differences.

## 10. Download the fragment for further work

Open **Downloads** in the page navigation. Under `Gastrointestinal_5`, select
only `CAMI_OTU_97_68_0`. Choose **Complete dataset**, click **Review selection**,
and inspect the inventory before **Download selected results (ZIP)**.

```{figure} assets/vignette/09-download.png
:alt: Downloads page with the discovered MAG selected and its generated file inventory.

The MAG archive includes its mapped sequences, annotations, and pathway tables,
with available GFF/GenBank subsets. The review describes the packaged files.
```

Use **EDA tables only** when you want CSV tables and sequence records without
additional retained source files. You can also select whole samples or combine
MAGs across samples. These choices are independent of explorer text filters:
searching PanC does not restrict a MAG archive to PanC. Read `manifest.json`
inside the download for provenance and unavailable subsets.

The teaching bundle deliberately excludes original PGDB archives and BAMs.
A full run you produce yourself can retain additional products. An archive from
this example does not promise those absent files.

## 11. Interpret what you found

A defensible conclusion is: **this Roseburia-associated fragment carries a
plausible part of a vitamin B5-related biosynthetic toolkit, supported by
neighboring gene annotations and partial pathway inference.**

The next biological question is whether the remaining steps occur elsewhere in
the complete genome and whether suitable transporters are present. A complete
assembly or high-quality MAG would be the next computational input; nutritional
experiments would be needed to test actual dependence or production. This crop
does not establish vitamin secretion, cross-feeding, or host benefit.

You have followed function → contig → bin → ribosomal evidence → pathway →
supporting genes, then exported a reproducible subset for downstream work.
Use the same pattern with your own target function and samples, retaining the
identifiers and evidence boundaries throughout.

## References and provenance

- CAMI input selection and attribution: [test bundle](test-bundle.md) and its included `provenance.json`.
- Biochemical entry point: [ENZYME EC 6.3.2.1](https://enzyme.expasy.org/EC/6.3.2.1).
- Biological context for follow-up: Hillman et al. (2020), [comparative genomics of Roseburia](https://doi.org/10.1099/mgen.0.000399). Literature motivates the question; the local fragment does not establish a full phenotype.
- [Pathway Tools licensing](https://bioinformatics.ai.sri.com/ptools/licensing/ptools-academic-license.shtml). Users obtain their own software; the example distributes analysis-derived tables, not the licensed software or reference database.
- Snapshot release records and original successful execution summary are inside the example. Historical absolute paths describe the source run and are not required on your machine.
