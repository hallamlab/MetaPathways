# Reviewer walkthrough: three tiny CAMI samples

[Home](index.md) · Previous: [getting started](getting-started.md) · Next: [your inputs](inputs.md)

This walkthrough uses ordinary MP commands and three distinct CAMI II samples included in the repository and installed package. No Pathway Tools license, installer or container is needed for the default commands. They run annotation, paired-read mapping and abundance, genome splitting, and the report/explorer. Pathway inference is explicitly skipped; licensed users can enable it in the optional section below.

The input bundle is about **2.4 MiB**, separate from software and reference-database downloads. Each sample contains three 50 kb assembly regions assigned to three CAMI source genomes:

| Sample | Assembly bases | Read pairs | Genome bins | Manifest |
| --- | ---: | ---: | ---: | --- |
| Urogenital_22 | 150,000 | 3,000 | 3 | `single.tsv` |
| Gastrointestinal_5 | 150,000 | 2,326 | 3 | `pair.tsv` |
| Skin_28 | 150,000 | 3,000 | 3 | `pair.tsv` |

These are real subsets of simulated CAMI reads and gold-standard assemblies, not duplicated samples. Bins contain small genome fragments, not complete genomes or recovered MAGs. We selected abundant contigs and matching reads to exercise the interfaces. Do not use these coverage-biased subsets to assess biological accuracy, abundance, genome completeness, pathway recovery or publication performance.

Protein taxonomy supports SwissProt (including the small reviewer database), UniRef and eggNOG. Each annotation row reports the taxon of its own reference hit in `taxonomy`; `lca_taxonomy` is computed from score-qualified hits in that same database, with support counted independently for each database. The primary ORF row follows its selected annotation's `reference_db` and target. No database priority or cross-database fallback is used. These describe reference-hit evidence, not a definitive organism assignment for the query. Unsupported databases report `Not computed`; missing or unknown taxon IDs in supported databases report `Unclassified`. An actual LCA of `root` remains `root`. SILVA rRNA taxonomy is separate. Existing outputs need annotation-table regeneration to gain these fields.

## 1. Install MP and prepare the small reference database

Install with [Conda/Mamba or GitHub source](installation.md), then activate the MP environment. Container users can follow the same test through [Docker or Apptainer](containers.md).

```bash
metapathways prepare_test -o ~/mp-reviewer
cd ~/mp-reviewer
metapathways build_db --test -d MPDB
```

`prepare_test` copies the included assemblies, paired reads, genome maps, manifests and reference FASTAs into your workspace. It leaves matching files in place and refuses to overwrite changed files. Use a new workspace for a fresh test. `build_db` indexes the small SwissProt/SILVA references and downloads enzyme/taxonomy support records into `MPDB`. Preparation requires internet access and writable workspace storage; it does not modify your installation.

The bundled `swissprot_test` contains 227 genuine SwissProt proteins: the 39 original K12 test references plus 188 reference targets matched by the three CAMI samples in a full-SwissProt run. The protein FASTA is about 109 KB. This gives the reviewer useful annotation, taxonomy and genome-splitting examples without a full protein database download. Selection and reference provenance are recorded in [`reviewer_reference.json`](https://github.com/hallamlab/MetaPathways/blob/HEAD/metapathways/regtests/test_db/reviewer_reference.json). The deliberately selected reference set is for workflow testing, not annotation accuracy or biological benchmarking. SILVA fixtures remain small; empty RNA tables can be legitimate.

When updating an existing installation, rerun `metapathways build_db --test -d MPDB` after reinstalling MP to rebuild the reference indexes and mappings, and use a fresh output directory for the reviewer acceptance run.

## Run all three samples first

The fixture is included in both the repository and installed package; no separate input download is needed. Keep the prepared folder structure intact so manifest paths resolve correctly.

```bash
cd ~/mp-reviewer
metapathways analysis_wf \
  --manifest cami-reviewer/all.tsv \
  -o all -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools \
  --threads 4 --memory '4 GB' --max_tasks 2

metapathways report -o all --serve --no-browser --port 8765
```

The analysis uses the references in your workspace’s `MPDB` directory. It requires neither a full SwissProt/UniRef download nor Pathway Tools.

The explorer starts at **Samples**. On a remote machine, use the [SSH tunnel instructions](reports-tutorial.md#view-a-remote-report-through-ssh) and open the token-bearing URL printed by the report command. The single- and two-sample commands below are additional checks; they are not prerequisites for the three-sample run.

## 2. Run the single-sample workflow

```bash
cd ~/mp-reviewer
metapathways analysis_wf \
  --manifest cami-reviewer/single.tsv \
  -o single -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools \
  --threads 4 --max_cpus 8 --memory '4 GB' --max_memory '16 GB'
```

This uses `Urogenital_22`. These explicit resource settings suit a small interface test; use resources actually available on your machine. They are not recommended memory estimates for full metagenomes. Append `--dryrun` first if you want to inspect the plan, then repeat without that flag to execute it. Do not add `--test` to `analysis_wf`: that mode belongs to the separate bundled installation check and selects its own inputs.

## 3. Run two distinct samples together

```bash
metapathways analysis_wf \
  --manifest cami-reviewer/pair.tsv \
  -o pair -d MPDB \
  --annotation_dbs swissprot_test \
  --rRNA_refdbs SILVA_SSU_test SILVA_LSU_test \
  --skip_ptools \
  --threads 4 --max_cpus 8 --memory '4 GB' --max_memory '16 GB'
```

This uses `Gastrointestinal_5` and `Skin_28`. Ready tasks can overlap within the CPU and memory budgets; dependencies still determine when each task may start. Repeat the same command once to check reuse: successful unchanged tasks should report `ALREADY_COMPUTED`. That repeat is not a fresh runtime benchmark.

`all.tsv` selects all three samples. The manifests use paths relative to the manifest, so you can move the entire `cami-reviewer` folder without editing them. The `sample_id` column controls output names independently of assembly filenames.

To exercise automatic discovery on all three samples, replace `--manifest cami-reviewer/pair.tsv` with `-i cami-reviewer/inputs`, and choose a new output such as `-o all-auto`. The discovery root contains only `assemblies/`, `reads/`, and `mag_maps/`. Manifests and provenance live outside it so strict discovery does not reject extra files.

## 4. Check outputs and explore tables

| Check | Single-sample location |
| --- | --- |
| Resolved sample ID and paired read paths | `single/inputs.resolved.tsv` |
| Functional/taxonomic annotation tables | `single/Urogenital_22/results/annotation_table/` |
| Contig and ORF abundance | `single/Urogenital_22/results/rpkm/` |
| Genome assignments and split inputs | `single/Urogenital_22/magsplitter/` |
| Task outcomes and resource records | `single/logs/analysis_wf/RUN_ID/` |
| Report and explorer | `single/reports/` |

`RUN_ID` means the generated invocation directory, not a literal directory name. Inspect `summary.json` and the linked logs. Required tasks must succeed; a report existing alone does not prove success. Both paired-read paths must be distinct, all three genome bins should have split inputs, and the pair report should contain both sample IDs. Pathway tables are absent by design with `--skip_ptools`.

```bash
metapathways report -o single --serve --no-rebuild
```

In `EDA_portal.html`, select a sample, search an annotation, follow an ORF to its related records, and export a filtered CSV. Check that sample/entity identifiers remain in the export. Stop the server with Ctrl-C and use `-o pair` to inspect the two-sample output. Follow [the report tutorial](reports-tutorial.md) for remote SSH access and table relationships.

## 5. Optional: include licensed Pathway Tools

Follow [the license and installer guide](pathway-tools.md#get-the-installer). Build and register your own image:

```bash
metapathways build_pt \
  -i ~/Downloads/pathway-tools-29.5-linux-64-tier1-install \
  -o ~/mp-reviewer/containers
```

Repeat either workflow command with a new output directory, remove `--skip_ptools`, and add `--taxprune --taxonomic_scope all`. MP uses the registered image. This adds community and per-bin PGDB inference/export. It does not require MetaCyc as a sequence annotation database: inference uses the licensed MetaCyc inside the image. See the Pathway Tools guide if you also want to build MetaCyc annotation references.

Tiny genome fragments and sparse annotation references may yield few or no pathways. Check task outcomes and exports rather than requiring a universal pathway count. A full-reference test is more informative biologically.

## Provenance and validation

The bundle's [README](https://github.com/hallamlab/MetaPathways/blob/HEAD/metapathways/regtests/cami_reviewer/README.md), `provenance.json`, and `validation.json` describe selection, original sample paths, retained contig regions, matching read counts, hashes, and checks performed. The preparation script is [scripts/prepare_cami_reviewer.py](https://github.com/hallamlab/MetaPathways/blob/HEAD/scripts/prepare_cami_reviewer.py). No Pathway Tools software or MetaCyc sequence/database material is included.

Input validation covers all manifests, automatic discovery, intact mate pairing, complete contig-to-genome maps and gene prediction on every retained contig. It does not substitute for the reviewer executing the full workflow above. See [benchmarking](benchmarking.md) before collecting publication statistics from the complete CAMI samples.

## Recorded three-sample validation

The 2026-10-02 local SwissProt reviewer run completed all 51 tasks successfully. The [audit record](validation/reviewer-2026-10-02.json) confirms all 393 CDS records were retained, paired read identifiers matched, all 421 gene/RNA abundance records and 9 contig records matched their source measurements, and taxonomy remained tied to each reference database and target. Each sample produced three CAMI genome bins. The 28 report placeholders were RNA features, not missing CDS annotations. Skin_28 had no rRNA queries, explaining its two BLAST empty-query warnings. This run used `--skip_ptools`; it does not establish successful PGDB construction or validate full-scale/HPC performance.
