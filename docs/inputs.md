# Organize assemblies, reads and genome assignments

[Home](../README.md) · [Command cookbook](commands.md) · [Manifest reference](../README.md#custom-analysis-manifest)

## Decide what each sample means

A sample is one assembly with its own reads and optional contig-to-genome assignments. The assembly is required. Reads add abundance measurements. A genome map adds genome-level splitting and PGDBs. MP does not assemble reads or infer bins itself.

Choose sample IDs that begin with a letter and contain only letters, digits and underscores, for example `Oral_15`. Keep them stable between reruns. IDs such as `sample`, `Sample` and `SAMPLE` are distinct. Do not use reserved directory names such as `reports`, `logs`, `inputs`, `assemblies`, `reads` or `mag_maps`.

`run -i ASSEMBLY_DIRECTORY` annotates multiple files, but it is not the automatic per-sample read/map matcher. Use `analysis_wf` for that complete workflow. It requires nucleotide FASTA; standalone `run` also has a protein-input mode.

## Automatic discovery

```text
my_inputs/
├── assemblies/
│   ├── SampleA.fasta
│   └── SampleB.fasta.gz
├── reads/
│   ├── SampleA_R1.fastq.gz
│   ├── SampleA_R2.fastq.gz
│   └── SampleB_interleaved.fastq.gz
└── mag_maps/
    ├── SampleA.tsv
    └── SampleB.tsv
```

Use `.fa`, `.fna` or `.fasta` for assemblies, optionally `.gz`. For reads use `.fq` or `.fastq`, optionally `.gz`. The suffixes `_R1`, `_R2`, `_interleaved` and `_single` specify the read layout. Matching is by sample ID, not file ordering or fuzzy similarity. Maps end in `.tsv`.

```bash
metapathways analysis_wf -i my_inputs -o analysis -d /path/to/MPDB \
  --annotation_dbs swissprot metacyc \
  --taxprune --taxonomic_scope all \
  --threads 8 --max_cpus 32 --max_memory '64 GB'
```

Every sample needs one unambiguous layout and one map by default. MP stops on missing, unmatched or ambiguous inputs instead of guessing. Do not place both paired files and an interleaved file for the same sample in this discovery directory. Use `--no_reads` or `--no_mags` to omit that input type for the entire discovered dataset. `--no_mags` still allows community pathway inference; `--skip_ptools` skips all PGDB construction.

If your directories are already separate:

```bash
metapathways analysis_wf \
  -i /project/assemblies --reads_dir /project/reads --mag_maps_dir /project/maps \
  -o /project/analysis -d /project/MPDB \
  --taxprune --taxonomic_scope all
```

## Organize with links, without copying large reads

A symbolic link gives an existing file another path without duplicating its data. Use absolute targets so moving the working directory does not change what they point to:

```bash
mkdir -p my_inputs/assemblies my_inputs/reads my_inputs/mag_maps
ln -s /data/original/anonymous_gsa.fasta my_inputs/assemblies/SampleA.fasta
ln -s /data/original/anonymous_reads.fq.gz my_inputs/reads/SampleA_interleaved.fastq.gz
ln -s /data/original/contig_to_genome.tsv my_inputs/mag_maps/SampleA.tsv
```

The original files must remain readable throughout the workflow. The same assembly or read file cannot be reused under multiple manifest samples via hard/symbolic aliases; this catches accidental duplicated inputs. Independent reviewer fixture copies are intentional and documented separately.

## Custom manifest: keep every file where it is

Use a tab-separated file with this exact header:

```text
sample_id	assembly	read_layout	reads_1	reads_2	mag_map
```

The columns are real tabs, not the literal characters `\t`. A spreadsheet can save “tab-delimited text”; do not save an Excel workbook and rename it `.tsv`. Blank cells must remain blank, not `None` or `NA`.

| Column | What to put in it |
| --- | --- |
| `sample_id` | Your chosen stable sample name; independent of assembly basename |
| `assembly` | Path to that sample's nucleotide FASTA |
| `read_layout` | `paired`, `interleaved`, `single` or `none` |
| `reads_1` | R1, the interleaved file, or the single-end file; blank for `none` |
| `reads_2` | R2 only for `paired`; blank otherwise |
| `mag_map` | Headerless contig-to-bin TSV, or blank to omit bins for that sample |

Relative paths are resolved against the manifest's directory. Absolute paths make files easier to locate across terminal sessions; on a cluster, they must also work on compute nodes. Different samples can use different read layouts and omit different optional inputs through the manifest.

```bash
metapathways analysis_wf --manifest /project/samples.tsv \
  -o /project/analysis -d /project/MPDB \
  --annotation_dbs swissprot metacyc \
  --taxprune --taxonomic_scope all \
  --threads 8 --max_cpus 32 --max_memory '64 GB'
```

Do not combine `--manifest` with discovery flags, `-1`, `-2` or `--interleaved`. The manifest carries those choices per sample. The resolved manifest and planned entities are saved in the output root as `inputs.resolved.tsv` and `inputs.entities.json`.

## Contig-to-genome map format

MAGSplitter needs exactly two tab-separated columns and **no header**:

```text
contig_001	GenomeA
contig_002	GenomeA
contig_003	GenomeB
```

The first column is the original FASTA record ID: the token after `>` up to the first whitespace. Do not use MP's renamed `Sample-C123` IDs here. The second column is the bin ID. Contigs may be absent if unbinned, but every supplied contig must belong to the assembly and must occur at most once in the map. Periods in bin names become underscores; avoid names that collide after normalization. Use names beginning with a letter and consisting of letters, numbers or underscores; `community` and names containing `non_binned` are reserved.

A contig removed by assembly QC cannot contribute features downstream. A mapped genome with no usable annotated features may have no generated PGDB input. Therefore the number of genome IDs in a map is an input inventory, not a promise of that many successful PGDBs.

## CAMI ground-truth bins versus recovered MAGs

CAMI `gsa_mapping.tsv` files provide original contig-to-source-genome assignments, often with taxon, source-sequence and coordinate columns. Convert only the anonymous contig ID and genome ID columns to the two-column format; keep the original file and an ID lookup if bin names are normalized.

These maps are suitable for demonstrating MP's splitting and per-genome analysis features. Describe them as **ground-truth genome bins**, not experimentally recovered MAGs. They do not test a binning algorithm or realistic contamination/recovery errors. A source genome appearing in several samples creates several sample-specific bins, not several distinct reference genomes.

For custom maps, validate duplicate IDs and completeness explicitly before a large benchmark. MP's workflow validation checks supplied IDs against assemblies, but assignment correctness remains a property of your source data.
