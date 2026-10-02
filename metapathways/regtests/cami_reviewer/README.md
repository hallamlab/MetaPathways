# Tiny CAMI II reviewer inputs

Use the [reviewer walkthrough](../../../docs/reviewer-test.md) for normal MP commands.

This bundle contains three different samples from the CAMI II human-associated short-read dataset: Urogenital_22, Gastrointestinal_5 and Skin_28. Each has three 50,000-base assembly regions from three source genomes, intact paired simulated reads, and a headerless contig-to-genome map. Total input size is about 2.4 MiB compressed. No Pathway Tools installer, image, MetaCyc reference or licensed PGDB is included.

- `single.tsv`: Urogenital_22.
- `pair.tsv`: Gastrointestinal_5 and Skin_28.
- `all.tsv`: all three.
- `inputs/`: strict automatic-discovery layout containing only assemblies, reads and genome maps.
- `provenance.json`: selection method, source files, original GSA mapping rows, retained coordinates, aligned-pair counts and SHA256 hashes.
- `validation.json`: input validation and predicted gene counts; the complete MP workflow has not been executed as part of preparing this bundle.
- `*.selection.log`: minimap2 selection diagnostics.

Manifest paths are relative to the manifest file. Copy the whole bundle when moving it. Original contig and read identifiers are preserved; contigs are cropped to their first 50,000 bases. Source genome IDs are prefixed with `CAMI_` and punctuation becomes underscores for MP-compatible entity names. Coordinates in `provenance.json` explicitly distinguish original source mapping intervals from retained contig-relative intervals.

For each sample, the preparation script ranks contigs at least 20 kb long by source read density, chooses one contig from each of three different genomes, and retains at most 50 kb per contig. It scans the first 250,000 original interleaved read pairs, aligns them to the retained regions with minimap2, and keeps at most 3,000 complete pairs with at least one primary alignment at MAPQ 20 or higher. Both original mates and qualities are retained even if only one mate maps. Every retained contig must have at least ten selected pairs aligning to it. The cap and selection deliberately bias coverage; this is interface-test data, not an accuracy or abundance benchmark. Three short genome fragments are not three complete MAGs.

Maintainers can reproduce the data from an MP manifest pointing at the original local CAMI downloads:

```bash
python scripts/prepare_cami_reviewer.py \
  --source-manifest /path/to/full-cami-manifest.tsv \
  --output /path/to/new-cami-reviewer \
  --minimap2 /path/to/minimap2
```

Run from the source checkout. The output directory must not already exist. Source assemblies must have a sibling `gsa_mapping.tsv`, and the selected sample read inputs must be the original interleaved FASTQs. Python and minimap2 are required. Temporary read subsets and alignments are cleaned up after preparation. Reproduction requires the original large CAMI files; reviewers use the bundled small files directly.

## Source attribution

Derived from the CAMI II multi-sample human microbiome dataset, [DOI 10.4126/FRL01-006425518](https://doi.org/10.4126/FRL01-006425518). Cite the CAMI II study, Meyer et al., *Critical Assessment of Metagenome Interpretation: the second round of challenges*, Nature Methods (2022), [DOI 10.1038/s41592-022-01431-4](https://doi.org/10.1038/s41592-022-01431-4). Cropping, read selection, mate separation and MP-format maps are the modifications made here. The data retain their source terms; the MP software license does not relicense third-party data.
