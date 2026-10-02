# MAGs and pathway inference

## Build Pathway Tools once

Follow the [Pathway Tools licensing and installer guide](pathway-tools.md#get-the-installer) to obtain your licensed Linux x86-64 installer, then provide it to MP:

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

## MetaCyc from Pathway Tools

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

## Community pathways

```bash
metapathways ptools -o results/sample --taxprune --taxonomic_scope all
```

These examples use broad cellular-life taxonomic pruning. Choose the scope to match your analysis and record it in your methods; see [scope choices](pathway-tools.md#choose-a-taxonomic-scope).

Each SIF task gets private Pathway Tools data, home and temporary state, allowing concurrent isolated instances. Pathway Tools uses one CPU per PGDB. A community PGDB failure fails the command. Native Pathway Tools is still supported when no image is selected, but serialized to protect shared state. The legacy `--container` flag retains its original meaning and bypasses automatic SIF selection.

## Add MAG pathways without reannotating

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
