# Prepare reference databases

For licensed pathway inference, [choose your Pathway Tools setup](pathway-tools.md#choose-your-setup-before-running-commands)
before building an image. If you also want MetaCyc protein annotations, use the
combined image-and-reference build after preparing public references, then select
MetaCyc when annotating. The bundled installation test skips this licensed setup.

For a minimal reference database:

```bash
metapathways build_db -d /data/MPDB --func swissprot -a fast
```

SILVA and supporting enzyme/taxonomy files are built alongside functional references. Select multiple functional databases with `--func swissprot cazy uniref50`; `uniref90` is also supported and substantially larger. `-a blast` builds BLAST indexes instead of FAST indexes. Use `--dryrun` to inspect the planned downloads and indexing jobs first.

The builder downloads public references from their configured upstream locations. These are not all version-pinned. eggNOG requires a local FASTA at `MPDB/functional/eggnog`; the old builder had no working eggNOG acquisition rule. Licensed MetaCyc reference acquisition is separate from this public builder. Existing compatible MPDB installations can be used directly with `-d`.

Choose annotation databases actually present in your MPDB. The FAST/BLAST choice in `run` must match its indexes. Keep reference release records and checksums; changing a database can change biological results.
