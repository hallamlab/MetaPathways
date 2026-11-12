# Maintain the vitamin B5 vignette

The user walkthrough is `docs/vitamin-b5-vignette.md`. Its biological snapshot
was produced with MP 4.0.8; reporting/screenshots use 4.0.13. Do not relabel an
old analysis as a new one. Input provenance and reference release records travel
with the downloadable archive.

To reproduce from the author's private `mp-demo` folder, choose a new staging
path outside the source. The preparation script refuses existing destinations:

```bash
PYTHONPATH=. python scripts/vignette/prepare_example.py /path/to/mp-demo /tmp/mp-vignette-v1
python scripts/vignette/capture.py /tmp/mp-vignette-v1 --browser /path/to/chromium
cd /tmp/mp-vignette-v1/exports
python /path/to/MetaPathways/docs/assets/vignette/summarize_evidence.py
cd /path/to/MetaPathways
python scripts/vignette/pack_example.py /tmp/mp-vignette-v1
```

Capture requires Playwright and Chromium, plus dependencies for the MP report
server. It serves localhost only and saves screenshots without the address bar
or access token. It checks the two browser CSV exports and a MAG ZIP containing
FASTA and GenBank. The packer checks database integrity, all indexed source
paths, key biological findings, and exported row counts before making checksums.
Extract the archive somewhere else, verify manifest hashes, and rebuild its
report before publishing. Compare the walkthrough's clicks and selected columns
with the actual UI, not just SQLite queries.

`terminal-steps.json` distinguishes command references from historical evidence
and actual downstream output. These are rendered terminal transcript/reference
images, not screenshots of a new Mamba/Pathway Tools installation. Never invent
success messages or present old output as newly executed. Explorer images are
browser captures of the staged biological results. Update expected counts only
after investigating a real change. Preserve the original dataset.

The allowlisted bundle excludes licensed installers, SIFs, MPDB, PGDB databases
and archives, diagnostic workspaces, BAMs and caches. It includes derived
pathway/gene tables and preserves source attribution. Review applicable upstream
terms before adding any new reference or licensed material. The public SRI
academic agreement distinguishes user-created PGDBs from reference databases:
https://bioinformatics.ai.sri.com/ptools/licensing/ptools-academic-license.shtml

Upload the archive/checksum as new assets on the existing v4.0.13 GitHub release;
never overwrite a published snapshot. For changed data, choose a new snapshot
name and update the documentation URLs. The docs change alone does not require
rebuilding the Conda package or containers.

Validate documentation with:

```bash
python -m sphinx -W --keep-going -b html docs /tmp/mp-vignette-html
```
