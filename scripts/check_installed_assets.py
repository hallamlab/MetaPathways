#!/usr/bin/env python3
"""Check the installed MP payload outside the source checkout (no biology run)."""
import hashlib
import json
from pathlib import Path
import metapathways
from metapathways.analysis_workflow import read_manifest, validate

root = Path(metapathways.__file__).resolve().parent
for name in ('nextflow.py', 'nf_worker.py', 'nf_databases.py', 'analysis_workflow.py',
             'reporting.py', 'report_server.py', 'protein_taxonomy.py', 'pt_container.py', 'pt_sequences.py',
             'bin/fastal', 'bin/fastdb', 'bin/metacount'):
    path = root / name
    assert path.is_file() and path.stat().st_size, f'Missing package asset: {path}'
assets = root / 'report_assets'
assert list(assets.glob('*.html')) and list(assets.glob('*.js')), assets
fixture = root / 'regtests/cami_reviewer'
for path, digest in json.loads((fixture / 'provenance.json').read_text())['sha256'].items():
    assert hashlib.sha256((fixture / path).read_bytes()).hexdigest() == digest, path
for label, count in [('single', 1), ('pair', 2), ('all', 3)]:
    assert len(validate(read_manifest(fixture / f'{label}.tsv'))) == count, label
print(f'Installed workflow modules, binaries, report assets and CAMI fixtures verified: {root}')
