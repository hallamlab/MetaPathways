"""Stage the vitamin B5 teaching results without software/reference databases.

Run with the repository on PYTHONPATH. Original data is never modified.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('source', type=Path)
parser.add_argument('destination', type=Path)
args = parser.parse_args()
source = args.source.resolve()
dest = args.destination.resolve()
if dest.exists():
    parser.error('Choose a new destination; existing files are never overwritten.')
dest.mkdir(parents=True)
shutil.copytree(source/'cami-test', dest/'cami-test')
samples = ['Gastrointestinal_5', 'Skin_28', 'Urogenital_22']
for sample in samples:
    original = source/'all'/sample
    for folder in ['preprocessed', 'orf_prediction', 'genbank', 'run_statistics']:
        shutil.copytree(original/folder, dest/'all'/sample/folder)
    for folder in ['annotation_table', 'rRNA', 'tRNA', 'rpkm']:
        shutil.copytree(original/'results'/folder, dest/'all'/sample/'results'/folder)
    for path in [original/'magsplitter/contig_to_mag.tsv', *original.glob('results/pgdb/community/*_pwy*.tsv'), *original.glob('results/pgdb/MAGs/*/*_pwy*.tsv')]:
        target = dest/'all'/path.relative_to(source/'all')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
# Preserve the successful analysis record, not unrelated report invocations or tokens.
summary = source/'all/logs/analysis_wf/a6ed138e64294ad5876fb2eab85f283c/summary.json'
record = json.loads(summary.read_text())
assert record['status'] == 'SUCCESS'
assert all(t['status'] in ('SUCCESS', 'ALREADY_COMPUTED') for t in record['tasks'])
target = dest/'all'/summary.relative_to(source/'all')
target.parent.mkdir(parents=True)
shutil.copy2(summary, target)
provenance = dest/'provenance'
provenance.mkdir()
for path in (source/'MPDB/functional_categories').glob('*reldate.txt'):
    shutil.copy2(path, provenance/path.name)
shutil.copy2(source/'MPDB/functional_categories/MetaCyc_provenance.json', provenance/'MetaCyc_provenance.json')
(provenance/'README.txt').write_text('Historical paths in provenance refer to the original run; they are not prerequisites.\nAnalysis: MP 4.0.8; Pathway Tools/MetaCyc 29.5.\nCommand: metapathways analysis_wf --manifest cami-test/all.tsv -o all -d MPDB --annotation_dbs swissprot metacyc --threads 4 --memory "4 GB" --max_tasks 2\n')
(dest/'README.txt').write_text('''MetaPathways vitamin B5 teaching example, snapshot 1

Inputs: cropped simulated CAMI II assemblies, paired reads and genome maps.
See cami-test/README.md and provenance.json for selection and attribution.
Results: MP 4.0.8 with SwissProt, MetaCyc and Pathway Tools 29.5.
This is a teaching dataset, not a biological benchmark or complete MAG collection.

With MP 4.0.13 or later installed, run from this directory:
  metapathways report -o all --serve
No Pathway Tools or reference MPDB is required to explore existing results.

The bundle retains annotations, RNA evidence, abundance, sequence/GFF/GenBank
files, contig-to-MAG membership and inferred pathway/gene TSV exports.
It excludes the licensed installer/SIF, MPDB, PGDB databases/archives, diagnostic
workspaces, alignments, and workflow caches. Complete dataset downloads contain
only retained files; they do not reconstruct excluded products.
Original execution paths in provenance are historical, not local requirements.
The successful workflow summary includes tasks reused after an earlier attempt.

New organism-specific inference results were generated using Pathway Tools and
MetaCyc. They are not a redistribution of the MetaCyc reference database.
https://bioinformatics.ai.sri.com/ptools/licensing/ptools-academic-license.shtml
Reference software/data remain subject to their respective terms.
''')
from metapathways.reporting import build_report
build_report(dest/'all')
files = []
for path in sorted(dest.rglob('*')):
    if path.is_symlink():
        raise ValueError(f'Unexpected symlink: {path}')
    if path.is_file():
        files.append({'path':path.relative_to(dest).as_posix(), 'bytes':path.stat().st_size, 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(dest/'manifest.json').write_text(json.dumps({'snapshot':1,'analysis_mp_version':'4.0.8','report_mp_version':'4.0.13','files':files}, indent=2)+'\n')
print(f'Staged {len(files)} files at {dest}')
