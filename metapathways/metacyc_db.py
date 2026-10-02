"""Prepare MP references from a user's licensed MetaCyc export or Pathway Tools SIF."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import sys

from metapathways.pt_container import digest, exec_command, save_json


REQUIRED = ('protseq.fsa', 'proteins.dat', 'enzrxns.dat', 'reactions.dat',
            'pathways.dat', 'compounds.dat', 'classes.dat')
TABLES = ('MetaCyc-monomer-rxn-pairs.tsv', 'MetaCyc-PWY-RXN-CMP-map.tsv',
          'MetaCyc_PWY_Ontology.tsv')


def source_path(source):
    p = Path(source).expanduser().resolve()
    if p.name == 'protseq.fsa':
        p = p.parent
    if p.is_file() and p.suffix.lower() == '.sif':
        return p
    if p.is_dir():
        missing = [name for name in REQUIRED if not (p/name).is_file()]
        if not missing:
            return p
        raise ValueError('Incomplete MetaCyc data directory; missing ' + ', '.join(missing)
                         + '. Use the Pathway Tools SIF to export the matching flat files; protseq.fsa alone is insufficient. See README.md: MetaCyc from Pathway Tools.')
    raise ValueError(f'MetaCyc source must be a local SIF or complete data directory: {p}; copy remote SFTP data locally first')


def make_tables(source, output):
    """Use the repository's MetaCyc builders on a private, same-release copy."""
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fasta_ids = []
    with (source/'protseq.fsa').open() as stream:
        for line in stream:
            if line.startswith('>'):
                match = re.match(r'>gnl\|META\|([^\s]+)', line)
                if not match:
                    raise ValueError('Unexpected MetaCyc FASTA identifier: ' + line.strip())
                fasta_ids.append(match.group(1))
    fasta_set = set(fasta_ids)
    if not fasta_ids or len(fasta_ids) != len(fasta_set):
        raise ValueError('MetaCyc FASTA has no sequences or contains duplicate identifiers')
    versions = set()
    identifiers = {}
    for filename in REQUIRED[1:]:
        with (source/filename).open(errors='replace') as stream:
            ids = set()
            for line in stream:
                if line.startswith('# Version:'):
                    versions.add(line.partition(':')[2].strip())
                elif line.startswith('UNIQUE-ID - '):
                    ids.add(line.removeprefix('UNIQUE-ID - ').strip())
        if not ids:
            raise ValueError('Empty MetaCyc flat file: ' + filename)
        identifiers[filename] = ids
    if len(versions) > 1:
        raise ValueError('Mixed MetaCyc releases in input flat files: ' + ', '.join(sorted(versions)))
    missing = fasta_set - identifiers['proteins.dat']
    if missing:
        raise ValueError('FASTA proteins missing from proteins.dat: ' + ', '.join(sorted(missing)[:10]))
    scripts = Path(__file__).parent / 'build_DBs'
    subprocess.run([sys.executable, str(scripts/'metacyc_mapping_build.py'), str(source)], check=True)
    subprocess.run([sys.executable, str(scripts/'metacyc_build_ont.py'), str(source), str(output)], check=True)
    for name in TABLES[:2]:
        shutil.copy2(source/name, output/name)
    counts = {}
    for name in TABLES:
        with (output/name).open() as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            counts[name] = sum(1 for row in reader)
        if not counts[name]:
            raise ValueError('MetaCyc builder produced an empty table: ' + name)
    with (output/TABLES[0]).open() as stream:
        pairs = list(csv.DictReader(stream, delimiter='\t'))
    if any(r['MC'] not in fasta_set or r['RXN'] not in identifiers['reactions.dat'] for r in pairs):
        raise ValueError('MetaCyc mapping table references proteins or reactions absent from its source')
    counts['proteins'] = len(fasta_ids)
    counts['proteins_with_reaction'] = len({r['MC'] for r in pairs})
    counts['proteins_without_reaction'] = len(fasta_ids) - counts['proteins_with_reaction']
    return counts


def export_image(image, state):
    state = Path(state)
    (state/'export').mkdir(parents=True)
    script = r'''set -eu
cp -a /opt/ptools-local-template /data/ptools-local
if test -f /opt/mp-ncbirc; then cp /opt/mp-ncbirc /data/.ncbirc; fi
meta_root=/opt/pathway-tools/aic-export/pgdbs/biocyc/metacyc
version=$(cat "$meta_root/default-version")
test -n "$version"
cp "$meta_root/$version/data/protseq.fsa" /data/export/protseq.fsa
printf '%s\n' "$version" > /data/export/version.txt
xvfb-run -a /opt/pathway-tools/pathway-tools -no-patch-download -no-cel-overview -disable-metadata-saving -nologfile -lisp -eval '(progn (with-organism (:org-id '\''META) (dump-frames-to-attribute-value-files "/data/export/")) (format t "~%MP-METACYC-EXPORTED~%") (exit))'
'''
    # This exports the bundled reference; it does not create a sample PGDB.
    log = state/'export.log'
    with log.open('w') as out:
        result = subprocess.run(exec_command(image, state, ['sh', '-ec', script]),
                                stdout=out, stderr=subprocess.STDOUT, timeout=1800)
    text = log.read_text(errors='replace')
    print(text, flush=True)
    if result.returncode or 'MP-METACYC-EXPORTED' not in text.splitlines():
        raise RuntimeError('MetaCyc reference export failed; see ' + str(log))
    return source_path(state/'export')


def prepare(source, root, aligner):
    source, root = source_path(source), Path(root).resolve()
    if aligner not in ('fast', 'blast'):
        raise ValueError('MetaCyc aligner must be fast or blast')
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.metacyc-build-', dir=root) as scratch:
        work = Path(scratch)
        if source.is_file():
            raw = export_image(source, work/'state')
        else:
            # Existing scripts write tables beside the FASTA; never modify the source.
            raw = work/'source'
            raw.mkdir()
            for name in REQUIRED:
                shutil.copy2(source/name, raw/name)
            (raw/'version.txt').write_text((source/'version.txt').read_text() if (source/'version.txt').is_file() else source.parent.name)

        staged = work/'reference'
        functional, categories = staged/'functional', staged/'functional_categories'
        formatted = functional/'formatted'
        formatted.mkdir(parents=True)
        categories.mkdir()
        counts = make_tables(raw, categories)
        fasta = functional/'metacyc'
        shutil.copy2(raw/'protseq.fsa', fasta)
        prefix = formatted/'metacyc'
        with Path(str(prefix)+'-names.txt').open('w') as out, fasta.open() as stream:
            for line in stream:
                if line.startswith('>'):
                    out.write(line)
        command = (['fastdb', '-p', str(prefix), str(fasta)] if aligner == 'fast' else
                   ['makeblastdb', '-in', str(fasta), '-dbtype', 'prot', '-parse_seqids', '-out', str(prefix)])
        subprocess.run(command, check=True)
        sentinel = Path(str(prefix) + ('.prj' if aligner == 'fast' else '.pdb'))
        if not sentinel.is_file():
            raise RuntimeError('MetaCyc indexer did not produce ' + sentinel.name)
        version_file = raw/'version.txt'
        version = version_file.read_text().strip() if version_file.is_file() else raw.parent.name
        manifest = dict(source=str(source), release=version, aligner=aligner, counts=counts,
                        created_at=datetime.now(timezone.utc).isoformat(),
                        source_sha256={name: digest(raw/name) for name in REQUIRED})
        if source.is_file():
            manifest['image_sha256'] = digest(source)
        (categories/'MetaCyc_reldate.txt').write_text('Release: ' + version + '\n')
        save_json(categories/'MetaCyc_provenance.json', manifest)
        # All export, mapping and indexing checks finish before publication.
        new_indexes = {file.name for file in formatted.glob('metacyc.*')}
        for old in (root/'functional/formatted').glob('metacyc.*'):
            if old.name not in new_indexes:
                old.unlink()
        for file in sorted(staged.rglob('*')):
            if file.is_file():
                destination = root/file.relative_to(staged)
                destination.parent.mkdir(parents=True, exist_ok=True)
                file.replace(destination)
        print('Prepared MetaCyc ' + version + ': ' + json.dumps(counts), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--aligner', choices=('fast', 'blast'), required=True)
    args = parser.parse_args()
    prepare(args.source, args.output, args.aligner)
