"""Opt-in, terminal per-sample cleanup for analysis_wf.

Keep report source files rather than a snapshot alone, so reports remain rebuildable.
Never follow directory symlinks when deleting generated output.
"""
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import sys

from metapathways.nf_worker import atomic_json

MARKER = 'compact-results.json'


def resume_key(args, row):
    options = vars(args).copy()
    # A compact sample cannot be silently reused under changed implementation.
    code = hashlib.sha256()
    for source in sorted(Path(__file__).parent.rglob('*.py')):
        code.update(source.relative_to(Path(__file__).parent).as_posix().encode())
        code.update(source.read_bytes())
    options['implementation_sha256'] = code.hexdigest()
    for key in ('dryrun', 'force_redo'):
        options.pop(key, None)
    paths = [Path(row[k]) for k in ('assembly', 'reads_1', 'reads_2', 'mag_map') if row[k]]
    from metapathways.pt_container import registered_image
    image = getattr(args, 'image', None) or (None if args.skip_ptools else registered_image())
    if image:
        paths.append(Path(image).expanduser())
    db = Path(args.refdb_dir).expanduser()
    for folder in ('functional', 'functional_categories', 'taxonomic', 'ncbi_tree'):
        paths.extend(p for p in (db/folder).rglob('*') if p.is_file())
    evidence = [(str(p.resolve()), p.stat().st_size, p.stat().st_mtime_ns) for p in sorted(set(paths))]
    payload = json.dumps([options, row, evidence], sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def marker_state(base, key, force=False):
    marker = base/MARKER
    if not marker.exists():
        return None
    record = json.loads(marker.read_text())
    if force or record.get('resume_key') != key:
        raise ValueError(f'{base} contains compact results; changed inputs/settings or --force_redo require a new output directory')
    if record.get('state') not in ('compacting', 'complete'):
        raise ValueError(f'Invalid compact-results marker: {marker}')
    for relative in record['keep']:
        if not (base/relative).is_file():
            raise ValueError(f'Compact result is missing {base/relative}; use a new output directory')
    return record['state']


def compact(base, key):
    from metapathways.reporting import Importer, SCHEMA
    base = Path(base)
    if base.is_symlink() or not base.is_dir():
        raise ValueError(f'Expected a real sample output directory: {base}')
    base = base.resolve()
    state = marker_state(base, key)
    if state == 'complete':
        return
    if state == 'compacting':
        record = json.loads((base/MARKER).read_text())
        keep = set(record['keep'])
    else:
        for pattern in ('preprocessed/*.mapping.txt', 'results/annotation_table/*.functional_and_taxonomic_table.txt'):
            if not list(base.glob(pattern)):
                raise ValueError(f'Missing required report source {pattern} in {base}; cleanup cancelled')
        # Validate all relationships and discover exactly what the explorer reads
        # before deleting anything. Use disk-backed SQLite for large samples.
        import tempfile
        with tempfile.TemporaryDirectory(prefix='mp-compact-', dir=base) as tmp:
            with sqlite3.connect(str(Path(tmp)/'check.sqlite')) as db:
                db.execute('PRAGMA journal_mode=MEMORY')
                db.execute('PRAGMA synchronous=OFF')
                db.execute('PRAGMA foreign_keys=ON')
                db.executescript(SCHEMA)
                importer = Importer(db, base)
                # Current-run summary is published only after the DAG finishes.
                # Read completed PGDB receipts to exclude failed partial tables.
                manifests = list((base.parent/'logs/analysis_wf').glob('*/tasks.json'))
                if manifests:
                    manifest = max(manifests, key=lambda p: p.stat().st_mtime_ns)
                    for task in json.loads(manifest.read_text()):
                        if task.get('sample') == base.name and task.get('entity'):
                            receipt = Path(task['invocation_receipt'])
                            status = json.loads(receipt.read_text())['status']
                            importer.pgdb_status[(base.name, task['entity'])] = status
                importer.sample(base)
                if db.execute('PRAGMA foreign_key_check').fetchone():
                    raise ValueError(f'Invalid report relationships in {base}; cleanup cancelled')
                keep = {p.relative_to(base).as_posix() for p in importer.source_stats}
        # Retain final supporting tables and diagnostics, including failed MAG logs.
        for directory, dirs, files in os.walk(base, followlinks=False):
            dirs[:] = [d for d in dirs if not (Path(directory)/d).is_symlink()]
            for name in files:
                p = Path(directory)/name
                rel = p.relative_to(base)
                if (rel.parts[0] in ('logs', 'run_statistics', 'reports') or
                    p.suffix == '.log' or name.endswith('_log.txt') or
                    (rel.parts[0] == 'results' and p.suffix in ('.tsv', '.txt'))):
                    keep.add(rel.as_posix())
        for relative in keep:
            p = base/relative
            if p.is_symlink() or not p.is_file():
                raise ValueError(f'Report source must be a regular file before cleanup: {p}')
        record = dict(state='compacting', resume_key=key, keep=sorted(keep), removed_bytes=0, removed_files=0)
        atomic_json(base/MARKER, record)
    keep.add(MARKER)
    for directory, dirs, files in os.walk(base, topdown=False, followlinks=False):
        for name in files:
            p = Path(directory)/name
            if p.relative_to(base).as_posix() not in keep:
                if not p.is_symlink():
                    record['removed_bytes'] += p.stat().st_size
                p.unlink()
                record['removed_files'] += 1
        for name in dirs:
            p = Path(directory)/name
            if p.is_symlink():
                p.unlink()
            elif not any(p.iterdir()):
                p.rmdir()
    record['state'] = 'complete'
    atomic_json(base/MARKER, record)
    print(f"Compact results: {base.name}; removed {record['removed_files']} files ({record['removed_bytes']} bytes)", flush=True)


if __name__ == '__main__':
    compact(Path(sys.argv[1]), sys.argv[2])
