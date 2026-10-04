"""Strict input discovery and one shared DAG for complete multi-sample analyses."""
import argparse
import copy
import csv
import fcntl
import gzip
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sys

from metapathways import nextflow

FIELDS = ('sample_id', 'assembly', 'read_layout', 'reads_1', 'reads_2', 'mag_map')
GUIDE = ('https://hallamlab-metapathways.readthedocs.io/en/latest/analysis.html#analysis-wf-input-layout '
         'and https://hallamlab-metapathways.readthedocs.io/en/latest/analysis.html#custom-analysis-manifest')
FASTA = re.compile(r'(.+)\.(?:fasta|fna|fa)(?:\.gz)?$')
READS = re.compile(r'(.+)_(R1|R2|interleaved|single)\.(?:fastq|fq)(?:\.gz)?$')
SAFE = re.compile(r'[A-Za-z][A-Za-z0-9_]*\Z')
RESERVED = {'logs', 'reports', 'inputs', 'assemblies', 'reads', 'mag_maps'}


def fail(message):
    raise ValueError(f'{message}\nSee {GUIDE}')


def parser():
    from metapathways.pipeline import runParser
    p = runParser('analysis_wf')
    sub = next(a for a in p._actions if isinstance(a, argparse._SubParsersAction)).choices['analysis_wf']
    sub.description = 'Annotate multiple metagenomes, split MAGs, build PGDBs and generate a combined report.'
    for action in sub._actions:
        if action.dest in ('fwd_fastq', 'rev_fastq', 'interleaved', 'samples', 'test'):
            action.help = argparse.SUPPRESS
        elif action.dest == 'input_file':
            action.help = 'dataset root, or assemblies directory with --reads_dir/--mag_maps_dir; alternative: --manifest'

    sub.add_argument('--manifest', help='TSV with sample_id, assembly, read_layout, reads_1, reads_2, mag_map')
    sub.add_argument('--reads_dir', help='Flat reads directory; -i then names the assemblies directory')
    sub.add_argument('--mag_maps_dir', help='Flat contig-to-MAG map directory; -i then names the assemblies directory')
    sub.add_argument('--no_reads', action='store_true', help='Explicitly omit read mapping during automatic discovery')
    sub.add_argument('--no_mags', action='store_true', help='Explicitly omit MAG splitting during automatic discovery')
    sub.add_argument('--skip_ptools', action='store_true', help='Omit community and MAG PGDB construction')
    sub.add_argument('--compact_results', action='store_true',
                     help='Use task scratch, archive PGDBs/diagnostics, and remove completed sample intermediates')
    sub.add_argument('--scratch_dir', help='Worker-local scratch directory for compact mode [Slurm: SLURM_TMPDIR; local: system temporary directory]')
    sub.add_argument('--image', help='Pathway Tools SIF [registered by build_pt]')
    from metapathways.pt_taxonomy import add_taxonomy_options
    add_taxonomy_options(sub)
    sub.add_argument('--no_transport_inference', action='store_true', help='Disable TIP transport inference')
    sub.add_argument('--taxprune', action='store_true', help='Enable Pathway Tools taxonomic pruning')
    sub.add_argument('--ptools_memory', type=nextflow.memory, default=None, help='Optional PGDB memory override [same as --memory]')
    return p


def files(directory):
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        fail(f'Missing input directory: {directory}')
    result = []
    for p in sorted(directory.iterdir()):
        if p.name.startswith('.'):
            continue
        if not p.is_file():
            fail(f'Expected a flat input directory; unexpected entry: {p}')
        result.append(p)
    return result


def discover(args):
    root = Path(args.input_file).expanduser().resolve()
    separate = bool(args.reads_dir or args.mag_maps_dir)
    assembly_dir = root if separate else root / 'assemblies'
    if not separate:
        if not root.is_dir():
            fail(f'Missing dataset directory: {root}')
        for p in root.iterdir():
            if not p.name.startswith('.') and p.name not in ('assemblies', 'reads', 'mag_maps'):
                fail(f'Unexpected dataset entry: {p}')
    rows = {}
    for p in files(assembly_dir):
        match = FASTA.fullmatch(p.name)
        if not match:
            fail(f'Unrecognized assembly filename: {p}')
        sample = match[1]
        if sample in rows:
            fail(f'Multiple assemblies match sample {sample}')
        rows[sample] = dict(zip(FIELDS, (sample, str(p), 'none', '', '', '')))
    if args.no_reads and args.reads_dir or args.no_mags and args.mag_maps_dir:
        fail('An input directory conflicts with its --no_reads/--no_mags flag')
    if not args.no_reads:
        groups = {}
        for p in files(args.reads_dir or root / 'reads'):
            match = READS.fullmatch(p.name)
            if not match or match[1] not in rows:
                fail(f'Unrecognized or unmatched reads: {p}')
            group = groups.setdefault(match[1], {})
            if match[2] in group:
                fail(f'Multiple read files match {match[1]} / {match[2]}')
            group[match[2]] = str(p)
        for sample, row in rows.items():
            group = groups.get(sample, {})
            if set(group) == {'R1', 'R2'}:
                row.update(read_layout='paired', reads_1=group['R1'], reads_2=group['R2'])
            elif set(group) in ({'single'}, {'interleaved'}):
                mode = next(iter(group))
                row.update(read_layout=mode, reads_1=group[mode])
            else:
                fail(f'Missing mates or ambiguous read layout for {sample}: {list(group)}')
    if not args.no_mags:
        for p in files(args.mag_maps_dir or root / 'mag_maps'):
            if p.suffix != '.tsv' or p.stem not in rows:
                fail(f'Unrecognized or unmatched MAG map: {p}')
            rows[p.stem]['mag_map'] = str(p)
        missing = [s for s, row in rows.items() if not row['mag_map']]
        if missing:
            fail(f'Missing MAG maps for: {", ".join(missing)}')
    return list(rows.values())


def read_manifest(filename):
    source = Path(filename).expanduser().resolve()
    with source.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != list(FIELDS):
            fail(f'Manifest header must be exactly: {", ".join(FIELDS)}')
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                fail(f'Manifest line {line} must contain exactly six tab-separated fields')
            for key in ('assembly', 'reads_1', 'reads_2', 'mag_map'):
                if row[key]:
                    p = Path(row[key]).expanduser()
                    row[key] = str((source.parent / p).resolve())
            rows.append(row)
    return rows


def map_entities(filename):
    entities, contigs = {}, set()
    with open(filename) as handle:
        for line, text in enumerate(handle, 1):
            fields = text.rstrip('\r\n').split('\t')
            if len(fields) != 2 or not all(fields) or any(x != x.strip() for x in fields):
                fail(f'{filename}:{line}: expected headerless contig_id<TAB>mag_id')
            contig, original = fields
            normalized = original.replace('.', '_')
            if not SAFE.fullmatch(normalized) or normalized == 'community' or 'non_binned' in normalized:
                fail(f'{filename}:{line}: unsupported/reserved MAG ID {original!r}')
            if normalized in entities and entities[normalized] != original:
                fail(f'{filename}: MAG IDs collide after replacing periods with underscores: {original}')
            if contig in contigs:
                fail(f'{filename}:{line}: repeated contig assignment {contig}')
            contigs.add(contig)
            entities[normalized] = original
    if not entities:
        fail(f'MAG map is empty: {filename}')
    return sorted(entities), contigs


def validate(rows):
    if not rows:
        fail('No samples were found')
    ids, used = set(), {}
    for row in rows:
        sample = row['sample_id']
        if not SAFE.fullmatch(sample) or sample.lower() in RESERVED or sample in ids:
            fail(f'Invalid, reserved, or duplicate sample ID: {sample!r}; use letters, digits and underscores, starting with a letter')
        ids.add(sample)
        mode, r1, r2 = row['read_layout'], row['reads_1'], row['reads_2']
        if mode not in ('paired', 'single', 'interleaved', 'none') or \
           (mode == 'paired' and not (r1 and r2)) or \
           (mode in ('single', 'interleaved') and (not r1 or r2)) or \
           (mode == 'none' and (r1 or r2)):
            fail(f'{sample}: read_layout and read paths disagree')
        for key in ('assembly', 'reads_1', 'reads_2', 'mag_map'):
            if not row[key]:
                if key == 'assembly':
                    fail(f'{sample}: missing assembly')
                continue
            p = Path(row[key]).expanduser().resolve()
            if not p.is_file() or p.stat().st_size == 0 or not os.access(p, os.R_OK):
                fail(f'{sample}: missing, empty or unreadable {key}: {p}')
            stat = p.stat()
            identity = (stat.st_dev, stat.st_ino)
            if identity in used:
                fail(f'{sample}: {key} reuses the same file as {used[identity]}')
            used[identity] = f'{sample}/{key}'
            row[key] = str(p)
        if not FASTA.fullmatch(Path(row['assembly']).name):
            fail(f'{sample}: assembly must end in .fa, .fna or .fasta, optionally .gz')
        row['entities'], mapped = map_entities(row['mag_map']) if row['mag_map'] else ([], set())
        # Check map membership before submitting any jobs. FASTQs are not fully scanned.
        found = set()
        opener = gzip.open if row['assembly'].endswith('.gz') else open
        with opener(row['assembly'], 'rt') as handle:
            first = handle.readline()
            if not first.startswith('>'):
                fail(f'{sample}: assembly is not FASTA')
            def header(text):
                tokens = text[1:].split()
                if not tokens:
                    fail(f'{sample}: empty FASTA identifier')
                if tokens[0] in found:
                    fail(f'{sample}: duplicate FASTA identifier {tokens[0]}')
                found.add(tokens[0])
            header(first)
            for text in handle:
                if text.startswith('>'):
                    header(text)
        unknown = mapped - found
        if unknown:
            fail(f'{sample}: MAG map contigs absent from assembly: {", ".join(sorted(unknown)[:5])}')
    return sorted(rows, key=lambda row: row['sample_id'])


def save_inputs(rows, output):
    text = io.StringIO()
    writer = csv.DictWriter(text, FIELDS, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    writer.writerows(rows)
    manifest = output / 'inputs.resolved.tsv'
    if manifest.exists() and manifest.read_text() != text.getvalue():
        fail(f'{manifest} describes different inputs; use a new output directory')
    entities_file = output / 'inputs.entities.json'
    entities = {row['sample_id']: row['entities'] for row in rows}
    if entities_file.exists() and json.loads(entities_file.read_text()) != entities:
        fail('The set of MAG IDs changed; use a new output directory')
    if not manifest.exists():
        # Refuse to combine an unrelated existing analysis with a new dataset.
        if any(p.is_dir() and not p.name.startswith('.') and p.name not in ('logs',) for p in output.iterdir()):
            fail('Use a new output directory for the first analysis_wf invocation')
        manifest.write_text(text.getvalue())
    if not entities_file.exists():
        entities_file.write_text(json.dumps(entities, indent=2) + '\n')
    staging = output / '.metapathways/analysis_wf/inputs'
    staging.mkdir(parents=True, exist_ok=True)
    staged = []
    for row in rows:
        item = dict(row)
        for key in ('assembly', 'reads_1', 'reads_2', 'mag_map'):
            if not row[key]:
                continue
            directory = staging / row['sample_id'] / key
            directory.mkdir(parents=True, exist_ok=True)
            suffix = ('.fasta' if key == 'assembly' else '.tsv' if key == 'mag_map' else '.fastq')
            if row[key].endswith('.gz'):
                suffix += '.gz'
            alias = directory / (row['sample_id'] + suffix)
            if alias.is_symlink():
                if alias.resolve() != Path(row[key]):
                    fail(f'Staged input changed unexpectedly: {alias}')
            elif alias.exists():
                fail(f'Unexpected file at staged input path: {alias}')
            else:
                alias.symlink_to(row[key])
            item[key] = str(alias)
        staged.append(item)
    return staged


def downstream(row, output, annotations, args, image):
    sample = row['sample_id']
    base = output / sample
    # Pathologic inputs follow annotation reports in the annotation DAG.
    parent = next(t['id'] for t in annotations if t.get('context', {}).get('name') == 'PATHOLOGIC_INPUT')
    tasks = []
    if row['mag_map']:
        ms = base / 'magsplitter'
        pf, orfmap = base / 'ptools/0.pf', base / 'ptools/orf_map.txt'
        table = base / f'results/annotation_table/{sample}.ORF_annotation_table.txt'
        feature_table = base / f'results/annotation_table/{sample}.ptinput.tsv'
        mapping = base / f'preprocessed/{sample}.mapping.txt'
        saved = ms / 'contig_to_mag.tsv'
        cmd = ['magsplitter', '-p', str(pf), '-r', str(orfmap), '-c', str(table),
               '-m', row['mag_map'], '-i', str(mapping), '-o', str(ms)]
        copy_cmd = [sys.executable, '-c', 'import shutil,sys; shutil.copyfile(*sys.argv[1:])', row['mag_map'], str(saved)]
        tasks.append(nextflow.task(f'{sample}:mag_split', f'{sample}:mag_split',
            [shlex.join(['mkdir', '-p', str(ms)]), shlex.join(cmd), shlex.join(copy_cmd)],
            [str(pf), str(orfmap), str(table), str(mapping), row['mag_map'], str(feature_table)],
            [str(ms / 'results'), str(saved)], [parent], memory=args.memory, sample=sample, adopt_existing=False,
            cache_version='authoritative-mag-coordinates-v1'))
    if args.skip_ptools:
        return tasks
    script = shutil.which('pgdb_build_wf.py') or str(Path(__file__).resolve().parents[1] / 'dev/pgdb_build_wf.py')
    if not Path(script).is_file():
        fail('pgdb_build_wf.py is missing; reinstall MetaPathways')
    for entity in ['community'] + row['entities']:
        community = entity == 'community'
        inputs = base / 'ptools' if community else base / 'magsplitter/results' / entity
        results = base / 'results/pgdb/community' if community else base / 'results/pgdb/MAGs' / entity
        tag = sample if community else entity
        cmd = [sys.executable, script, '--mp_out', str(base), '--tag', sample, '--entity', entity, '--image', image]
        if getattr(args, 'compact_results', False):
            cmd.append('--compact_results')
        if args.taxprune:
            cmd.append('--taxprune')
        if args.no_transport_inference:
            cmd.append('--no_transport_inference')
        from metapathways.pt_taxonomy import resolve_taxon
        taxon_id = resolve_taxon(args)
        if taxon_id is not None:
            cmd += ['--taxon_id', str(taxon_id)]
        tasks.append(nextflow.task(f'{sample}:pgdb:{entity}', f'{sample}:pgdb:{entity}', [shlex.join(cmd)],
            [image, str(inputs), str(base / f'results/annotation_table/{sample}.EC_RXN_map.tsv'),
             str(base / f'results/annotation_table/{sample}.ptinput.tsv'), str(base / f'preprocessed/{sample}.fasta')],
            [str(results / (tag + suffix)) for suffix in ('cyc.tar.bz2', '_pwy.tsv', '_pwy2orf.tsv')],
            [parent if community else f'{sample}:mag_split'], cpus=1, memory=args.ptools_memory or args.memory,
            sample=sample, entity=entity, allow_failure=not community, adopt_existing=False,
            cache_version='sequence-backed-pgdb-trna-names-v3',
            skip_if_missing=None if community else str(inputs / '0.pf')))
        tasks[-1]['fingerprint_inputs'] = tasks[-1]['inputs'] + [str(base / f'orf_prediction/{sample}.cds.gff')]
    return tasks


def main(argv=None):
    from metapathways.pipeline import prepare_annotation
    from metapathways.pt_container import registered_image
    p = parser()
    args = p.parse_args(['analysis_wf'] + list(sys.argv[1:] if argv is None else argv))
    if args.scratch_dir and not args.compact_results:
        fail('--scratch_dir requires --compact_results')
    if args.compact_results and any(getattr(args, k, None) for k in ('keep_work', 'work_dir', 'conda_cache')):
        fail('--compact_results cannot be combined with --keep_work, --work_dir or --conda_cache')
    if not args.output_dir or not args.refdb_dir or not (args.input_file or args.manifest):
        fail('Provide -o OUTPUT, -d MPDB and either -i INPUTS or --manifest FILE')
    if args.manifest and any((args.input_file, args.reads_dir, args.mag_maps_dir, args.no_reads, args.no_mags)):
        fail('--manifest cannot be combined with discovery flags')
    if args.fwd_fastq or args.rev_fastq or args.interleaved or args.test or args.samples:
        fail('Use per-sample reads and IDs from discovery/manifest, not -1/-2/--interleaved/--test/--samples')
    if args.input_format != 'fasta':
        fail('analysis_wf requires nucleotide FASTA assemblies')
    if any(getattr(args, key) == 'skip' for key in ('PREPROCESS_INPUT', 'ORF_PREDICTION', 'FILTER_AMINOS',
           'FUNC_SEARCH', 'PARSE_FUNC_SEARCH', 'ANNOTATE_ORFS', 'CREATE_ANNOT_REPORTS', 'PATHOLOGIC_INPUT')):
        fail('Required annotation stages cannot be skipped in analysis_wf; completed tasks are reused automatically')
    output = Path(args.output_dir).expanduser().resolve()
    refdb = Path(args.refdb_dir).expanduser().resolve()
    if any(re.search(r'[^A-Za-z0-9_./-]', str(x)) for x in (output, refdb)):
        fail('Output and MPDB paths must use letters, digits, underscores, hyphens, periods and slashes (legacy tool requirement)')
    from metapathways._version import __version__
    print(f'RUNNING MetaPathways: v{__version__}', flush=True)
    print(f'Output directory: {output}', flush=True)
    print('Validating sample inputs, read layouts and genome maps...', flush=True)
    rows = validate(read_manifest(args.manifest) if args.manifest else discover(args))
    from metapathways.compact_results import MARKER, resume_key, marker_state
    for row in rows:
        base = output/row['sample_id']
        if base.is_symlink():
            fail(f'Sample output must not be a symlink: {base}')
        if (base/MARKER).exists() and not args.compact_results:
            fail('This output contains compact samples; resume with --compact_results or use a new output directory')
        if args.compact_results:
            protected = [Path(r[k]).resolve() for r in rows
                         for k in ('assembly', 'reads_1', 'reads_2', 'mag_map') if r[k]]
            protected += [refdb]
            if args.image:
                protected.append(Path(args.image).expanduser().resolve())
            if any(p == base or base in p.parents for p in protected):
                fail(f'Compact mode requires inputs, MPDB and SIF outside sample output: {base}')
    image = None
    if not args.skip_ptools:
        image = args.image or registered_image()
        if not image or not Path(image).expanduser().is_file():
            fail('Build a Pathway Tools SIF with metapathways build_pt, pass --image, or explicitly --skip_ptools')
        image = str(Path(image).expanduser().resolve())
        if args.compact_results and any((output/r['sample_id']) in Path(image).parents for r in rows):
            fail('Compact mode requires the SIF outside sample output directories')
        if not shutil.which('apptainer'):
            fail('Apptainer is required for Pathway Tools')
    if any(row['mag_map'] for row in rows) and not shutil.which('magsplitter'):
        fail('MAGSplitter is required when MAG maps are supplied; see README installation instructions')
    output.mkdir(parents=True, exist_ok=True)
    control = output / '.metapathways/analysis_wf'
    control.mkdir(parents=True, exist_ok=True)
    with (control / 'planning.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fail(f'Another analysis_wf invocation is using {output}')
        print(f'Staging input links for {len(rows)} samples; references: {", ".join(args.annotation_dbs)}', flush=True)
        staged = save_inputs(rows, output)
        tasks = []
        for index, (original, row) in enumerate(zip(rows, staged), 1):
            sample = row['sample_id']
            prefix = f'Planning [{index}/{len(rows)}] {sample}'
            pgdbs = 'skipped' if args.skip_ptools else str(1 + len(row['entities']))
            print(f'{prefix}: reads={row["read_layout"]}; genome bins={len(row["entities"])}; '
                  f'PGDBs={pgdbs}; compact={"yes" if args.compact_results else "no"}', flush=True)
            first_task = len(tasks)
            key = resume_key(args, original) if args.compact_results else None
            state = marker_state(output/sample, key, args.force_redo) if args.compact_results else None
            if state == 'complete':
                print(f'{prefix}: completed compact results retained; no tasks scheduled', flush=True)
                continue
            if state == 'compacting':
                tasks.append(compact_task(sample, output, key, [], args.memory))
                print(f'{prefix}: resuming interrupted cleanup only', flush=True)
                continue
            sample_args = copy.copy(args)
            sample_args._analysis_planning = True
            sample_args.input_file = row['assembly']
            sample_args.output_dir, sample_args.refdb_dir = str(output), str(refdb)
            sample_args.fwd_fastq = row['reads_1'] or None
            sample_args.rev_fastq = row['reads_2'] or None
            sample_args.interleaved = row['read_layout'] == 'interleaved'
            annotation, _ = prepare_annotation(sample_args, p)
            # A combined analysis only reuses outputs whose receipts prove provenance.
            for task in annotation:
                task['adopt_existing'] = False
            tasks.extend(annotation)
            tasks.extend(downstream(row, output, annotation, args, image))
            if args.compact_results:
                dependencies = [t['id'] for t in tasks if t.get('sample') == sample]
                tasks.append(compact_task(sample, output, key, dependencies, args.memory))
            print(f'{prefix}: ready; {len(tasks) - first_task} tasks planned', flush=True)
        if args.force_redo:
            for task in tasks:
                if task['status'] != 'skip':
                    task['status'] = 'redo'
        print(f'Validated {len(rows)} samples. Resolved inputs: {output / "inputs.resolved.tsv"}', flush=True)
        print(f'Planning complete: {len(tasks)} tasks.', flush=True)
        if tasks:
            nextflow.launch(tasks, output, args, 'analysis_wf', dryrun=args.dryrun)
        if args.compact_results and not args.dryrun:
            # No sample needs staged input links or MP task receipts once complete.
            for name in ('inputs', 'receipts', 'tmp'):
                shutil.rmtree(control/name, ignore_errors=True)


def compact_task(sample, output, key, dependencies, memory):
    from metapathways.compact_results import MARKER
    command = shlex.join([sys.executable, '-m', 'metapathways.compact_results', str(output/sample), key])
    return nextflow.task(f'{sample}:compact_results', f'{sample}:compact_results', [command],
                         outputs=[str(output/sample/MARKER)], dependencies=dependencies,
                         sample=sample, memory=memory, adopt_existing=False)
