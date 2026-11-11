"""Preflight and restart records for explicit MPDB create/add/update operations."""
import hashlib
import json
from pathlib import Path
import shlex

SUPPORT = ('functional_categories/enzyme.dat', 'ncbi_tree/ncbi_taxonomy_tree.txt')
SILVA = ('SILVA_138.1_LSURef_NR99_tax_silva_trunc', 'SILVA_138.1_SSURef_NR99_tax_silva_trunc')


def support_ready(root):
    root = Path(root)
    paths = list(SUPPORT)
    paths += [f'taxonomic/formatted/{db}{suffix}' for db in SILVA for suffix in ('.ndb', '-names.txt')]
    return all((root / name).is_file() and (root / name).stat().st_size for name in paths)


def references(root):
    folder = Path(root) / 'functional'
    names = {p.name for p in folder.glob('*') if p.is_file()}
    names.update(p.name.rsplit('.', 1)[0] for p in (folder/'formatted').glob('*')
                 if p.suffix in ('.prj', '.pdb'))
    return names


def signature(tasks, args):
    from metapathways.nf_worker import fingerprint
    commands = []
    sources = []
    for task in tasks:
        values = task['commands']
        if task['id'] == 'screen_metacyc':
            # Runtime capacity can change between attempts; it is not a new DB selection.
            values = [shlex.join(shlex.split(c)[:-2]) for c in values]
        if task['id'] == 'prepare_metacyc':
            sources = fingerprint(task['inputs'])
        commands.append((task['id'], values))
    payload = {'mode': 'add' if args.add else 'update' if args.update else 'create',
               'databases': sorted(args.func), 'aligner': args.aligner,
               'commands': commands, 'source_identity': sources}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def preflight(root, args, token):
    root = Path(root).expanduser().resolve()
    if root.exists() and not root.is_dir():
        raise ValueError(f'MPDB destination is not a directory: {root}')
    record = root / '.metapathways/database-build.json'
    previous = json.loads(record.read_text()) if record.exists() else None
    if previous and previous['signature'] == token and (previous['status'] == 'running' or not (args.add or args.update)):
        return record  # Exact command repeats reuse task receipts, including partial builds.
    if previous and previous['status'] == 'running':
        raise ValueError('This MPDB has an unfinished build. Repeat its original command before changing references: '
                         + previous['command'])
    selected = set(args.func)
    found = references(root)
    command = ['metapathways', 'build_db', '-d', str(root), '--func', *args.func, '-a', args.aligner]
    for flag in ('metacyc_source', 'screen_image'):
        if getattr(args, flag, None):
            command += ['--' + flag, str(getattr(args, flag))]
    if getattr(args, 'skip_pt_screen', False):
        command += ['--skip_pt_screen']
    def suggestion(names, operation):
        tail = command[5 + len(args.func):]
        if 'metacyc' not in names:
            tail = tail[:2]  # Only the index format applies to public references.
        return shlex.join(command[:5] + sorted(names) + tail + [operation])

    if args.test and (args.add or args.update):
        raise ValueError('--test cannot be combined with --add or --update; use a separate test workspace.')
    if args.add or args.update:
        if not support_ready(root):
            raise ValueError('Existing MPDB is missing shared support files. Resume its original build, or choose a new directory for a complete build.')
        if args.add and selected & found:
            raise ValueError('Already present: ' + ', '.join(sorted(selected & found)) +
                             '. To replace selected references explicitly, use: ' + suggestion(selected & found, '--update') +
                             ('. Add the remaining references separately: ' + suggestion(selected - found, '--add') if selected - found else ''))
        if args.update and selected - found:
            raise ValueError('Not present: ' + ', '.join(sorted(selected - found)) +
                             '. To add references, use: ' + suggestion(selected - found, '--add'))
    else:
        # CLI logging alone does not turn an empty directory into an MPDB.
        entries = {p.name for p in root.iterdir()} - {'logs', '.metapathways'} if root.exists() else set()
        # prepare_test intentionally supplies only raw seeds before its first build.
        if args.test:
            from metapathways.test_data import SEEDS
            allowed = set(SEEDS)
            files = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()
                     and p.relative_to(root).parts[0] not in ('logs', '.metapathways')}
            if files <= allowed:
                entries = set()
        # eggNOG is supplied locally, rather than downloaded.
        if entries == {'functional'} and 'eggnog' in selected and found == {'eggnog'} and list((root/'functional').iterdir()) == [root/'functional/eggnog']:
            entries = set()
        if entries:
            if support_ready(root):
                option = '--update' if selected <= found else '--add'
                raise ValueError('Existing MPDB contains: ' + ', '.join(sorted(found)) +
                                 '. No references changed. Use an explicit operation: ' + suggestion(selected if option == '--update' else selected - found, option))
            raise ValueError(f'Destination is nonempty and is not a complete MPDB: {root}. Choose a new directory, or resume its original build.')
    return record


def execute(root, args, tasks, command, launch):
    """Keep restart identity and a nonblocking writer lock for the whole build."""
    import fcntl
    import uuid
    from metapathways.nf_worker import atomic_json
    token = signature(tasks, args)
    record = preflight(root, args, token)
    if args.dryrun:
        prepare_tasks(tasks, args, 'preview')
        return launch(tasks, root, args, 'build_db', dryrun=True)
    record.parent.mkdir(parents=True, exist_ok=True)
    with (record.parent / 'database-build.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another database builder is using this MPDB. Wait for it to finish.') from None
        preflight(root, args, token)
        previous = json.loads(record.read_text()) if record.exists() else {}
        generation = (previous.get('generation') if previous.get('status') == 'running'
                      and previous.get('signature') == token else None) or uuid.uuid4().hex
        state = dict(signature=token, status='running', command=command, generation=generation)
        prepare_tasks(tasks, args, generation)
        atomic_json(record, state)
        launch(tasks, root, args, 'build_db', dryrun=False)
        state['status'] = 'complete'
        atomic_json(record, state)


def prepare_tasks(tasks, args, generation):
    if args.update:
        print('Updating selected references. Stop analyses reading this MPDB until the build finishes.')
        for task in tasks:
            if task['id'] != 'directories':
                task['cache_version'] = generation
                task['adopt_existing'] = False
                if task['id'].startswith('index_'):
                    prefix = str(task['outputs'][0]).rsplit('.', 1)[0]
                    task['commands'].insert(0, 'rm -f ' + shlex.quote(prefix) + '.*')
