"""Nextflow scheduling behind the existing MetaPathways CLI.

Stage tools continue writing their established output paths. Each Nextflow task
validates a durable receipt against those files instead of trusting a cached
sentinel when a user has removed or changed published outputs.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import uuid
import queue
import signal
import threading
import time


def positive(value):
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError('must be a positive integer')
    return n


def memory(value):
    if not re.fullmatch(r'\d+(?:\.\d+)?\s*(?:KB|MB|GB|TB)', value, re.I):
        raise argparse.ArgumentTypeError('use a size such as "16 GB"')
    if float(re.match(r'[0-9.]+', value).group()) <= 0:
        raise argparse.ArgumentTypeError('memory must be positive')
    return value.upper()


def add_resources(parser, task_memory='16 GB'):
    group = parser.add_argument_group('Execution Resources')
    group.add_argument('--max_cpus', type=positive, default=None,
                       help='total CPU budget [local: available CPUs; Slurm: no aggregate cap]')
    group.add_argument('--memory', type=memory, default=task_memory,
                       help=f'memory reservation per task [{task_memory}]')
    group.add_argument('--max_memory', type=memory, default=None,
                       help='total memory budget [local: available memory; Slurm: no aggregate cap]')
    group.add_argument('--max_tasks', type=positive, default=None,
                       help='maximum submitted tasks, including queued/running [local: CPU budget; Slurm: 4]')


    group.add_argument('--executor', choices=['local', 'slurm'], default='local',
                       help='execution backend [local]; Slurm uses your logged-in cluster identity')
    group.add_argument('--account', type=slurm_token, help='Slurm allocation/account')
    group.add_argument('--partition', type=slurm_token, help='Slurm partition [cluster default]')
    group.add_argument('--qos', type=slurm_token, help='Slurm quality of service')
    group.add_argument('--reservation', type=slurm_token, help='Slurm reservation')
    group.add_argument('--time_limit', type=duration, default='24h', help='Slurm walltime per task [24h]')
    group.add_argument('--submit_rate', type=positive, default=6, help='maximum Slurm submissions per minute [6]')
    group.add_argument('--work_dir', help='custom Nextflow work directory; retained after the run')
    group.add_argument('--conda_cache', help='custom Conda cache directory; retained after the run')
    group.add_argument('--keep_work', action='store_true', help='retain automatically allocated work/cache directories')


def slurm_token(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise argparse.ArgumentTypeError('use letters, digits, underscore, period or hyphen')
    return value


def duration(value):
    if not re.fullmatch(r'[1-9][0-9]*(?:s|m|h|d)', value):
        raise argparse.ArgumentTypeError('use a positive duration such as 30m, 12h or 2d')
    return value


def memory_bytes(value):
    number, unit = re.fullmatch(r'([0-9.]+)\s*(KB|MB|GB|TB)', memory(value)).groups()
    return int(float(number) * 1024 ** {'KB': 1, 'MB': 2, 'GB': 3, 'TB': 4}[unit])


def ordered(tasks):
    by_id = {t['id']: t for t in tasks}
    if len(by_id) != len(tasks):
        raise ValueError('Duplicate task identifiers')
    result, visiting, seen = [], set(), set()
    def visit(identifier):
        if identifier not in by_id:
            raise ValueError(f'Unknown dependency: {identifier}')
        if identifier in visiting:
            raise ValueError(f'Cyclic task dependency: {identifier}')
        if identifier in seen:
            return
        visiting.add(identifier)
        for dep in by_id[identifier]['dependencies']:
            visit(dep)
        visiting.remove(identifier)
        seen.add(identifier)
        result.append(by_id[identifier])
    for t in tasks:
        visit(t['id'])
    return result


def task(identifier, label, commands, inputs=(), outputs=(), dependencies=(), cpus=1,
         memory='16 GB', **kwargs):
    return dict(id=identifier, label=label, commands=commands, inputs=list(inputs),
                outputs=list(outputs), dependencies=list(dependencies), cpus=cpus,
                memory=memory, status='yes', **kwargs)


def groovy(value):
    return "'" + str(value).replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n') + "'"


def process_definition(t, name, manifest):
    package_root = str(Path(__file__).resolve().parent.parent)
    command = shlex.join([sys.executable, '-m', 'metapathways.nf_worker', str(manifest), t['id']])
    body = f'export PYTHONPATH={shlex.quote(package_root)}${{PYTHONPATH:+:$PYTHONPATH}}\n{command}\n'
    count = len(t['dependencies'])
    inputs = '\n'.join(f'    val dependency_{i}' for i in range(1 if count > 64 else max(1, count)))
    return f'''process {name} {{
    tag {groovy(t['label'])}
    cpus {t['cpus']}
    memory {groovy(t['memory'])}
    cache false
    input:
{inputs}
    output:
    path 'receipt.json'
    script:
    {groovy(body)}
}}'''


def render_modules(tasks, manifest, batch_size=64):
    """Bound compiled script size without adding scheduling barriers.

    Each module exposes individual task completion channels. A consumer waits
    only for its own dependencies, even when their producers share a module
    with other tasks. All processes still use the same Nextflow executor.
    """
    if batch_size < 1 or batch_size > 64:
        raise ValueError('Workflow module size must be between 1 and 64')
    tasks = ordered(tasks)
    # Sample completion/cleanup has a large fan-in. Flattening every sample's
    # channels into main.nf can exceed the JVM's 64 KiB method limit even though
    # process definitions are batched. Keep independent sample DAGs in named
    # subworkflows so both process code AND channel wiring stay bounded.
    samples = list(dict.fromkeys(t.get('sample') for t in tasks))
    by_id = {t['id']: t for t in tasks}
    if len(samples) > 1 and None not in samples and all(
            by_id[dep].get('sample') == t['sample'] for t in tasks for dep in t['dependencies']):
        files, includes, calls = {}, [], []
        for i, sample in enumerate(samples):
            name = f'SAMPLE_{i:04d}'
            subset = [t for t in tasks if t['sample'] == sample]
            for relative, content in render_modules(subset, manifest, batch_size).items():
                if relative == 'main.nf':
                    content = content.replace('\nworkflow {\n', f'\nworkflow {name} {{\n')
                files[f'samples/{name}/{relative}'] = content
            includes.append(f"include {{ {name} }} from './samples/{name}/main'")
            calls.append(f'    {name}()')
        files['main.nf'] = 'nextflow.enable.dsl=2\n\n' + '\n'.join(includes) + '\n\nworkflow {\n' + '\n'.join(calls) + '\n}\n'
        return files
    names = {t['id']: f'TASK_{i:04d}' for i, t in enumerate(tasks)}
    batches = [tasks[i:i + batch_size] for i in range(0, len(tasks), batch_size)]
    owner = {t['id']: i for i, batch in enumerate(batches) for t in batch}
    exports = {dep for t in tasks for dep in t['dependencies'] if owner[dep] != owner[t['id']]}
    files, includes, calls = {}, [], []
    for i, batch in enumerate(batches):
        workflow = f'BATCH_{i:04d}'
        external = list(dict.fromkeys(dep for t in batch for dep in t['dependencies'] if owner[dep] != i))
        bundled = len(external) > 64
        ports = {dep: ('upstream.' if bundled else '') + f'upstream_{j}' for j, dep in enumerate(external)}
        blocks = [process_definition(t, names[t['id']], manifest) for t in batch]
        lines = [f'workflow {workflow} {{']
        if external:
            lines += ['    take:'] + (['    upstream'] if bundled else [f'    {ports[dep]}' for dep in external])
        lines.append('    main:')
        for t in batch:
            channels = [ports[dep] if owner[dep] != i else names[dep] + '.out' for dep in t['dependencies']]
            args = ', '.join(channels)
            if len(channels) > 64:
                # One completion gate, with bounded operator argument counts.
                # collect emits only after every prerequisite channel closes.
                args = 'Channel.empty()'
                for start in range(0, len(channels), 32):
                    args += '.mix(' + ', '.join(channels[start:start+32]) + ')'
                args += '.collect()'
            lines.append(f"    {names[t['id']]}({args or 'Channel.value(true)'})")
        emitted = [t['id'] for t in batch if t['id'] in exports]
        if emitted:
            lines += ['    emit:'] + [f'    done_{names[dep]} = {names[dep]}.out' for dep in emitted]
        lines.append('}')
        files[f'modules/{workflow}.nf'] = '\n\n'.join(blocks) + '\n\n' + '\n'.join(lines) + '\n'
        includes.append(f"include {{ {workflow} }} from './modules/{workflow}'")
        args = ', '.join(f'BATCH_{owner[dep]:04d}.out.done_{names[dep]}' for dep in external)
        if bundled:
            args = '[' + ', '.join(f'upstream_{j}: BATCH_{owner[dep]:04d}.out.done_{names[dep]}'
                                  for j, dep in enumerate(external)) + ']'
        calls.append(f'    {workflow}({args})')
    files['main.nf'] = 'nextflow.enable.dsl=2\n\n' + '\n'.join(includes) + '\n\nworkflow {\n' + '\n'.join(calls) + '\n}\n'
    return files


def local_capacity():
    cpus = len(os.sched_getaffinity(0))
    available = None
    try:
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith('MemAvailable:'):
                available = int(line.split()[1]) * 1024
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max':
            cpus = min(cpus, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        pass
    try:
        limit = Path('/sys/fs/cgroup/memory.max').read_text().strip()
        if limit != 'max':
            remaining = max(1, int(limit) - int(Path('/sys/fs/cgroup/memory.current').read_text()))
            available = min(available, remaining) if available is not None else remaining
    except (OSError, ValueError):
        pass
    return cpus, (f'{max(1, available // (1024 * 1024))} MB' if available is not None else None)


def configuration(tasks, args, conda_cache):
    backend = getattr(args, 'executor', 'local')
    local_cpus, local_memory = local_capacity() if backend == 'local' else (None, None)
    cpus = getattr(args, 'max_cpus', None) or local_cpus
    limit_memory = getattr(args, 'max_memory', None) or local_memory
    if any(t['cpus'] < 1 or (cpus is not None and t['cpus'] > cpus) for t in tasks):
        raise ValueError('Task threads must be positive and no greater than --max_cpus')
    largest_memory = max(memory_bytes(t['memory']) for t in tasks)
    if limit_memory and largest_memory > memory_bytes(limit_memory):
        raise ValueError('A task requests more memory than --max_memory')
    config = [f'process.executor = {groovy(backend)}', "process.errorStrategy = 'finish'",
              'process.maxRetries = 0', f'conda.cacheDir = {groovy(str(conda_cache / "envs"))}']
    if backend == 'local':
        if any(getattr(args, key, None) for key in ('account', 'partition', 'qos', 'reservation')):
            raise ValueError('Slurm flags require --executor slurm')
        config.append(f'executor.cpus = {cpus}')
        if limit_memory:
            config.append(f'executor.memory = {groovy(limit_memory)}')
        # Default Nextflow queueSize (100) can otherwise underutilize large hosts.
        jobs = getattr(args, 'max_tasks', None) or cpus
    else:
        if not getattr(args, 'account', None):
            raise ValueError('Slurm requires --account')
        # Slurm does not use executor.cpus/memory as aggregate limits. Bound
        # submitted jobs conservatively using the largest task request instead.
        jobs = getattr(args, 'max_tasks', None) or 4
        if cpus is not None:
            jobs = min(jobs, cpus // max(t['cpus'] for t in tasks))
        if limit_memory is not None:
            jobs = min(jobs, memory_bytes(limit_memory) // largest_memory)
        options = []
        for key in ('account', 'qos', 'reservation'):
            if getattr(args, key, None):
                options.append('--' + key + '=' + slurm_token(getattr(args, key)))
        if getattr(args, 'partition', None):
            config.append(f'process.queue = {groovy(slurm_token(args.partition))}')
        config += [f'process.clusterOptions = {groovy(" ".join(options))}',
                   f'process.time = {groovy(duration(getattr(args, "time_limit", "24h")))}',
                   f'executor.submitRateLimit = {groovy(str(getattr(args, "submit_rate", 6)) + "/1min")}',
                   "executor.queueStatInterval = '1min'", "executor.pollInterval = '10sec'"]
    config.append(f'executor.queueSize = {jobs}')
    return '\n'.join(config) + '\n', dict(executor=backend, max_cpus=cpus, max_memory=limit_memory, max_tasks=jobs)


def stream_run(command, cwd, env, tasks, console):
    """Mirror controller and worker output to the terminal and a durable log."""
    lines = queue.Queue()
    process = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, errors='replace', bufsize=1)
    def read_controller():
        for line in process.stdout:
            lines.put(line)
    reader = threading.Thread(target=read_controller, daemon=True)
    reader.start()
    offsets = {t['id']: 0 for t in tasks}
    def emit(line):
        print(line, end='', flush=True)
        console.write(line)
        console.flush()
    def drain():
        while True:
            try:
                emit(lines.get_nowait())
            except queue.Empty:
                break
        for t in tasks:
            path = Path(t['log'])
            if path.exists():
                with path.open(errors='replace') as f:
                    f.seek(offsets[t['id']])
                    for line in f:
                        emit(f"[{t['label']}] {line}")
                    offsets[t['id']] = f.tell()
    try:
        while process.poll() is None:
            drain()
            time.sleep(0.2)
        reader.join()
        drain()
    except BaseException:
        process.send_signal(signal.SIGINT)
        try:
            process.wait(timeout=30)
        except subprocess.TimeoutExpired:
            process.terminate()
        raise
    finally:
        process.stdout.close()
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command)


def launch(tasks, output_dir, args, name, dryrun=False):
    if not tasks:
        raise ValueError('No tasks were selected')
    tasks = ordered(tasks)
    output_dir = Path(output_dir).resolve()
    state = output_dir / '.metapathways' / name
    state.mkdir(parents=True, exist_ok=True)
    with (state / 'controller.lock').open('w') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Another {name} command is already using {output_dir}')
        run_id = uuid.uuid4().hex
        run_dir = output_dir / 'logs' / name / run_id
        run_dir.mkdir(parents=True)
        scratch = state / 'tmp' / run_id
        work_override = getattr(args, 'work_dir', None)
        cache_override = getattr(args, 'conda_cache', None)
        work = (Path(work_override).expanduser().resolve() / (name + '-' + run_id) if work_override else scratch / 'work')
        cache = Path(cache_override).expanduser().resolve() if cache_override else scratch / 'conda'
        config_text, resources = configuration(tasks, args, cache)
        for t in tasks:
            key = hashlib.sha256(t['id'].encode()).hexdigest()
            t['receipt'] = str(state / 'receipts' / (key + '.json'))
            t['log'] = str(run_dir / 'tasks' / (key + '.log'))
            t['invocation_receipt'] = str(run_dir / 'tasks' / (key + '.json'))
        manifest = run_dir / 'tasks.json'
        manifest.write_text(json.dumps(tasks, indent=2) + '\n')
        nf = run_dir / 'main.nf'
        for relative, content in render_modules(tasks, manifest).items():
            target = run_dir / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content)
        config = run_dir / 'nextflow.config'
        config.write_text(config_text)
        from metapathways._version import __version__
        summary = dict(mp_version=__version__, resources=resources, work_dir=str(work), conda_cache=str(cache), status='PLANNED')
        summary_path = run_dir / 'summary.json'
        summary_path.write_text(json.dumps(summary, indent=2) + '\n')
        if dryrun:
            for t in tasks:
                print(f"{t['label']}: cpus={t['cpus']} memory={t['memory']} after={','.join(t['dependencies']) or '-'}")
                for cmd in t['commands']:
                    print('  ' + cmd)
            print(f'Execution plan: {manifest}')
            return run_dir
        executable = shutil.which('nextflow')
        if not executable:
            raise RuntimeError('Nextflow is required; activate the MetaPathways environment.')
        if resources['executor'] == 'slurm':
            missing = [c for c in ('sbatch', 'squeue', 'scancel') if not shutil.which(c)]
            if missing:
                raise RuntimeError('Run from your Slurm login node; missing commands: ' + ', '.join(missing))
        session = scratch / 'session'
        session.mkdir(parents=True)
        work.mkdir(parents=True, exist_ok=True)
        cache.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        env.update(NXF_ANSI_LOG='false', PYTHONUNBUFFERED='1',
                   NXF_CONDA_CACHEDIR=str(cache / 'envs'), CONDA_PKGS_DIRS=str(cache / 'packages'))
        # Isolate resource limits from ambient ~/.nextflow/config overrides.
        command = [executable, '-C', str(config), '-log', str(run_dir / 'nextflow.log'),
                   'run', str(nf), '-work-dir', str(work),
                   '-with-trace', str(run_dir / 'trace.tsv'), '-with-report', str(run_dir / 'report.html'),
                   '-with-timeline', str(run_dir / 'timeline.html')]
        cpu_budget = resources['max_cpus'] if resources['max_cpus'] is not None else 'per-job requests'
        print(f"MetaPathways: {name}; executor={resources['executor']}; CPU budget={cpu_budget}; task limit={resources['max_tasks']}", flush=True)
        print(f'Logs and resource reports: {run_dir}', flush=True)
        interrupted = False
        try:
            with (run_dir / 'console.log').open('w') as console:
                stream_run(command, session, env, tasks, console)
            summary['status'] = 'SUCCESS'
        except BaseException as exc:
            interrupted = isinstance(exc, (KeyboardInterrupt, SystemExit))
            summary['status'] = 'INTERRUPTED' if interrupted else 'FAILED'
            summary['error'] = str(exc)
            raise
        finally:
            records = []
            for t in tasks:
                receipt = Path(t['invocation_receipt'])
                if receipt.exists():
                    record = json.loads(receipt.read_text())
                    record['label'] = t['label']
                    records.append(record)
                else:
                    records.append(dict(task=t['id'], label=t['label'], status='NOT_STARTED'))
            for t, record in zip(tasks, records):
                for key in ('sample', 'entity'):
                    if key in t:
                        record[key] = t[key]
            summary['tasks'] = records
            # Preserve wrapper diagnostics even for scheduler/bootstrap failures
            # before a worker could open its own log.
            for path in work.rglob('.command.*'):
                if path.is_file():
                    target = run_dir / 'nextflow_tasks' / path.relative_to(work)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(path, target)
            summary_path.write_text(json.dumps(summary, indent=2) + '\n')
            if getattr(args, 'compact_results', False) and summary['status'] == 'SUCCESS':
                # Keep trace, task plans and diagnostics; discard generated NF code.
                nf.unlink(missing_ok=True)
                config.unlink(missing_ok=True)
                shutil.rmtree(run_dir / 'modules', ignore_errors=True)
                shutil.rmtree(run_dir / 'samples', ignore_errors=True)
                for path in (run_dir / 'nextflow_tasks').rglob('*'):
                    if path.is_file() and path.name not in ('.command.log', '.command.out', '.command.err'):
                        path.unlink()
            if not getattr(args, 'keep_work', False) and not interrupted:
                shutil.rmtree(scratch)
            else:
                print(f'Work/cache retained: {scratch}', flush=True)
        failed = [r['label'] for r in records if r.get('status') == 'FAILED']
        if failed:
            print(f'Completed with {len(failed)} optional PGDB failures: ' + ', '.join(failed))
        print(f'Finished {name}. Logs: {run_dir}', flush=True)
        return run_dir


def annotation_tasks(samples, params, configs, memory_value):
    from metapathways import jobscreator
    creator = jobscreator.JobCreator(params, configs)
    threaded = {'ORF_PREDICTION', 'FUNC_SEARCH', 'SCAN_rRNA', 'SCAN_tRNA', 'COMPUTE_TPM'}
    tasks = []
    for key in sorted(samples):
        s = samples[key]
        creator.addJobs(s, block_mode=True)
        s.writeParamsToRunLogs(type('Parameters', (), {'params': params})())
        previous, group_ids, last_group = [], [], None
        for block, contexts in enumerate(s.getContextBlocks()):
            for i, c in enumerate(contexts):
                group = c.name if c.name == 'SCAN_rRNA:barrnap' else c.name.split(':')[0]
                if group != last_group:
                    previous, group_ids, last_group = group_ids, [], group
                identifier = f'{s.sample_name}:{block}:{i}:{c.name}'
                t = task(identifier, f'{s.sample_name}:{c.name}', c.commands,
                         c.inputs.values(), c.outputs.values(), previous,
                         int(configs['NUM_CPUS']) if c.name.split(':')[0] in threaded else 1,
                         memory_value, kind='annotation', context=vars(c),
                         sample_output=s.output_dir, sample=s.sample_name, block=block)
                t['status'] = c.status
                t['adopt_existing'] = c.name != 'COMPUTE_TPM'
                if c.name.startswith('FUNC_SEARCH:') or c.name == 'COMPUTE_REFSCORES':
                    t['cache_version'] = 'private-fast-scratch-v1'
                    t['adopt_existing'] = False
                if c.name == 'PATHOLOGIC_INPUT':
                    t['cache_version'] = 'normalized-ec-input-v1'
                    t['adopt_existing'] = False
                if c.name == 'COMPUTE_TPM':
                    # bwaFolder is a destination for this task's intermediate
                    # BAM/count files, not a biological input. Fingerprinting
                    # it makes a successful run invalidate its own receipt.
                    t['inputs'] = [value for name, value in c.inputs.items()
                                   if name != 'bwaFolder']
                    if c.temps.get('rev_fq'):
                        t['inputs'].append(c.temps['rev_fq'])
                # Auxiliary inputs affect resumability but may be optional references.
                t['fingerprint_inputs'] = t['inputs'] + list(getattr(c, 'inputs1', {}).values()) + list(getattr(c, 'inputs_optional', {}).values())
                if c.name.startswith('FUNC_SEARCH:'):
                    db = c.name.split(':', 1)[1]
                    t['fingerprint_inputs'].append(str(Path(configs['REFDBS']) / 'functional/formatted' / db) + '.*')
                if c.name.startswith('SCAN_rRNA:') and c.name != 'SCAN_rRNA:barrnap':
                    t['fingerprint_inputs'].append(c.inputs1['dbpath'] + '.*')
                if c.name == 'PATHOLOGIC_INPUT':
                    t['outputs'] += [str(Path(s.output_dir) / 'ptools/0.pf'), str(Path(s.output_dir) / 'ptools/orf_map.txt'),
                                     str(Path(s.output_dir) / 'results/annotation_table' / f'{s.sample_name}.EC_RXN_map.tsv')]
                if c.name == 'CREATE_ANNOT_REPORTS':
                    t['outputs'].append(str(Path(s.output_dir) / 'results/annotation_table' / f'{s.sample_name}.ORF_annotation_table.txt'))
                t['outputs'] = [p for p in t['outputs'] if not Path(p).is_dir()]

                # Read abundance is optional, not a missing-input failure.
                if c.name == 'COMPUTE_TPM' and not c.inputs.get('fwd_fq'):
                    t['status'] = 'skip'
                if c.name == 'COMPUTE_TPM':
                    t['outputs'] += [str(Path(s.output_dir) / 'results/rpkm' / f'{s.sample_name}.contig_counts.tsv'),
                                     str(Path(s.output_dir) / 'results/rpkm' / f'{s.sample_name}.orf_counts.tsv')]
                tasks.append(t)
                group_ids.append(identifier)
    # Some legacy stages intentionally update an earlier stage's output (GBK).
    # The final producer owns its fingerprint; earlier stages still require it.
    owner = {p: t['id'] for t in tasks for p in t['outputs']}
    for t in tasks:
        t['cache_outputs'] = [p for p in t['outputs'] if owner[p] == t['id']]
    return tasks
