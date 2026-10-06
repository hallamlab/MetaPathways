"""Execute one planned task; preserve MP outputs and check resumability."""
import glob
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
import tempfile
from contextlib import ExitStack


def fingerprint(paths):
    result = []
    for value in paths:
        if not value:
            result.append([value, None])
            continue
        matches = sorted(glob.glob(str(value))) or [str(value)]
        for name in matches:
            p = Path(name)
            if not p.exists() and Path(name + '.gz').exists():
                p = Path(name + '.gz')
            if p.exists():
                st = p.stat()
                result.append([str(p), st.st_size, st.st_mtime_ns])
                if p.is_dir():
                    # A directory mtime alone does not detect edits to its files.
                    for child in sorted(p.rglob('*')):
                        if child.is_file():
                            st = child.stat()
                            result.append([str(child), st.st_size, st.st_mtime_ns])
            else:
                result.append([str(p), None])
    return result


def available(paths):
    return bool(paths) and all(value and (Path(value).exists() or Path(str(value)+'.gz').exists()) for value in paths)


def atomic_json(p, obj):
    p = Path(p)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Concurrent nodes must never share a staging filename, even when a
    # filesystem does not provide cross-node advisory locking.
    temp = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8',
                                         prefix=f'.{p.name}.', suffix='.tmp',
                                         dir=p.parent, delete=False) as handle:
            temp = Path(handle.name)
            json.dump(obj, handle, indent=2)
            handle.write('\n')
        temp.replace(p)
    finally:
        if temp is not None:
            temp.unlink(missing_ok=True)


def execute(t):
    print(f"Starting {t['label']} (cpus={t['cpus']}, memory={t['memory']})", flush=True)
    start = time.monotonic()
    receipt_path = Path(t['receipt'])
    previous = json.loads(receipt_path.read_text()) if receipt_path.exists() else None
    signature_data = {k: t.get(k) for k in ['commands', 'context', 'cpus', 'kind']}
    if 'cache_version' in t:
        signature_data['cache_version'] = t['cache_version']
    if signature_data['context']:
        signature_data['context'] = {k: v for k, v in signature_data['context'].items() if k not in ('status', 'message')}
    signature = hashlib.sha256(json.dumps(signature_data, sort_keys=True).encode()).hexdigest()
    input_state = fingerprint(t.get('fingerprint_inputs', t['inputs']))
    output_state = fingerprint(t.get('cache_outputs', t['outputs']))
    status = 'FAILED'
    error = ''
    def log(status, duration):
        if t.get('kind') == 'annotation':
            p = Path(t['sample_output']) / 'metapathways_steps_log.txt'
            with p.open('a') as f:
                f.write(f"{t['context']['name']}\t{status} - Time elapsed: {duration:.2f} seconds\n")
    try:
        if t['status'] == 'skip' or (t.get('skip_if_missing') and not Path(t['skip_if_missing']).is_file()):
            status = 'SKIPPED'
        elif (t['status'] != 'redo' and available(t['outputs']) and
              ((previous and previous['status'] == 'SUCCESS' and previous['signature'] == signature
                and previous['inputs'] == input_state and previous['outputs'] == output_state)
               or (previous is None and t.get('adopt_existing', True)))):
            status = 'ALREADY_COMPUTED'
        else:
            if any(not value or (not glob.glob(str(value)) and not Path(str(value)+'.gz').exists()) for value in t['inputs']):
                raise RuntimeError(f"Missing required input: {t['inputs']}")
            # Invalidate the previous success before touching outputs, even if killed.
            atomic_json(receipt_path, dict(status='RUNNING'))
            for command in t['commands']:
                print('Command: ' + command, flush=True)
            if t.get('kind') == 'annotation':
                from types import SimpleNamespace
                from metapathways.context import Context
                from metapathways.metapathways_utils import WorkflowLogger
                from metapathways.execution import execute as run_stage
                c = Context()
                c.__dict__.update(t['context'])
                c.removeOutput()
                base = Path(t['sample_output'])
                s = SimpleNamespace(errorlogger=WorkflowLogger(str(base/'errors_warnings_log.txt'), open_mode='a'),
                    runstatslogger=WorkflowLogger(str(base/'run_statistics'/f"{t['sample']}.run.stats.txt"), open_mode='a'))
                result = run_stage(s, c)
                if result is None or result[0] != 0:
                    raise RuntimeError(f'Stage failed: {result}')
            else:
                for cmd in t['commands']:
                    subprocess.run(['bash', '-euo', 'pipefail', '-c', cmd], check=True)
            if not available(t['outputs']):
                raise RuntimeError(f"Task did not produce expected outputs: {t['outputs']}")
            status = 'SUCCESS'
    except Exception as exc:
        error = str(exc)
        traceback.print_exc()
    duration = time.monotonic() - start
    log(status, duration)
    record = dict(task=t['id'], status=status, elapsed_seconds=duration, signature=signature,
                  inputs=input_state, outputs=fingerprint(t.get('cache_outputs', t['outputs'])), error=error)
    # Preserve the cache's SUCCESS on a cache hit; record this invocation separately.
    if status != 'SKIPPED':
        saved = dict(record, status='SUCCESS' if status == 'ALREADY_COMPUTED' else status)
        atomic_json(receipt_path, saved)
    atomic_json('receipt.json', record)
    if t.get('invocation_receipt'):
        atomic_json(t['invocation_receipt'], record)
    print(f"{t['label']}: {status} ({duration:.2f}s)", flush=True)
    if status == 'FAILED' and not t.get('allow_failure', False):
        return 1
    if status == 'FAILED':
        print(f"Optional-entity task failed; inspect its diagnostics: {t['label']}: {error}", file=sys.stderr)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    direct = argv and argv[0] == '--execute'
    if direct:
        argv.pop(0)
    manifest, identifier = argv
    tasks = json.loads(Path(manifest).read_text())
    t = next(t for t in tasks if t['id'] == identifier)
    if direct:
        with ExitStack() as contexts:
            if t.get('compact_results') and not t['id'].endswith(':compact_results'):
                from metapathways.compact_storage import task_scratch
                work = contexts.enter_context(task_scratch(t.get('scratch_dir')))
                os.environ['TMPDIR'] = str(work)
                os.environ['METAPATHWAYS_COMPACT_SCRATCH'] = str(work)
                tempfile.tempdir = None
            if t.get('host_serial'):
                lock = contexts.enter_context(open(f'/tmp/metapathways-ptools-{os.getuid()}.lock', 'a'))
                fcntl.flock(lock, fcntl.LOCK_EX)
            return execute(t)
    logfile = Path(t['log'])
    logfile.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["METAPATHWAYS_STREAM_TOOLS"] = "1"
    # Prevent hidden BLAS/OpenMP pools from exceeding the scheduled CPU request.
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = str(t['cpus'])
    with logfile.open('w') as log:
        process = subprocess.Popen([sys.executable, '-u', '-m', 'metapathways.nf_worker',
                                    '--execute', manifest, identifier], env=env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, errors='replace', bufsize=1)
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end='', flush=True)
            return process.wait()
        finally:
            process.stdout.close()


if __name__ == '__main__':
    sys.exit(main())
