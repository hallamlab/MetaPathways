#!/usr/bin/env python3
"""Opt-in licensed integration check: MP SIF -> Pathway Tools SIF/build.

Run on the target host. Nothing is uploaded; installer and SIF remain private.
"""
import argparse
import json
import os
import signal
from pathlib import Path
import shutil
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mp-image', required=True, type=Path)
    p.add_argument('--pt-image', required=True, type=Path)
    p.add_argument('--installer', type=Path, help='also test full nested build_pt')
    p.add_argument('--output', required=True, type=Path, help='new diagnostic directory')
    p.add_argument('--apptainer', default=shutil.which('apptainer'))
    p.add_argument('--timeout', type=int, default=900, help='seconds allowed per check [900]')
    a = p.parse_args()
    if a.timeout <= 0:
        p.error('--timeout must be positive')
    if not a.apptainer:
        p.error('Host Apptainer is required')
    images = [a.mp_image.resolve(), a.pt_image.resolve()]
    if a.installer:
        images.append(a.installer.resolve())
    for path in images:
        if not path.is_file() or any(c in str(path) for c in ':,\n'):
            p.error(f'Invalid or missing input: {path}')
    output = a.output.resolve()
    if any(c in str(output) for c in ':,\n'):
        p.error('Output path cannot contain colons, commas or newlines')
    output.mkdir(parents=True, exist_ok=False)
    command = [a.apptainer, 'exec', '--cleanenv', '--bind', f'{output}:{output}',
               '--pwd', str(output)]
    for path in images[1:]:
        command += ['--bind', f'{path}:{path}:ro']
    for name, value in [('XDG_DATA_HOME', output / 'config'),
                        ('NXF_HOME', output / 'nextflow'),
                        ('APPTAINER_CACHEDIR', output / 'cache'),
                        ('APPTAINER_TMPDIR', output / 'tmp')]:
        Path(value).mkdir(exist_ok=True)
        command += ['--env', f'{name}={value}']
    command += [str(images[0])]
    tests = [('runtime', ['python', '-c',
             'import sys; from metapathways.pt_container import validate; '
             'print(validate(sys.argv[1], sys.argv[2]))',
             str(images[1]), str(output / 'runtime')])]
    if a.installer:
        tests.append(('build', ['metapathways', 'build_pt', '-i', str(images[2]),
                     '-o', str(output / 'built'), '--threads', '2',
                     '--max_cpus', '2', '--memory', '4 GB', '--max_memory', '4 GB']))
    results = {'mp_image': str(images[0]), 'pt_image': str(images[1]), 'checks': {}}
    for name, args in tests:
        print(f'Testing nested {name}; log: {output / (name + ".log")}', flush=True)
        with (output / (name + '.log')).open('w') as log:
            process = subprocess.Popen(command + args, stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
            timed_out = False
            try:
                code = process.wait(timeout=a.timeout)
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                if isinstance(error, KeyboardInterrupt):
                    raise
                code, timed_out = 124, True
        results['checks'][name] = {'exit_code': code, 'passed': code == 0,
                                  'timed_out': timed_out}
        (output / 'result.json').write_text(json.dumps(results, indent=2) + '\n')
        print(f'{name}: {"PASS" if code == 0 else "FAIL"}', flush=True)
    return int(any(not r['passed'] for r in results['checks'].values()))


if __name__ == '__main__':
    sys.exit(main())
