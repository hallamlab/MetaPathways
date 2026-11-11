"""Unprivileged readiness checks and cached public Pathway Tools dependencies."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def checksum(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def environment():
    env = os.environ.copy()
    for key in ('APPTAINER_BIND', 'APPTAINER_BINDPATH', 'SINGULARITY_BIND', 'SINGULARITY_BINDPATH'):
        env.pop(key, None)
    return env


def host_check(executable, output):
    """Build and execute a tiny local shell container, without network or sudo."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    log = output/'apptainer-check.log'
    log.write_text(f'Apptainer: {executable}\n')
    try:
        with tempfile.TemporaryDirectory(prefix='.pt-check-', dir=output) as directory:
            work = Path(directory)
            root = work/'rootfs'
            for name in ('bin', 'etc', 'tmp', 'var/tmp', 'proc', 'sys', 'dev', 'root'):
                (root/name).mkdir(parents=True, exist_ok=True)
            shell = Path('/bin/sh')
            result = subprocess.run(['ldd', str(shell)], check=True, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
            libraries = {Path(p) for p in re.findall(r'(/[^\s()]+)', result.stdout)}
            for source in [shell, *sorted(libraries)]:
                destination = root/str(source).lstrip('/')
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
                destination.chmod(source.stat().st_mode & 0o777)
            (root/'etc/passwd').write_text('root:x:0:0:root:/root:/bin/sh\n')
            (root/'etc/group').write_text('root:x:0:\n')
            # Relative names avoid quoting user-selected paths inside the definition.
            (work/'probe.def').write_text('Bootstrap: localimage\nFrom: rootfs\n\n%post\n    echo MP-CHECK-READY > /mp-check-ready\n')
            commands = [
                [executable, 'build', '--fakeroot', '--mksquashfs-args', '-processors 1', 'probe.sif', 'probe.def'],
                [executable, 'exec', '--containall', '--cleanenv', '--workdir', str(work/'scratch'),
                 'probe.sif', '/bin/sh', '-ec', 'read value < /mp-check-ready; test "$value" = MP-CHECK-READY'],
            ]
            (work/'scratch').mkdir()
            with log.open('w') as stream:
                stream.write(f'Apptainer: {executable}\n')
                stream.flush()
                for command in commands:
                    subprocess.run(command, cwd=work, env=environment(), check=True,
                                   stdout=stream, stderr=subprocess.STDOUT, timeout=120)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as error:
        detail = log.read_text(errors='replace') if log.exists() else str(error)
        namespace = any(word in detail.lower() for word in ('namespace', 'operation not permitted', 'permission denied', 'newuidmap', 'newgidmap', 'subuid', 'subgid'))
        advice = ('The host denied container setup. Ask the administrator to allow this Mamba Apptainer installation '
                  'to use unprivileged user/mount namespaces and check fakeroot UID mappings/helpers if reported in the log. '
                  'MP cannot grant host permissions and has not changed host security settings.' if namespace else
                  'The local build/execution probe failed. Check the saved log for missing helpers, storage or runtime errors.')
        with log.open('a') as stream:
            stream.write(f'\nProbe error: {error}\n')
        raise RuntimeError(f'Pathway Tools setup is not ready. {advice}\n'
                           f'Apptainer: {executable}\nDiagnostic log: {log}\n'
                           'No licensed patches were downloaded. Ordinary annotation remains available.\n'
                           'Setup guide: https://hallamlab-metapathways.readthedocs.io/en/latest/pathway-tools.html#host-prerequisites') from error
    print(f'Apptainer readiness check passed. Log: {log}', flush=True)


def prepared_base(executable, source=None):
    """Pull a digest-pinned public SIF once; verify cached bytes before reuse."""
    if source:
        path = Path(source).expanduser().resolve()
        if not path.is_file():
            raise ValueError(f'Dependency base must be an existing public dependency SIF: {path}')
        return path
    config = json.loads(Path(__file__).with_name('resources').joinpath('ptools_base.json').read_text())
    reference = config['image']
    if not re.fullmatch(r'docker://quay\.io/hallamlab/metapathways@sha256:[a-f0-9]{64}', reference):
        raise ValueError('Installed public dependency image reference must contain a verified OCI digest')
    root = Path(os.environ.get('XDG_CACHE_HOME', Path.home()/'.cache'))/'metapathways/ptools-base'
    root.mkdir(parents=True, exist_ok=True)
    key = reference.rsplit(':', 1)[1]
    image, receipt = root/(key+'.sif'), root/(key+'.json')
    with (root/(key+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if image.is_file() and receipt.is_file():
            try:
                recorded = json.loads(receipt.read_text())
            except (ValueError, OSError):
                recorded = {}
            if not isinstance(recorded, dict):
                recorded = {}
            if recorded.get('image') == reference and recorded.get('sha256') == checksum(image):
                return image
        with tempfile.TemporaryDirectory(prefix='.download-', dir=root) as directory:
            partial = Path(directory)/'base.sif'
            print('Preparing public container dependencies (cached for later builds).', flush=True)
            subprocess.run([executable, 'pull', str(partial), reference], check=True, env=environment())
            metadata = dict(image=reference, sha256=checksum(partial))
            partial.replace(image)
            from metapathways.nf_worker import atomic_json
            atomic_json(receipt, metadata)
        return image


def validate_base(executable, image):
    # This is a public base, not a licensed image; no source data is copied here.
    subprocess.run([executable, 'exec', '--cleanenv', str(image), '/bin/sh', '-ec',
                    'test -s /opt/mp-base-packages.tsv; test ! -e /opt/pathway-tools; '
                    'test -x /usr/bin/xvfb-run; test -x /usr/bin/makeblastdb; test -x /usr/bin/blastp'],
                   check=True, env=environment(), timeout=120)
