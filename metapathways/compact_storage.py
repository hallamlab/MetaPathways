"""Task-local scratch and atomic archive publication for compact workflows."""
from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import tarfile
import tempfile
import uuid


def scratch_root(override=None):
    # Resolve on the worker, never on the controller/login node. TMPDIR alone
    # is not evidence of node-local storage on an arbitrary Slurm cluster.
    candidate = override or os.environ.get('SLURM_TMPDIR')
    if not candidate and os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Compact tasks require job-local scratch: set --scratch_dir '
                           'to your cluster node-local directory when SLURM_TMPDIR is unavailable')
    root = Path(candidate or tempfile.gettempdir()).expanduser().resolve()
    if not root.is_dir():
        raise RuntimeError(f'Scratch directory does not exist on this worker: {root}')
    usage = shutil.disk_usage(root)
    stats = os.statvfs(root)
    print(f'Compact scratch: {root}; {usage.free / 2**30:.1f} GiB free; '
          f'{stats.f_favail} available inodes', flush=True)
    if usage.free < 1024**3 or (stats.f_files and stats.f_favail < 1000):
        raise RuntimeError(f'Insufficient scratch headroom at {root} (need at least 1 GiB and 1000 inodes)')
    return root


@contextmanager
def task_scratch(override=None):
    with tempfile.TemporaryDirectory(prefix='mp-task-', dir=scratch_root(override)) as work:
        yield Path(work)


def publish_file(source, destination):
    """Cross-filesystem copy, then atomic rename; never expose partial output."""
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name('.' + destination.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        shutil.copyfile(source, temporary)
        if temporary.stat().st_size != source.stat().st_size:
            raise RuntimeError(f'Incomplete copy of {source}')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def archive_directory(source, destination):
    """Bundle without deleting the source; callers delete only after success."""
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name('.' + destination.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with tarfile.open(temporary, 'w:gz', compresslevel=1, dereference=False) as archive:
            archive.add(source, arcname=source.name)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def compact_pgdb(output, tag, build):
    """Build and extract locally, publish report products, retain failure evidence."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    root = os.environ.get('METAPATHWAYS_COMPACT_SCRATCH')
    if not root:
        raise RuntimeError('Compact PGDB execution requires worker scratch setup')
    with tempfile.TemporaryDirectory(prefix='pgdb-', dir=root) as temporary:
        working = Path(temporary)/'result'
        working.mkdir()
        try:
            build(str(working))
            names = [tag + suffix for suffix in ('cyc.tar.bz2', '_pwy.tsv', '_pwy2orf.tsv')]
            for name in names:
                if not (working/name).is_file() or not (working/name).stat().st_size:
                    raise RuntimeError(f'Missing PGDB result: {working/name}')
            if (working/'diagnostics').exists():
                archive_directory(working/'diagnostics', output/'diagnostics.tar.gz')
            for name in names:
                publish_file(working/name, output/name)
        except BaseException:
            # Includes interrupted attempts when Python can still execute cleanup.
            # A hard kill/node loss can only preserve already published records.
            archive_directory(Path(temporary), output/f'failed-attempt-{uuid.uuid4().hex}.tar.gz')
            raise
