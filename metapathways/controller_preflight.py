"""Advisory controller and file-count checks; never change cleanup policy."""
from pathlib import Path
import resource
import sys


def cgroup_limits(proc=Path('/proc'), root=None):
    """Resolve this process's cgroups through mountinfo, including ancestors.

    root is a filesystem prefix used by isolated tests. Namespace-relative
    cgroup paths and non-root mount roots are both supported.
    """
    groups = {}
    try:
        for line in (proc/'self/cgroup').read_text().splitlines():
            _, controllers, path = line.split(':', 2)
            for controller in controllers.split(','):
                groups[controller] = path
        mounts = (proc/'self/mountinfo').read_text().splitlines()
    except (OSError, ValueError):
        return []
    result = []
    for line in mounts:
        before, separator, after = line.partition(' - ')
        fields, details = before.split(), after.split()
        if not separator or len(fields) < 5 or len(details) < 3:
            continue
        if details[0] not in ('cgroup', 'cgroup2'):
            continue
        mount_root, mount = Path(fields[3]), Path(fields[4])
        # Linux mountinfo encodes spaces, tabs, newlines and backslashes.
        import re
        decode = lambda p: Path(re.sub(r'\\([0-7]{3})', lambda m: chr(int(m[1], 8)), str(p)))
        mount_root, mount = decode(mount_root), decode(mount)
        controllers = [''] if details[0] == 'cgroup2' else details[2].split(',')
        for controller in controllers:
            if controller not in groups:
                continue
            group = Path(groups[controller])
            try:
                relative = group.relative_to(mount_root)
            except ValueError:
                if mount_root == Path('/'):
                    continue
                # A cgroup namespace can report a path relative to its root.
                relative = group.relative_to('/')
            current = mount/relative
            while current == mount or mount in current.parents:
                for kind in ('pids', 'memory'):
                    if controller not in ('', kind):
                        continue
                    limit_name = kind+'.max' if details[0] == 'cgroup2' or kind == 'pids' else 'memory.limit_in_bytes'
                    usage_name = kind+'.current' if details[0] == 'cgroup2' or kind == 'pids' else 'memory.usage_in_bytes'
                    path = (root/str(current).lstrip('/')) if root else current
                    try:
                        limit = int((path/limit_name).read_text())
                        usage = int((path/usage_name).read_text())
                    except (OSError, ValueError):
                        continue
                    result.append(dict(kind=kind, limit=limit, current=usage, path=str(current)))
                if current == mount:
                    break
                current = current.parent
    return result


def report(tasks, rows=None, compact=False, file_quota=None):
    limits = cgroup_limits()
    soft, _ = resource.getrlimit(resource.RLIMIT_NPROC)
    print(f'Controller preflight: {len(tasks)} tasks; RLIMIT_NPROC={soft if soft != resource.RLIM_INFINITY else "unlimited"}', flush=True)
    for item in limits:
        print(f"Controller cgroup {item['kind']}: {item['current']}/{item['limit']} ({item['path']})", flush=True)
    if any((x['kind'] == 'pids' and x['limit'] - x['current'] < 512) or
           (x['kind'] == 'memory' and x['limit'] - x['current'] < 8*1024**3) for x in limits) or \
            (soft != resource.RLIM_INFINITY and soft < 1024):
        print('WARNING: constrained controller headroom; use a suitably sized Slurm controller job. '
              '--max_tasks does not cap JVM threads or memory. Limits may be shared with other processes.',
              file=sys.stderr, flush=True)
    if rows is not None:
        contigs = sum(row.get('contig_count', 0) for row in rows)
        mags = sum(len(row['entities']) for row in rows)
        estimate = 4*(mags+len(rows)) if compact else contigs + 150*mags
        print(f'PGDB file-count estimate: approximately {estimate:,} '
              f'{"retained compact products" if compact else "unpacked PGDB files"}; '
              'excludes other outputs and concurrent temporary builds. User quota/usage is not auto-detected.', flush=True)
        if (file_quota and estimate >= file_quota*0.8) or (not compact and estimate >= 100000):
            print('WARNING: PGDB outputs may exhaust your file quota. Check remaining quota with your '
                  'cluster tools; consider --compact_results before starting and size worker scratch '
                  'for concurrent community builds.', file=sys.stderr, flush=True)
    return limits
