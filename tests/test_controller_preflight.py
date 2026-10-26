import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from metapathways import controller_preflight as preflight


class PreflightTests(unittest.TestCase):
    def test_v2_resolves_nested_mount_and_ancestor_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            proc = root/'proc/self'
            proc.mkdir(parents=True)
            (proc/'cgroup').write_text('0::/users/ryan/job\n')
            (proc/'mountinfo').write_text('1 0 0:1 /users /sys/fs/cgroup rw - cgroup2 cgroup rw\n')
            for folder, limit, current in [('ryan/job','max','20'), ('ryan','512','100'), ('','1024','500')]:
                base = root/'sys/fs/cgroup'/folder
                base.mkdir(parents=True, exist_ok=True)
                (base/'pids.max').write_text(limit)
                (base/'pids.current').write_text(current)
            result = preflight.cgroup_limits(root/'proc', root)
            self.assertEqual([(x['limit'],x['current']) for x in result], [(512,100),(1024,500)])

    def test_v1_and_unknown_controllers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            proc = root/'proc/self';proc.mkdir(parents=True)
            (proc/'cgroup').write_text('2:pids:/job\n3:memory:/job\n')
            (proc/'mountinfo').write_text('1 0 0:1 / /cgroup/pids rw - cgroup cgroup rw,pids\n'
                                         '2 0 0:2 / /cgroup/memory rw - cgroup cgroup rw,memory\n')
            for controller, limit_file, usage_file in [('pids','pids.max','pids.current'),
                    ('memory','memory.limit_in_bytes','memory.usage_in_bytes')]:
                base=root/'cgroup'/controller/'job';base.mkdir(parents=True)
                (base/limit_file).write_text('512');(base/usage_file).write_text('100')
            self.assertEqual({x['kind'] for x in preflight.cgroup_limits(root/'proc',root)}, {'pids','memory'})

    def test_estimate_warns_without_changing_compact_policy(self):
        out, err = io.StringIO(), io.StringIO()
        with patch.object(preflight, 'cgroup_limits', return_value=[]), \
             contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            preflight.report([], [dict(contig_count=90000,entities=['MAG1'])], file_quota=100000)
        self.assertIn('90,150',out.getvalue())
        self.assertIn('--compact_results',err.getvalue())
