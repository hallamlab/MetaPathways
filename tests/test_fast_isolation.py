import concurrent.futures
from pathlib import Path
import shlex
import subprocess
import tempfile
import threading
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from metapathways import MetaPathways_func_search as search


class FastIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mp fast ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def options(self, name):
        db = self.root/name
        Path(str(db)+'.prj').write_text('numofsequences=1\n')
        return SimpleNamespace(last_db=str(db), last_o=str(self.root/(name+'.out')),
            last_query=str(self.root/'query.faa'), last_executable='fastal', last_f='2',
            num_threads='2', num_hits='10', run_mode='pervol')

    def test_simultaneous_searches_cannot_share_sort_files(self):
        barrier = threading.Barrier(2)
        roots = []
        def fake(command):
            args = shlex.split(command)
            root = Path(args[args.index('-X')+1])
            roots.append(root)
            # Model FAST processes selecting exactly the same internal filename.
            collision = root/'same-second-LASTtemp0'
            collision.write_text(args[-2])
            barrier.wait(timeout=5)
            Path(args[args.index('-o')+1]).write_text(collision.read_text())
            return 0, ''
        options = [self.options('metacyc'), self.options('swissprot')]
        with patch.object(search.sysutils, 'getstatusoutput', side_effect=fake), concurrent.futures.ThreadPoolExecutor(2) as pool:
            results = list(pool.map(search._execute_FAST, options))
        self.assertTrue(all(r[0] == 0 for r in results))
        self.assertEqual(len(set(roots)), 2)
        self.assertTrue(all(not p.exists() for p in roots))
        for opt in options:
            self.assertEqual(Path(opt.last_o).read_text(), opt.last_db)

    def test_failed_first_volume_stops_without_publishing_partial_results(self):
        opt = self.options('multi')
        Path(opt.last_db+'.prj').write_text('volumes=2\n')
        Path(opt.last_o).write_text('previous output')
        with patch.object(search.sysutils, 'getstatusoutput', return_value=(1, 'failure')) as run:
            self.assertEqual(search._execute_FAST(opt), (1, 'failure'))
        self.assertEqual(run.call_count, 1)
        self.assertEqual(Path(opt.last_o).read_text(), 'previous output')
        self.assertFalse(list(self.root.glob('.fast-*')))

    def test_failure_propagates_without_hidden_retry(self):
        with patch.object(search, '_execute_FAST', return_value=(7, 'broken')) as run:
            self.assertNotEqual(search.main(['--algorithm', 'FAST']), 0)
            self.assertEqual(run.call_count, 1)
            self.assertNotEqual(search.MetaPathways_func_search(['--algorithm', 'FAST'])[0], 0)
            self.assertEqual(run.call_count, 2)

    def test_blast_failure_does_not_launch_fast(self):
        with patch.object(search, '_execute_BLAST', return_value=(3, 'blast failed')), patch.object(search, '_execute_FAST') as fast:
            self.assertEqual(search.main(['--algorithm', 'BLAST']), 3)
            fast.assert_not_called()

    def test_staging_failure_preserves_existing_output(self):
        opt = self.options('multi')
        opt.stage_db = True
        Path(opt.last_o).write_text('previous output')
        with patch.dict(search.os.environ, {'METAPATHWAYS_FAST_SCRATCH': str(self.root/'missing')}), patch.object(search.sysutils, 'getstatusoutput') as run:
            self.assertNotEqual(search._execute_FAST(opt)[0], 0)
            run.assert_not_called()
        self.assertEqual(Path(opt.last_o).read_text(), 'previous output')

    def test_module_exit_status(self):
        result = subprocess.run([search.sys.executable, '-m', 'metapathways.MetaPathways_func_search',
                                 '--algorithm', 'FAST', '--last_o', str(self.root/'missing.out'),
                                 '--last_db', str(self.root/'missing'), '--run_mode', 'pervol'],
                                capture_output=True)
        self.assertNotEqual(result.returncode, 0)

    def test_stages_each_volume_and_removes_previous_volume(self):
        opt = self.options('multi')
        Path(opt.last_db+'.prj').write_text('volumes=2\n')
        for v in range(2):
            Path(opt.last_db+str(v)+'.prj').write_text('index metadata')
            Path(opt.last_db+str(v)+'.suf').write_bytes(b'index data')
        # A different database and volume 10 must not be copied with volume 1.
        Path(opt.last_db+'10.suf').write_text('unrelated')
        opt.stage_db = True
        seen = []
        def fake(command):
            args = shlex.split(command)
            db = Path(args[-2])
            self.assertEqual(db.parent.parent.parent, self.root)
            self.assertEqual(sorted(p.name for p in db.parent.iterdir()), [db.name+'.prj', db.name+'.suf'])
            self.assertEqual(Path(str(db)+'.suf').read_bytes(), b'index data')
            seen.append(db)
            Path(args[args.index('-o')+1]).write_text('')
            return 0, ''
        with patch.dict(search.os.environ, {'METAPATHWAYS_FAST_SCRATCH': str(self.root)}), patch.object(search.sysutils, 'getstatusoutput', side_effect=fake):
            self.assertEqual(search._execute_FAST(opt), (0, ''))
        self.assertEqual(len(seen), 2)
        self.assertTrue(Path(opt.last_o).exists())
        self.assertFalse(list(self.root.glob('.fast-*')))

    def test_bundled_fast_parallel_searches_return_only_their_own_targets(self):
        binary = Path(search.__file__).parent/'bin'
        sequence = 'MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE'
        query = self.root/'query.faa'
        query.write_text('>query\n'+sequence+'\n')
        options = []
        for name in ('meta_fixture', 'swiss_fixture'):
            opt = self.options(name)
            fasta = self.root/(name+'.faa')
            fasta.write_text('>'+name+'\n'+sequence+'\n')
            subprocess.run([str(binary/'fastdb'), '-p', opt.last_db, str(fasta)], check=True, capture_output=True)
            opt.last_executable = str(binary/'fastal')
            opt.stage_db = True
            options.append(opt)
        with patch.dict(search.os.environ, {'METAPATHWAYS_FAST_SCRATCH': str(self.root)}), concurrent.futures.ThreadPoolExecutor(2) as pool:
            results = list(pool.map(search._execute_FAST, options))
        self.assertTrue(all(r[0] == 0 for r in results), results)
        for opt in options:
            targets = {line.split('\t')[1] for line in Path(opt.last_o).read_text().splitlines()}
            self.assertEqual(targets, {Path(opt.last_db).name})


if __name__ == '__main__':
    unittest.main()
