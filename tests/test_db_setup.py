import argparse
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

from metapathways import db_setup as db
from metapathways.nf_databases import plan
from metapathways.nextflow import ordered


class SetupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/'MPDB'
        self.args = argparse.Namespace(add=False, update=False, test=False, dryrun=False,
                                       func=['swissprot'], aligner='fast')

    def put(self, name):
        path = self.root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('fixture')

    def existing(self):
        for name in db.SUPPORT:
            self.put(name)
        for name in db.SILVA:
            for suffix in ('.ndb', '-names.txt'):
                self.put('taxonomic/formatted/'+name+suffix)
        self.put('functional/swissprot')

    def test_new_empty_and_unrelated_destinations(self):
        db.preflight(self.root, self.args, 'new')
        self.assertFalse(self.root.exists())
        self.root.mkdir()
        db.preflight(self.root, self.args, 'new')
        self.put('unrelated.txt')
        with self.assertRaisesRegex(ValueError, 'nonempty'):
            db.preflight(self.root, self.args, 'new')

    def test_existing_requires_explicit_operation(self):
        self.existing()
        with self.assertRaisesRegex(ValueError, '--update'):
            db.preflight(self.root, self.args, 'new')
        self.args.func = ['cazy']
        with self.assertRaisesRegex(ValueError, '--add'):
            db.preflight(self.root, self.args, 'new')
        self.args.add = True
        db.preflight(self.root, self.args, 'new')
        self.args.func = ['swissprot']
        with self.assertRaisesRegex(ValueError, 'Already present'):
            db.preflight(self.root, self.args, 'new')
        self.args.add, self.args.update = False, True
        db.preflight(self.root, self.args, 'new')
        self.args.func = ['cazy']
        with self.assertRaisesRegex(ValueError, 'Not present'):
            db.preflight(self.root, self.args, 'new')

    def test_incomplete_database_cannot_be_extended(self):
        self.put('functional/swissprot')
        self.args.add = True
        self.args.func = ['cazy']
        with self.assertRaisesRegex(ValueError, 'missing shared support'):
            db.preflight(self.root, self.args, 'new')

    def test_fresh_metacyc_only_and_mixed_have_shared_support(self):
        image = Path(self.temp.name)/'pt.sif'
        image.write_text('fixture')
        for selection in (['metacyc'], ['swissprot', 'metacyc']):
            tasks = plan(self.root, selection, 'fast', metacyc_source=image, skip_pt_screen=True)
            ordered(tasks)
            ids = {t['id'] for t in tasks}
            self.assertTrue({'prepare_metacyc', 'fetch_enzyme', 'taxonomy_tree'} <= ids)
            self.assertIn('index_'+db.SILVA[0], ids)
            self.assertEqual('fetch_swissprot' in ids, 'swissprot' in selection)
        self.assertFalse(self.root.exists())

    def test_add_meta_and_public_do_not_modify_shared_or_other_references(self):
        image = Path(self.temp.name)/'pt.sif'
        image.write_text('fixture')
        for selection in (['metacyc'], ['cazy']):
            tasks = plan(self.root, selection, 'fast', metacyc_source=image,
                         skip_pt_screen=True, include_support=False)
            ordered(tasks)
            outputs = [o for t in tasks if t['id'] != 'directories' for o in t['outputs']]
            self.assertFalse(any('swissprot' in o or '/taxonomic/' in o or 'enzyme.dat' in o for o in outputs))
        with patch('metapathways.pt_container.registered_image', side_effect=AssertionError('license lookup')):
            plan(self.root, ['cazy'], 'fast', include_support=False)

    def test_missing_licensed_source_stops_without_creating_destination(self):
        with patch('metapathways.pt_container.registered_image', return_value=None):
            with self.assertRaisesRegex(ValueError, 'licensed source'):
                plan(self.root, ['swissprot', 'metacyc'], 'fast')
        self.assertFalse(self.root.exists())

    def test_missing_screen_image_is_preflight_error(self):
        image = Path(self.temp.name)/'pt.sif'
        image.write_text('fixture')
        with self.assertRaisesRegex(ValueError, 'existing SIF'):
            plan(self.root, ['metacyc'], 'fast', metacyc_source=image,
                 screen_image=Path(self.temp.name)/'absent.sif')
        self.assertFalse(self.root.exists())

    def test_resume_and_changed_plan(self):
        tasks = plan(self.root, ['swissprot'], 'fast')
        failed = Mock(side_effect=RuntimeError('interrupted'))
        with self.assertRaises(RuntimeError):
            db.execute(self.root, self.args, tasks, 'original command', failed)
        self.put('functional/swissprot')
        token = db.signature(tasks, self.args)
        db.preflight(self.root, self.args, token)
        with self.assertRaisesRegex(ValueError, 'original command'):
            db.preflight(self.root, self.args, 'changed')
        done = Mock()
        db.execute(self.root, self.args, tasks, 'original command', done)
        record = json.loads((self.root/'.metapathways/database-build.json').read_text())
        self.assertEqual(record['status'], 'complete')

    def test_update_restart_reuses_generation_and_does_not_force_successes(self):
        self.existing()
        self.args.update = True
        def tasks():
            return plan(self.root, ['swissprot'], 'blast', include_support=False)
        generations = []
        def launch(items, *args, **kwargs):
            index = next(t for t in items if t['id'] == 'index_swissprot')
            generations.append(index['cache_version'])
            self.assertNotEqual(index['status'], 'redo')
            self.assertFalse(index['adopt_existing'])
            self.assertTrue(index['commands'][0].startswith('rm -f '))
            if len(generations) == 1:
                raise RuntimeError('interrupted')
        with self.assertRaises(RuntimeError):
            db.execute(self.root, self.args, tasks(), 'update', launch)
        db.execute(self.root, self.args, tasks(), 'update', launch)
        self.assertEqual(generations[0], generations[1])

    def test_dryrun_does_not_write_restart_record(self):
        self.args.dryrun = True
        db.execute(self.root, self.args, [], 'dryrun', Mock())
        self.assertFalse(self.root.exists())

    def test_concurrent_builder_is_rejected(self):
        import fcntl
        folder = self.root/'.metapathways'
        folder.mkdir(parents=True)
        with (folder/'database-build.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            launch = Mock()
            with self.assertRaisesRegex(ValueError, 'Another database builder'):
                db.execute(self.root, self.args, [], 'build', launch)
            launch.assert_not_called()

    @patch('metapathways.nextflow.local_capacity', return_value=(4, '32 GB'))
    def test_registered_sif_and_exported_directory_routes(self, _):
        from metapathways.metacyc_db import REQUIRED
        image = Path(self.temp.name)/'pt.sif'
        image.write_text('fixture')
        with patch('metapathways.pt_container.registered_image', return_value=str(image)):
            tasks = plan(self.root, ['metacyc'], 'fast')
            self.assertIn(str(image), next(t for t in tasks if t['id']=='prepare_metacyc')['inputs'])
        source = Path(self.temp.name)/'data'
        source.mkdir()
        for name in REQUIRED:
            (source/name).write_text('fixture')
        with patch('metapathways.pt_container.registered_image', return_value=None):
            with self.assertRaisesRegex(ValueError, 'screening requires'):
                plan(self.root, ['metacyc'], 'fast', metacyc_source=source)
            tasks = plan(self.root, ['metacyc'], 'fast', metacyc_source=source, screen_image=image)
            self.assertIn('screen_metacyc', [t['id'] for t in tasks])
            tasks = plan(self.root, ['metacyc'], 'fast', metacyc_source=source, skip_pt_screen=True)
            self.assertNotIn('screen_metacyc', [t['id'] for t in tasks])
        (source/'proteins.dat').unlink()
        with self.assertRaisesRegex(ValueError, 'Incomplete MetaCyc'):
            plan(self.root, ['metacyc'], 'fast', metacyc_source=source, skip_pt_screen=True)

    def test_cli_preflight_never_launches_on_invalid_source_or_existing_database(self):
        import contextlib
        import io
        from metapathways import pipeline
        self.existing()
        before = {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        for options in (['--func', 'cazy'], ['--func', 'metacyc', '--metacyc_source', '/missing/source']):
            with patch('sys.argv', ['metapathways', 'build_db', '-d', str(self.root), *options]), \
                 patch('metapathways.nextflow.launch') as launch, contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    pipeline.build_db()
                launch.assert_not_called()
        after = {str(p): p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(before, after)

    def test_add_completed_reference_is_not_silently_replaced(self):
        self.existing()
        self.args.add = True
        self.args.func = ['cazy']
        tasks = plan(self.root, ['cazy'], 'fast', include_support=False)
        db.execute(self.root, self.args, tasks, 'add', lambda *a, **kw: self.put('functional/cazy'))
        with self.assertRaisesRegex(ValueError, 'Already present'):
            db.execute(self.root, self.args, tasks, 'add', Mock())

    def test_changed_metacyc_source_is_not_an_implicit_resume(self):
        image = Path(self.temp.name)/'pt.sif'
        image.write_text('old source')
        self.args.func = ['metacyc']
        tasks = plan(self.root, ['metacyc'], 'fast', metacyc_source=image, skip_pt_screen=True)
        before = db.signature(tasks, self.args)
        image.write_text('new licensed source with different contents')
        self.assertNotEqual(before, db.signature(tasks, self.args))

    def test_update_worker_reuses_completed_tasks_on_resume(self):
        from metapathways.nf_worker import execute
        from metapathways.nextflow import task
        import os
        self.existing()
        self.args.update = True
        output = self.root/'counter'
        def tasks():
            return [task('fixture', 'counter', [f'echo run >> {output}'], outputs=[str(output)])]
        attempt = []
        def launch(items, *a, **kw):
            items[0]['receipt'] = str(self.root/'worker.json')
            # The worker writes its per-attempt receipt in its working directory.
            previous = os.getcwd()
            try:
                os.chdir(self.temp.name)
                result = execute(items[0])
                self.assertEqual(result, 0)
            finally:
                os.chdir(previous)
            attempt.append(True)
            if len(attempt) == 1:
                raise RuntimeError('later task failed')
        with self.assertRaises(RuntimeError):
            db.execute(self.root, self.args, tasks(), 'update', launch)
        db.execute(self.root, self.args, tasks(), 'update', launch)
        self.assertEqual(output.read_text(), 'run\n')
