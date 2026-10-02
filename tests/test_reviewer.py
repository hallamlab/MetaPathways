import tempfile
from pathlib import Path
import unittest
from metapathways.reviewer import prepare, SEEDS


class ReviewerTests(unittest.TestCase):
    def test_bundled_data_copies_idempotently_without_indexes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = prepare(directory)
            self.assertEqual(len(list((root/'cami-reviewer/inputs').rglob('*.gz'))), 9)
            self.assertEqual(len(list((root/'cami-reviewer/inputs/mag_maps').glob('*.tsv'))), 3)
            self.assertEqual(len(list((root/'MPDB').rglob('*.*'))), 1)
            self.assertTrue(all((root/'MPDB'/name).is_file() for name in SEEDS))
            file = root/'MPDB/functional/swissprot_test'
            mtime = file.stat().st_mtime_ns
            prepare(root)
            self.assertEqual(file.stat().st_mtime_ns, mtime)

    def test_modified_data_is_preserved_before_any_copy(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'cami-reviewer').mkdir()
            file = root/'cami-reviewer/all.tsv'
            file.write_text('my manifest')
            with self.assertRaisesRegex(ValueError, 'differs'):
                prepare(root)
            self.assertEqual(file.read_text(), 'my manifest')
            self.assertFalse((root/'MPDB').exists())

    def test_reference_symlinks_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'external').mkdir()
            (root/'MPDB').symlink_to(root/'external', target_is_directory=True)
            with self.assertRaisesRegex(ValueError, 'symbolic link'):
                prepare(root)
            self.assertEqual(list((root/'external').iterdir()), [])

    def test_build_db_uses_explicit_writable_destination(self):
        from unittest.mock import patch
        from metapathways import pipeline
        with tempfile.TemporaryDirectory() as directory:
            root = prepare(directory)
            destination = str(root/'MPDB')
            with patch('sys.argv', ['metapathways', 'build_db', '--test', '-d', destination]), \
                 patch('metapathways.nextflow.launch') as launch:
                pipeline.build_db()
            tasks, output, args, command = launch.call_args.args
            self.assertEqual(output, destination)
            self.assertEqual(command, 'build_db')
            self.assertTrue(any(destination in x for t in tasks for x in t['inputs']))
            self.assertFalse(any('/regtests/' in x for t in tasks for x in t['outputs']))
