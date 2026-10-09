import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from metapathways.compact_storage import compact_pgdb, publish_file, scratch_root
from metapathways.compact_results import resume_key
from metapathways import nextflow as nf


class StorageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.scratch = self.root/'scratch'
        self.scratch.mkdir()
        self.output = self.root/'shared'
        self.env = patch.dict(os.environ, METAPATHWAYS_COMPACT_SCRATCH=str(self.scratch))
        self.env.start()
        self.addCleanup(self.env.stop)

    def build(self, directory):
        directory = Path(directory)
        self.assertTrue(directory.is_relative_to(self.scratch))
        for suffix in ('cyc.tar.bz2', '_pwy.tsv', '_pwy2orf.tsv'):
            (directory/('test'+suffix)).write_text('result')
        for folder in ('diagnostics', '1.0/data'):
            (directory/folder).mkdir(parents=True)
            for i in range(100):
                (directory/folder/f'{i}.txt').write_text('evidence')

    def test_success_publishes_tables_archive_diagnostics_without_extracted_tree(self):
        compact_pgdb(self.output, 'test', self.build)
        self.assertEqual(len(list(self.output.iterdir())), 4)
        self.assertFalse(list(self.scratch.iterdir()))
        with tarfile.open(self.output/'diagnostics.tar.gz') as archive:
            self.assertEqual(archive.extractfile('diagnostics/0.txt').read(), b'evidence')

    def test_failed_build_archives_recovery_and_does_not_publish_partial_tables(self):
        def fail(directory):
            self.build(directory)
            raise RuntimeError('inference failed')
        with self.assertRaisesRegex(RuntimeError, 'inference failed'):
            compact_pgdb(self.output, 'test', fail)
        files = list(self.output.iterdir())
        self.assertEqual(len(files), 1)
        with tarfile.open(files[0]) as archive:
            self.assertTrue(any(n.endswith('/1.0/data/0.txt') for n in archive.getnames()))
        self.assertFalse(list(self.scratch.iterdir()))

    def test_failed_copy_never_overwrites_valid_destination(self):
        source = self.root/'source'; source.write_text('new')
        dest = self.root/'dest'; dest.write_text('old')
        def fail(src, target):
            Path(target).write_text('partial')
            raise OSError('quota')
        with patch('metapathways.compact_storage.shutil.copyfile', side_effect=fail):
            with self.assertRaises(OSError):
                publish_file(source, dest)
        self.assertEqual(dest.read_text(), 'old')
        self.assertFalse(list(self.root.glob('.*.tmp')))

    def test_slurm_requires_known_scratch_and_accepts_override(self):
        with patch.dict(os.environ, {'SLURM_JOB_ID':'1'}, clear=True):
            with self.assertRaisesRegex(RuntimeError, '--scratch_dir'):
                scratch_root()
            self.assertEqual(scratch_root(str(self.scratch)), self.scratch)
        with patch.dict(os.environ, {'SLURM_JOB_ID':'1', 'SLURM_TMPDIR':str(self.scratch)}, clear=True):
            self.assertEqual(scratch_root(), self.scratch)

    def test_resource_changes_preserve_compact_sample_key(self):
        assembly = self.root/'in.fasta'; assembly.write_text('>a\nACGT\n')
        row = dict(assembly=str(assembly), reads_1='', reads_2='', mag_map='')
        args = SimpleNamespace(skip_ptools=True, image=None, refdb_dir=str(self.root/'db'),
                               threads=8, max_tasks=100, max_cpus=None, max_memory=None,
                               partition='a', scratch_dir=None)
        key = resume_key(args, row)
        args.max_tasks = 1; args.max_cpus = 16; args.max_memory = '64 GB'
        args.partition = 'b'; args.scratch_dir = '/local'
        self.assertEqual(key, resume_key(args, row))
        args.threads = 4
        self.assertNotEqual(key, resume_key(args, row))

    def test_worker_scratch_is_cleaned_and_checkpoint_survives(self):
        product = self.root/'result.txt'
        task = nf.task('sample:stage', 'stage', [f'printf result > {product}'], outputs=[str(product)])
        task.update(receipt=str(self.root/'checkpoint.json'), compact_results=True,
                    scratch_dir=str(self.scratch), invocation_receipt=str(self.root/'invocation.json'))
        manifest = self.root/'tasks.json'; manifest.write_text(json.dumps([task]))
        env = dict(os.environ, PYTHONPATH=str(Path(nf.__file__).resolve().parent.parent))
        for expected in ('SUCCESS', 'ALREADY_COMPUTED'):
            result = subprocess.run([sys.executable, '-m', 'metapathways.nf_worker', '--execute',
                                     str(manifest), task['id']], cwd=self.root, env=env,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(json.loads((self.root/'invocation.json').read_text())['status'], expected)
            self.assertFalse(list(self.scratch.iterdir()))


if __name__ == '__main__':
    unittest.main()
