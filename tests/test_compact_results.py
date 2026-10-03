"""Destructive cleanup is tested only on temporary fixtures."""
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from metapathways.compact_results import compact, marker_state, MARKER
from metapathways.reporting import build_report
from metapathways.analysis_workflow import compact_task, parser
import test_reporting


class CompactTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        fixture = test_reporting.ReportTests()
        fixture.root = self.root
        fixture.fixture('alpha')
        self.base = self.root/'alpha'
        for name in ('bwa/reads.sorted.bam', 'preprocessed/alpha.fasta',
                     'orf_prediction/alpha.faa', 'blast_results/raw.FASTout',
                     'results/pgdb/community/alphacyc.tar.bz2', 'ptools/0.fasta'):
            path = self.base/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('large intermediate')

    def tables(self):
        build_report(self.root)
        with sqlite3.connect(self.root/'reports/results.sqlite') as db:
            self.assertEqual(db.execute('pragma integrity_check').fetchone()[0], 'ok')
            self.assertFalse(db.execute('pragma foreign_key_check').fetchall())
            return {t:db.execute(f'select * from {t} order by 1,2').fetchall() for t in (
                'contigs','orfs','annotations','annotation_taxonomy','entities',
                'contig_mags','entity_orfs','orf_groups','pathways','pathway_orfs',
                'abundance','abundance_explorer')}

    def test_report_rebuild_equivalence_and_external_symlink_safety(self):
        outside = self.root/'external'; outside.mkdir(); (outside/'important').write_text('keep')
        (self.base/'external_link').symlink_to(outside, target_is_directory=True)
        before = self.tables()
        compact(self.base, 'key')
        self.assertEqual(before, self.tables())
        self.assertEqual((outside/'important').read_text(), 'keep')
        self.assertFalse((self.base/'external_link').exists())
        self.assertFalse((self.base/'bwa').exists())
        self.assertTrue((self.base/'results/pgdb/community/alphacyc.tar.bz2').exists())
        self.assertTrue((self.base/'logs/ptools/run1/summary.json').is_file())
        self.assertEqual(marker_state(self.base, 'key'), 'complete')
        compact(self.base, 'key')
        with self.assertRaisesRegex(ValueError, 'new output'):
            marker_state(self.base, 'different')
        with self.assertRaisesRegex(ValueError, 'new output'):
            marker_state(self.base, 'key', force=True)

    def test_validation_failure_does_not_delete_data(self):
        (self.base/'preprocessed/test.mapping.txt').unlink()
        with self.assertRaisesRegex(ValueError, 'Missing required'):
            compact(self.base, 'key')
        self.assertTrue((self.base/'bwa/reads.sorted.bam').is_file())
        self.assertFalse((self.base/MARKER).exists())

    def test_interrupted_deletion_can_resume(self):
        original = Path.unlink
        def interrupt(path, *a, **kw):
            if path.name == 'reads.sorted.bam':
                raise OSError('simulated interruption')
            return original(path, *a, **kw)
        with patch.object(Path, 'unlink', interrupt):
            with self.assertRaisesRegex(OSError, 'simulated'):
                compact(self.base, 'key')
        self.assertEqual(marker_state(self.base, 'key'), 'compacting')
        compact(self.base, 'key')
        self.assertEqual(marker_state(self.base, 'key'), 'complete')
        self.tables()

    def test_cleanup_is_opt_in_and_waits_for_all_dependencies(self):
        args = parser().parse_args(['analysis_wf'])
        self.assertFalse(args.compact_results)
        t = compact_task('alpha', self.root, 'key', ['alpha:tpm','alpha:pgdb:community','alpha:pgdb:MAG1'], '4 GB')
        self.assertEqual(t['dependencies'], ['alpha:tpm','alpha:pgdb:community','alpha:pgdb:MAG1'])
        self.assertFalse(t.get('allow_failure', False))
        self.assertEqual(t['cpus'], 1)


if __name__ == '__main__':
    unittest.main()
