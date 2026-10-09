import ast
import glob
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

import pandas as pd

from metapathways.pt_exports import (
    EMPTY_EXPORT_ERROR, PATHWAY_COLUMNS, REPORT_COLUMNS,
    write_verified_empty_pathways,
)


class EmptyPathwayTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.flat = self.root / '1.0/data'
        self.flat.mkdir(parents=True)
        self.reports = self.root / '1.0/reports'
        self.reports.mkdir()
        self.raw = self.flat / 'pathways.dat'
        self.raw.write_text('# Export\n' + EMPTY_EXPORT_ERROR + '\n')
        self.summary = self.reports / 'pathways-report_2026-10-06.txt'
        self.summary.write_text('# Report\n' + ' | '.join(REPORT_COLUMNS) + '\n')
        self.evidence = self.reports / 'pwy-evidence-list.dat'
        self.evidence.write_text(';;; This file contains all inferred pathways and super-pathways\n')
        self.out = self.root / 'sample_pwy.tsv'

    def test_verified_empty_export_preserves_source(self):
        before = self.raw.read_bytes()
        self.assertTrue(write_verified_empty_pathways(self.flat, self.out))
        self.assertEqual(self.raw.read_bytes(), before)
        self.assertEqual(self.out.read_text(), '\t'.join(PATHWAY_COLUMNS) + '\n')

    def test_comment_only_export_also_requires_evidence(self):
        self.raw.write_text('# No instances\n')
        self.assertTrue(write_verified_empty_pathways(self.flat, self.out))

    def test_other_errors_and_real_records_are_not_suppressed(self):
        for text in ('Error: something else\n', 'UNIQUE-ID - PWY-1\n//\n',
                     EMPTY_EXPORT_ERROR + '\nUNIQUE-ID - PWY-1\n//\n'):
            self.raw.write_text(text)
            self.assertFalse(write_verified_empty_pathways(self.flat, self.out))
            self.assertFalse(self.out.exists())

    def test_missing_or_conflicting_reports_fail(self):
        original = self.summary.read_text()
        for text in ('', original + 'pathway | PWY-1\n', 'incorrect header\n'):
            self.summary.write_text(text)
            with self.assertRaises(ValueError):
                write_verified_empty_pathways(self.flat, self.out)
            self.assertFalse(self.out.exists())
        self.summary.unlink()
        with self.assertRaises(ValueError):
            write_verified_empty_pathways(self.flat, self.out)

    def test_missing_or_nonempty_evidence_fails(self):
        for text in ('', ';;; This file contains all inferred pathways and super-pathways\n(PWY-1 RXN-1)\n'):
            self.evidence.write_text(text)
            with self.assertRaises(ValueError):
                write_verified_empty_pathways(self.flat, self.out)
        self.evidence.unlink()
        with self.assertRaises(FileNotFoundError):
            write_verified_empty_pathways(self.flat, self.out)
        self.assertFalse(self.out.exists())

    def test_both_entrypoints_emit_empty_tables_without_calling_camelot(self):
        # Load functions only: these legacy scripts parse CLI arguments at import.
        repo = Path(__file__).resolve().parents[1]
        (self.root / 'samplecyc.tar.bz2').touch()
        annotations = self.root / 'results/annotation_table'
        annotations.mkdir(parents=True)
        pd.DataFrame([dict(orf_id='orf1', EC='1.2.3.4', RXN='RXN-1',
                           **{'ref dbname': 'db', 'target': 'protein', 'product': 'enzyme',
                              'value': '1', 'trim_target': 'protein'})]).to_csv(
            annotations / 'sample.EC_RXN_map.tsv', sep='\t', index=False)
        for name in ('pgdb_build_wf.py', 'pgdb_build_single.py'):
            with self.subTest(name=name):
                tree = ast.parse((repo / 'dev' / name).read_text())
                functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                             and node.name in ('extract_pwy', 'map_orfs2pwys')]
                camelot = Mock(side_effect=AssertionError('Camelot should not parse empty exports'))
                namespace = dict(os=os, glob=glob, pd=pd, make_camelot_file=camelot)
                exec(compile(ast.Module(body=functions, type_ignores=[]), name, 'exec'), namespace)
                namespace['extract_pwy'](str(self.root))
                namespace['map_orfs2pwys'](str(self.root), str(self.root))
                camelot.assert_not_called()
                self.assertTrue(pd.read_csv(self.out, sep='\t').empty)
                mapping = pd.read_csv(self.root / 'sample_pwy2orf.tsv', sep='\t')
                self.assertTrue(mapping.empty)
                self.assertIn('orf_id', mapping.columns)
                self.assertIn('PWY_NAME', mapping.columns)


if __name__ == '__main__':
    unittest.main()
