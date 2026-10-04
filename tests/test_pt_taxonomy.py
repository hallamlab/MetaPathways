import argparse
import contextlib
import io
import unittest
from metapathways.pt_taxonomy import add_taxonomy_options, resolve_taxon
from metapathways import pipeline, analysis_workflow


class TaxonomyTests(unittest.TestCase):
    def parser(self):
        parser = argparse.ArgumentParser()
        add_taxonomy_options(parser)
        return parser

    def test_names_and_alias_resolve(self):
        for name, taxon in [('all', 131567), ('bacteria', 2), ('archaea', 2157),
                            ('eukaryotes', 2759), ('euks', 2759)]:
            with self.subTest(name=name):
                self.assertEqual(resolve_taxon(self.parser().parse_args(['--taxonomic_scope', name])), taxon)

    def test_omission_preserves_existing_behavior_and_numeric_ids_work(self):
        self.assertIsNone(resolve_taxon(self.parser().parse_args([])))
        self.assertEqual(resolve_taxon(self.parser().parse_args(['--taxon_id', '562'])), 562)

    def test_invalid_and_conflicting_options_fail_early(self):
        for argv in [['--taxon_id', '0'], ['--taxon_id', '-1'],
                     ['--taxonomic_scope', 'prokaryotes'], ['--taxonomic_scope', 'unknown'],
                     ['--taxon_id', '2', '--taxonomic_scope', 'all']]:
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.parser().parse_args(argv)

    def test_both_commands_expose_option(self):
        args = pipeline.ptParser().parse_args(['ptools', '-o', '/tmp/example', '--taxonomic_scope', 'all'])
        self.assertEqual(resolve_taxon(args), 131567)
        args = analysis_workflow.parser().parse_args(['analysis_wf', '-i', '/tmp/in', '-o', '/tmp/out',
                                                      '-d', '/tmp/db', '--taxonomic_scope', 'euks'])
        self.assertEqual(resolve_taxon(args), 2759)
