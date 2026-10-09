"""Licensed-reference preparation using synthetic MetaCyc records."""
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from metapathways import metacyc_db as mc
from metapathways.nf_databases import plan


class MetaCycTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'29.5/data'
        self.source.mkdir(parents=True)
        (self.source/'protseq.fsa').write_text('>gnl|META|M1 protein\nMKWVTFISLLFLFSSAYSRG\n')
        self.dat('proteins', [('M1', [('COMPONENT-OF', 'C1'), ('COMPONENT-OF', 'C2')]),
                             ('C1', [('CATALYZES', 'E1')]), ('C2', [('CATALYZES', 'E2')])])
        self.dat('enzrxns', [('E1', [('REACTION', 'R1')]), ('E2', [('REACTION', 'R2')])])
        self.dat('reactions', [('R1', []), ('R2', [])])
        self.dat('compounds', [('A', []), ('B', [])])
        self.dat('classes', [('Pathways', [('TYPES', 'Generalized-Reactions')]),
                            ('Generalized-Reactions', [('TYPES', 'FRAMES')]),
                            ('Class1', [('TYPES', 'Pathways')])])
        self.dat('pathways', [(p, [('TYPES', 'Class1'), ('REACTION-LAYOUT',
                  '(R1 (:LEFT-PRIMARIES A) (:DIRECTION :L2R) (:RIGHT-PRIMARIES B))')]) for p in ['P1', 'P2']])

    def dat(self, name, records):
        (self.source/(name+'.dat')).write_text(''.join(
            'UNIQUE-ID - '+identifier+'\nCOMMON-NAME - '+identifier+' name\n'
            + ''.join(k+' - '+v+'\n' for k, v in values) + '//\n' for identifier, values in records))

    def test_existing_scripts_keep_reactions_from_all_complexes(self):
        out = self.root/'tables'
        mc.make_tables(self.source, out)
        with (out/mc.TABLES[0]).open() as stream:
            pairs = {(r['MC'], r['RXN']) for r in csv.DictReader(stream, delimiter='\t')}
        self.assertEqual(pairs, {('M1', 'R1'), ('M1', 'R2')})
        with (out/mc.TABLES[2]).open() as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['MetaCyc_Ontology_IDs'], 'Class1|P1')
        self.assertEqual(rows[1]['MetaCyc_Ontology_IDs'], 'Class1|P2')

    def test_incomplete_fasta_only_source_is_rejected(self):
        (self.source/'proteins.dat').unlink()
        with self.assertRaisesRegex(ValueError, 'protseq.fsa alone is insufficient'):
            mc.source_path(self.source/'protseq.fsa')

    def test_screen_concurrency_and_aggregate_reservations(self):
        import argparse
        from metapathways.nf_databases import screen_task
        args = argparse.Namespace(executor='local', max_tasks=32, max_cpus=None, max_memory=None)
        with patch('metapathways.nextflow.local_capacity', return_value=(32, '256 GB')):
            task = screen_task(self.root/'db', self.root/'pt.sif', '8 GB', resources=args)
        self.assertIn('--max_tasks 32', task['commands'][0])
        self.assertEqual(task['cpus'], 32)
        from metapathways.nextflow import memory_bytes
        self.assertEqual(memory_bytes(task['memory']), memory_bytes('256 GB'))
        args.max_cpus, args.max_memory = 12, '32 GB'
        with patch('metapathways.nextflow.local_capacity', return_value=(32, '256 GB')):
            task = screen_task(self.root/'db', self.root/'pt.sif', '8 GB', resources=args)
        self.assertEqual(task['cpus'], 4)
        self.assertIn('--max_tasks 4', task['commands'][0])
        self.assertEqual(memory_bytes(task['memory']), memory_bytes('32 GB'))

    def test_slurm_screen_reserves_selected_concurrency(self):
        import argparse
        from metapathways.nf_databases import screen_task
        args = argparse.Namespace(executor='slurm', max_tasks=8, max_cpus=None, max_memory=None)
        task = screen_task(self.root/'db', self.root/'pt.sif', '4 GB', resources=args)
        self.assertEqual(task['cpus'], 8)
        from metapathways.nextflow import memory_bytes
        self.assertEqual(memory_bytes(task['memory']), memory_bytes('32 GB'))

    @patch('metapathways.nextflow.local_capacity', return_value=(4, '32 GB'))
    def test_metacyc_screen_default_and_opt_out(self, _):
        image = self.root/'pt.sif'
        image.write_bytes(b'licensed-image-fixture')
        tasks = plan(self.root/'db', ['metacyc'], 'fast', metacyc_source=self.source, screen_image=image)
        self.assertEqual([t['id'] for t in tasks], ['directories', 'prepare_metacyc', 'screen_metacyc'])
        self.assertIn('--publish', tasks[-1]['commands'][0])
        tasks = plan(self.root/'db', ['metacyc'], 'fast', metacyc_source=self.source, skip_pt_screen=True)
        self.assertEqual(len(tasks), 2)

    def test_screen_rejects_insufficient_memory(self):
        from metapathways.nf_databases import screen_resources
        with patch('metapathways.nextflow.local_capacity', return_value=(2, '7 GB')):
            with self.assertRaisesRegex(ValueError, 'available screening budget'):
                screen_resources(None, '16 GB')

    def test_metacyc_only_plan_does_not_download_unrelated_references(self):
        tasks = plan(self.root/'db', ['metacyc'], 'fast', metacyc_source=self.source, skip_pt_screen=True)
        self.assertEqual([t['id'] for t in tasks], ['directories', 'prepare_metacyc'])
        self.assertFalse(tasks[-1]['adopt_existing'])
        self.assertTrue(all(any(n in o for o in tasks[-1]['outputs']) for n in mc.TABLES))
        with patch('metapathways.pt_container.registered_image', return_value=None):
            with self.assertRaisesRegex(ValueError, 'licensed source'):
                plan(self.root, ['metacyc'], 'fast')

    def test_failed_index_does_not_replace_existing_reference(self):
        db = self.root/'db'
        (db/'functional').mkdir(parents=True)
        old = db/'functional/metacyc'
        old.write_text('previous reference')
        real_run = subprocess.run
        def run(command, **kwargs):
            if command[0] == 'fastdb':
                raise subprocess.CalledProcessError(1, command)
            return real_run(command, **kwargs)
        with patch.object(mc.subprocess, 'run', side_effect=run):
            with self.assertRaises(subprocess.CalledProcessError):
                mc.prepare(self.source, db, 'fast')
        self.assertEqual(old.read_text(), 'previous reference')
        self.assertFalse(list(db.glob('.metacyc-build-*')))
        self.assertFalse((self.source/mc.TABLES[0]).exists())

    def test_success_publishes_matching_tables_and_provenance(self):
        db = self.root/'db'
        real_run = subprocess.run
        def run(command, **kwargs):
            if command[0] == 'fastdb':
                Path(command[2]+'.prj').write_text('fixture index')
                return subprocess.CompletedProcess(command, 0)
            return real_run(command, **kwargs)
        with patch.object(mc.subprocess, 'run', side_effect=run):
            mc.prepare(self.source, db, 'fast')
        record = json.loads((db/'functional_categories/MetaCyc_provenance.json').read_text())
        self.assertEqual(record['release'], '29.5')
        self.assertEqual(record['source_sha256']['protseq.fsa'], mc.digest(db/'functional/metacyc'))
        self.assertFalse((self.source/mc.TABLES[0]).exists())


if __name__ == '__main__':
    unittest.main()
