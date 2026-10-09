import argparse
import json
from pathlib import Path
import tempfile
import unittest

from metapathways.pt_screen import Screen, write_inputs


class ScreenTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        mapping = self.root/'db/functional_categories'
        mapping.mkdir(parents=True)
        (mapping/'MetaCyc-monomer-rxn-pairs.tsv').write_text('MC\tRXN\nM1\tGOOD\nM2\tBAD\n')
        image = self.root/'image.sif'
        image.write_bytes(b'image')
        self.args = argparse.Namespace(output_dir=str(self.root/'out'), image=str(image),
            refdb_dir=str(self.root/'db'), reactions=None, batch_size=100,
            confirm_runs=2, timeout=30, scratch_dir=None, max_tasks=1)
        self.calls = []

    def runner(self, image, directory, reactions, timeout, scratch):
        self.calls.append(tuple(reactions))
        return {'status': 'FAIL' if 'BAD' in reactions else 'PASS'}

    def test_bisection_confirmation_and_resume(self):
        screen = Screen(self.args, self.runner)
        screen.run()
        self.assertEqual(set(screen.candidates), {'BAD'})
        self.assertEqual(self.calls.count(('BAD',)), 2)
        self.assertFalse(screen.interactions)
        self.calls.clear()
        resumed = Screen(self.args, self.runner)
        resumed.run()
        self.assertEqual(set(resumed.candidates), {'BAD'})
        self.assertEqual(self.calls, [])

    def test_batch_interaction_is_not_blacklisted(self):
        def runner(*args):
            return {'status': 'FAIL' if len(args[2]) > 1 else 'PASS'}
        screen = Screen(self.args, runner)
        screen.run()
        self.assertFalse(screen.candidates)
        self.assertEqual(len(screen.interactions), 1)

    def test_inconclusive_is_not_blacklisted(self):
        screen = Screen(self.args, lambda *args: {'status': 'INCONCLUSIVE'})
        with self.assertRaises(RuntimeError):
            screen.run()
        self.assertFalse(screen.candidates)

    def test_failed_control_prevents_candidate(self):
        screen = Screen(self.args, self.runner)
        screen.attempt([], 'baseline')
        screen.runner = lambda *args: {'status': 'FAIL'}
        screen.isolate(['BAD'])
        self.assertFalse(screen.candidates)

    def test_confirmation_must_repeat_failure(self):
        screen = Screen(self.args, self.runner)
        calls = 0
        def runner(*args):
            nonlocal calls
            calls += 1
            return {'status': 'FAIL' if calls == 1 else 'PASS'}
        screen.runner = runner
        screen.isolate(['BAD'])
        self.assertFalse(screen.candidates)

    def test_interrupted_attempt_is_preserved_and_retried(self):
        screen = Screen(self.args, self.runner)
        def interrupt(*args):
            (args[1]/'console.log').write_text('partial evidence')
            raise KeyboardInterrupt()
        screen.runner = interrupt
        with self.assertRaises(KeyboardInterrupt):
            screen.attempt(['BAD'], 'screen')
        screen.runner = self.runner
        self.assertEqual(screen.attempt(['BAD'], 'screen')['status'], 'FAIL')
        preserved = list((screen.output/'attempts').glob('*-interrupted-*'))
        self.assertEqual(len(preserved), 1)
        self.assertEqual((preserved[0]/'console.log').read_text(), 'partial evidence')

    def test_inconclusive_receipt_retried_on_resume(self):
        screen = Screen(self.args, lambda *args: {'status': 'INCONCLUSIVE'})
        screen.attempt(['BAD'], 'screen')
        resumed = Screen(self.args, self.runner)
        self.assertEqual(resumed.attempt(['BAD'], 'screen')['status'], 'FAIL')
        self.assertEqual(self.calls, [('BAD',)])

    def test_full_screen_publication(self):
        self.args.publish = True
        screen = Screen(self.args, self.runner)
        screen.run()
        file = Path(self.args.refdb_dir)/'functional_categories/ptools_reaction_compatibility.json'
        data = json.loads(file.read_text())
        self.assertEqual(set(data['reactions']), {'BAD'})
        self.assertEqual(data['mapping_sha256'], screen.identity['mapping_sha256'])

    def test_changed_image_cannot_resume(self):
        Screen(self.args, self.runner)
        Path(self.args.image).write_bytes(b'new image')
        with self.assertRaises(ValueError):
            Screen(self.args, self.runner)

    def test_inputs_have_sequences_and_exact_reaction_ids(self):
        directory = self.root/'input'
        write_inputs(directory, ['BAD', 'GOOD'], 'MPscreen')
        self.assertIn('METACYC\tBAD\n', (directory/'contig.pf').read_text())
        sequence = (directory/'contig.fasta').read_text().splitlines()[1]
        self.assertEqual(len(sequence), 480)
        self.assertIn('ENDBASE\t480', (directory/'contig.pf').read_text())
