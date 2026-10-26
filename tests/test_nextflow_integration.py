"""Opt-in real-runtime checks: MP_TEST_NEXTFLOW=/path/to/nextflow python -m unittest discover -s tests -p test_nextflow_integration.py."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from metapathways import nextflow as nf


@unittest.skipUnless(os.environ.get('MP_TEST_NEXTFLOW'), 'set MP_TEST_NEXTFLOW for real Nextflow checks')
class RuntimeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        parser = argparse.ArgumentParser()
        nf.add_resources(parser)
        self.args = parser.parse_args(['--max_cpus', '4', '--max_memory', '4 GB', '--memory', '256 MB'])
        self.runtime = os.environ['MP_TEST_NEXTFLOW']

    def launch(self, tasks):
        with patch.object(nf.shutil, 'which', return_value=self.runtime):
            return nf.launch(tasks, self.root, self.args, 'runtime')

    def job(self, name, dependencies=(), fail=False):
        target = self.root/name
        command = 'false' if fail else 'echo done > '+shlex.quote(str(target))
        return nf.task(name, name, [command], outputs=[str(target)], dependencies=dependencies,
                       memory='256 MB', adopt_existing=False)

    def test_fan_out_fan_in_optional_failure_and_resume(self):
        source = self.job('source')
        children = [self.job(f'mag{i}', ['source'], fail=i == 0) for i in range(24)]
        children[0]['allow_failure'] = True
        join = self.job('compact', [t['id'] for t in children] + ['source'])
        run = self.launch([source, *children, join])
        records = json.loads((run/'summary.json').read_text())['tasks']
        self.assertEqual(records[-1]['status'], 'SUCCESS')
        self.assertEqual(records[1]['status'], 'FAILED')
        self.assertTrue((self.root/'compact').is_file())
        run = self.launch([source, *children, join])
        records = json.loads((run/'summary.json').read_text())['tasks']
        self.assertEqual(records[0]['status'], 'ALREADY_COMPUTED')
        self.assertEqual(records[-1]['status'], 'ALREADY_COMPUTED')

    def test_fast_branch_can_finish_while_unrelated_root_is_running(self):
        fast = self.job('fast')
        child = self.job('child', ['fast'])
        slow = self.job('slow')
        # Slow cannot finish until the independent fast branch's child runs.
        # A level-wide barrier would deadlock and fail after the bounded wait.
        marker = shlex.quote(str(self.root/'child'))
        slow['commands'] = [f'for i in {{1..300}}; do test -f {marker} && break; sleep 0.1; done; '
                            f'test -f {marker}; echo done > '+shlex.quote(str(self.root/'slow'))]
        run = self.launch([fast, slow, child])
        records = json.loads((run/'summary.json').read_text())['tasks']
        self.assertTrue(all(t['status']=='SUCCESS' for t in records))

    def test_required_failure_blocks_dependents_and_preserves_work(self):
        bad, dependent = self.job('required', fail=True), self.job('dependent', ['required'])
        with self.assertRaises(subprocess.CalledProcessError):
            self.launch([bad, dependent])
        run = next((self.root/'logs/runtime').iterdir())
        summary = json.loads((run/'summary.json').read_text())
        self.assertEqual(summary['status'], 'FAILED')
        self.assertTrue(Path(summary['work_dir']).exists())
        self.assertFalse((self.root/'dependent').exists())

    def test_large_plan_compiles_with_constant_process_count(self):
        tasks = [self.job(f'mag{i}') for i in range(6470)]
        self.root.joinpath('tasks.json').write_text(json.dumps(tasks))
        files = nf.render_modules(tasks, self.root/'tasks.json')
        for name, content in files.items():
            path = self.root/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
        config = self.root/'nextflow.config'
        config.write_text('executor.cpus = 4\n')
        result = subprocess.run([self.runtime, '-C', str(config), '-log', str(self.root/'nextflow.log'),
                                 'run', str(self.root/'main.nf'), '-preview'], cwd=self.root,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertEqual(files['main.nf'].count('include { RUN_TASK as'), 1)
