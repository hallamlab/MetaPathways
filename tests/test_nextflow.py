import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from metapathways import nextflow as nf


class SchedulingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.task = nf.task('test', 'test', ['echo hello'], outputs=[str(self.root / 'result')], cpus=8)

    def args(self, **kw):
        p = argparse.ArgumentParser()
        nf.add_resources(p)
        a = p.parse_args([])
        for k, v in kw.items():
            setattr(a, k, v)
        return a

    def test_local_uses_available_cpus(self):
        with patch.object(nf, 'local_capacity', return_value=(48, '128 GB')):
            text, limits = nf.configuration([self.task], self.args(), self.root)
        self.assertEqual(limits['max_cpus'], 48)
        self.assertEqual(limits['max_tasks'], 48)
        self.assertIn('executor.cpus = 48', text)

    def test_slurm_bounds_submissions_by_cpu_and_memory(self):
        args = self.args(executor='slurm', account='lab', partition='compute')
        text, limits = nf.configuration([self.task], args, self.root)
        self.assertEqual(limits['max_tasks'], 4)
        self.assertIn('--account=lab', text)
        self.assertIn("executor.submitRateLimit = '6/1min'", text)
        args.max_cpus = 16
        self.assertEqual(nf.configuration([self.task], args, self.root)[1]['max_tasks'], 2)
        args.max_memory = '16 GB'
        self.assertEqual(nf.configuration([self.task], args, self.root)[1]['max_tasks'], 1)

    def test_local_job_ceiling_uses_detected_capacity_without_user_totals(self):
        args = self.args(max_tasks=100)
        with patch.object(nf, 'local_capacity', return_value=(32, '128 GB')):
            config, limits = nf.configuration([self.task], args, self.root)
        self.assertEqual(limits['max_tasks'], 100)
        self.assertEqual(limits['max_cpus'], 32)
        self.assertEqual(limits['max_memory'], '128 GB')
        self.assertIn('executor.queueSize = 100', config)
        self.assertIn('executor.cpus = 32', config)
        self.assertIn("executor.memory = '128 GB'", config)

    def test_slurm_job_limit_needs_no_implicit_aggregate_budgets(self):
        tasks = [nf.task('threaded', 'threaded', [], cpus=8, memory='64 GB'),
                 nf.task('serial', 'serial', [], cpus=1, memory='64 GB')]
        args = self.args(executor='slurm', account='lab', max_tasks=100)
        with patch.object(nf, 'local_capacity', side_effect=AssertionError('Slurm must not use headnode capacity')):
            config, limits = nf.configuration(tasks, args, self.root)
        self.assertEqual(limits['max_tasks'], 100)
        self.assertIsNone(limits['max_cpus'])
        self.assertIsNone(limits['max_memory'])
        self.assertIn('executor.queueSize = 100', config)
        self.assertNotIn('process.queue =', config)
        self.assertIn('--account=lab', config)
        args.max_memory = '128 GB'
        self.assertEqual(nf.configuration(tasks, args, self.root)[1]['max_tasks'], 2)
        args.max_memory, args.max_cpus = None, 24
        self.assertEqual(nf.configuration(tasks, args, self.root)[1]['max_tasks'], 3)
        args.max_cpus, args.max_tasks = None, None
        self.assertEqual(nf.configuration(tasks, args, self.root)[1]['max_tasks'], 4)

    def test_invalid_slurm_options_and_oversized_tasks_rejected(self):
        with self.assertRaisesRegex(ValueError, 'requires --account'):
            nf.configuration([self.task], self.args(executor='slurm'), self.root)
        with self.assertRaises(argparse.ArgumentTypeError):
            nf.configuration([self.task], self.args(executor='slurm', account='lab; echo unsafe', partition='compute'), self.root)
        with self.assertRaisesRegex(ValueError, 'more memory'):
            nf.configuration([self.task], self.args(max_memory='1 GB', max_cpus=8), self.root)
        with self.assertRaisesRegex(ValueError, 'threads'):
            nf.configuration([self.task], self.args(max_cpus=1), self.root)

    def test_dependency_order_and_cycle_check(self):
        a = nf.task('a','a',[])
        b = nf.task('b','b',[],dependencies=['a'])
        self.assertEqual([t['id'] for t in nf.ordered([b,a])], ['a','b'])
        a['dependencies'] = ['b']
        with self.assertRaisesRegex(ValueError, 'Cyclic'):
            nf.ordered([a,b])

    def test_large_workflow_is_split_into_bounded_modules(self):
        tasks = [nf.task(str(i), f'task {i}', ['true']) for i in range(6470)]
        files = nf.render_modules(tasks, self.root / 'tasks.json')
        modules = [text for name, text in files.items() if name != 'main.nf']
        self.assertEqual(len(modules), 102)
        self.assertEqual(sum(text.count('process TASK_') for text in modules), 6470)
        self.assertTrue(all(text.count('process TASK_') <= 64 for text in modules))
        self.assertNotIn('process TASK_', files['main.nf'])

    def test_module_boundaries_preserve_individual_dependencies(self):
        tasks = [nf.task('a', 'a', []), nf.task('b', 'b', []),
                 nf.task('c', 'c', [], dependencies=['a']),
                 nf.task('d', 'd', [], dependencies=['c', 'b'])]
        files = nf.render_modules(tasks, self.root / 'tasks.json', batch_size=2)
        self.assertIn('BATCH_0001(BATCH_0000.out.done_TASK_0000, BATCH_0000.out.done_TASK_0001)', files['main.nf'])
        module = files['modules/BATCH_0001.nf']
        self.assertIn('TASK_0002(upstream_0)', module)
        self.assertIn('TASK_0003(TASK_0002.out, upstream_1)', module)
        self.assertNotIn('collect', module)
        tasks[2]['dependencies'] = ['missing']
        with self.assertRaisesRegex(ValueError, 'Unknown dependency'):
            nf.render_modules(tasks, self.root / 'tasks.json')

    def test_compact_sample_fan_in_is_wired_inside_sample_modules(self):
        tasks = []
        for sample in range(49):
            ids = [f'{sample}:pgdb:{i}' for i in range(132)]
            tasks.extend(nf.task(i, i, ['true'], sample=str(sample)) for i in ids)
            tasks.append(nf.task(f'{sample}:compact', 'compact', ['true'],
                                 sample=str(sample), dependencies=ids))
        files = nf.render_modules(tasks, self.root/'tasks.json')
        self.assertLess(len(files['main.nf']), 10000)
        self.assertNotIn('.out.', files['main.nf'])
        self.assertEqual(files['main.nf'].count('include { SAMPLE_'), 49)
        for i in range(49):
            module = files[f'samples/SAMPLE_{i:04d}/main.nf']
            self.assertIn(f'workflow SAMPLE_{i:04d}', module)
            self.assertIn('BATCH_0000.out.done_TASK_0000', module)
            self.assertIn('BATCH_0002([upstream_0:', module)
            cleanup = files[f'samples/SAMPLE_{i:04d}/modules/BATCH_0002.nf']
            self.assertIn('Channel.empty().mix(', cleanup)
            self.assertIn('.collect()', cleanup)
            self.assertNotIn('val dependency_1', cleanup)
            self.assertIn('    take:\n    upstream\n', cleanup)
        # Never hide a cross-sample dependency behind independent workflows.
        tasks[-1]['dependencies'].append('0:pgdb:0')
        files = nf.render_modules(tasks, self.root/'tasks.json')
        self.assertNotIn('samples/SAMPLE_0000/main.nf', files)

    def test_database_restart_restores_directories_without_invalidating_downloads(self):
        from metapathways.nf_databases import plan
        tasks = plan(self.root, ['swissprot', 'cazy'], 'fast')
        nf.ordered(tasks)
        directories = tasks[0]
        (self.root / '.metapathways').mkdir()
        subprocess.run(['bash', '-ec', directories['commands'][0]], check=True)
        ready = Path(directories['outputs'][0])
        stamp = ready.stat().st_mtime_ns
        (self.root / 'functional/formatted').rmdir()
        subprocess.run(['bash', '-ec', directories['commands'][0]], check=True)
        self.assertTrue((self.root / 'functional/formatted').is_dir())
        self.assertEqual(ready.stat().st_mtime_ns, stamp)
        self.assertEqual(directories['status'], 'redo')
        self.assertTrue({'release_silva', 'release_cazy'} <= {t['id'] for t in tasks})

    def test_legacy_tool_output_is_streamed_and_returned(self):
        from metapathways.sysutil import getstatusoutput
        output = io.StringIO()
        with patch.dict(os.environ, METAPATHWAYS_STREAM_TOOLS='1'), contextlib.redirect_stdout(output):
            status, captured = getstatusoutput("printf 'stdout\\n'; printf 'stderr\\n' >&2; exit 3")
        self.assertNotEqual(status, 0)
        self.assertEqual(captured, 'stdout\nstderr')
        self.assertEqual(output.getvalue(), 'stdout\nstderr\n')

    def test_annotation_searches_share_a_barrier(self):
        from metapathways.context import Context
        from metapathways import jobscreator
        contexts = []
        for name in ('FILTER_AMINOS', 'FUNC_SEARCH:swissprot', 'FUNC_SEARCH:cazy',
                     'PARSE_FUNC_SEARCH:swissprot', 'PARSE_FUNC_SEARCH:cazy', 'ANNOTATE_ORFS'):
            c = Context()
            c.name, c.status = name, 'yes'
            c.outputs = {'result': str(self.root/name)}
            contexts.append(c)
        sample = SimpleNamespace(sample_name='sample', output_dir=str(self.root),
            getContextBlocks=lambda: [contexts], writeParamsToRunLogs=lambda _: None)
        with patch.object(jobscreator, 'JobCreator'):
            tasks = nf.annotation_tasks({'sample': sample}, {},
                {'NUM_CPUS': 8, 'REFDBS': str(self.root)}, '16 GB')
        self.assertEqual([t['cpus'] for t in tasks], [1, 8, 8, 1, 1, 1])
        self.assertEqual(tasks[1]['dependencies'], [tasks[0]['id']])
        self.assertEqual(tasks[2]['dependencies'], tasks[1]['dependencies'])
        self.assertEqual(tasks[3]['dependencies'], [tasks[1]['id'], tasks[2]['id']])
        self.assertEqual(tasks[4]['dependencies'], tasks[3]['dependencies'])
        self.assertEqual(tasks[5]['dependencies'], [tasks[3]['id'], tasks[4]['id']])

    def test_all_thread_capable_stages_reserve_the_tool_budget(self):
        from metapathways.context import Context
        from metapathways import jobscreator
        names = ('ORF_PREDICTION', 'SCAN_rRNA:barrnap', 'SCAN_tRNA',
                 'FUNC_SEARCH:swissprot', 'COMPUTE_TPM', 'ANNOTATE_ORFS')
        for threads in (1, 3, 8):
            contexts = []
            for name in names:
                c = Context()
                c.name, c.status = name, 'yes'
                c.inputs = {'fwd_fq': 'reads.fq'}
                contexts.append(c)
            sample = SimpleNamespace(sample_name='sample', output_dir=str(self.root),
                getContextBlocks=lambda: [contexts], writeParamsToRunLogs=lambda _: None)
            with self.subTest(threads=threads), patch.object(jobscreator, 'JobCreator'):
                tasks = nf.annotation_tasks({'sample': sample}, {},
                    {'NUM_CPUS': threads, 'REFDBS': str(self.root)}, '16 GB')
                self.assertEqual([t['cpus'] for t in tasks], [threads]*5 + [1])
                config, limits = nf.configuration(tasks, self.args(max_cpus=16, max_memory='128 GB'), self.root)
                self.assertEqual(limits['max_tasks'], 16)
                self.assertIn('executor.cpus = 16', config)

    def fake_run(self, command, cwd, env, tasks, console):
        work = Path(command[command.index('-work-dir')+1])
        (work / 'fixture').mkdir()
        (work / 'fixture/.command.log').write_text('scheduler and tool output\n')
        console.write('all terminal output\n')
        self.assertTrue(env['CONDA_PKGS_DIRS'].startswith(str(self.root)))
        for t in tasks:
            Path(t['log']).parent.mkdir(parents=True, exist_ok=True)
            Path(t['log']).write_text('tool output\n')
            Path(t['invocation_receipt']).write_text(json.dumps(dict(status='SUCCESS')))

    def test_tpm_resume_ignores_own_scratch_but_checks_inputs_and_results(self):
        from metapathways.context import Context
        from metapathways import jobscreator, nf_worker
        for layout in ('paired', 'interleaved', 'single'):
            with self.subTest(layout=layout):
                root = self.root / layout
                root.mkdir()
                (root / 'run_statistics').mkdir()
                scratch = root / 'bwa'
                scratch.mkdir()
                c = Context()
                c.name, c.status = 'COMPUTE_TPM', 'yes'
                c.inputs = {key: str(root / name) for key, name in
                            [('fwd_fq', 'R1.fq'), ('output_gff', 'genes.gff'),
                             ('output_fas', 'assembly.fa')]}
                for value in c.inputs.values():
                    Path(value).write_text('original input')
                c.inputs['bwaFolder'] = str(scratch)
                c.temps = {'rev_fq': None, 'inter': layout == 'interleaved'}
                if layout == 'paired':
                    mate = root / 'R2.fq'
                    mate.write_text('original mate')
                    c.temps['rev_fq'] = str(mate)
                c.outputs = {'stats_file': str(root / 'stats.txt')}
                sample = SimpleNamespace(sample_name='sample', output_dir=str(root),
                    getContextBlocks=lambda: [[c]], writeParamsToRunLogs=lambda _: None)
                with patch.object(jobscreator, 'JobCreator'):
                    t, = nf.annotation_tasks({'sample': sample}, {},
                        {'NUM_CPUS': 4, 'REFDBS': str(root)}, '16 GB')
                t['receipt'] = str(root / 'checkpoint.json')

                def run_stage(*_):
                    (scratch / 'sorted.bam').write_text('generated intermediate')
                    for value in t['outputs']:
                        p = Path(value)
                        p.parent.mkdir(parents=True, exist_ok=True)
                        p.write_text('result')
                    return (0, '')

                def execute():
                    self.assertEqual(nf_worker.execute(t), 0)
                    return json.loads(Path('receipt.json').read_text())['status']

                cwd = Path.cwd()
                try:
                    os.chdir(root)
                    with patch('metapathways.execution.execute', side_effect=run_stage) as run, contextlib.redirect_stdout(io.StringIO()):
                        self.assertEqual(execute(), 'SUCCESS')
                        self.assertEqual(execute(), 'ALREADY_COMPUTED')
                        (scratch / 'sorted.bam').write_text('changed scratch')
                        (scratch / 'extra.txt').write_text('new scratch file')
                        self.assertEqual(execute(), 'ALREADY_COMPUTED')
                        self.assertEqual(run.call_count, 1)
                        for value in t['inputs']:
                            with Path(value).open('a') as handle:
                                handle.write('changed true input')
                            self.assertEqual(execute(), 'SUCCESS')
                            self.assertEqual(execute(), 'ALREADY_COMPUTED')
                        for value in t['outputs']:
                            Path(value).unlink()
                            self.assertEqual(execute(), 'SUCCESS')
                finally:
                    os.chdir(cwd)

    @patch.object(nf.shutil, 'which', return_value='/bin/nextflow')
    def test_default_cleanup_preserves_logs_and_results(self, _):
        sentinel = self.root / 'precious.txt'
        sentinel.write_text('keep')
        with patch.object(nf, 'stream_run', side_effect=self.fake_run):
            run = nf.launch([self.task], self.root, self.args(max_cpus=8, max_memory='64 GB'), 'test')
        summary = json.loads((run/'summary.json').read_text())
        self.assertFalse(Path(summary['work_dir']).exists())
        self.assertFalse(Path(summary['conda_cache']).exists())
        self.assertEqual(sentinel.read_text(), 'keep')
        self.assertEqual((run/'console.log').read_text(), 'all terminal output\n')
        self.assertEqual(len(list((run/'nextflow_tasks').rglob('.command.log'))), 1)

    @patch.object(nf.shutil, 'which', return_value='/bin/nextflow')
    def test_explicit_directories_preserved(self, _):
        args = self.args(max_cpus=8, max_memory='64 GB', work_dir=str(self.root/'custom-work'), conda_cache=str(self.root/'custom-cache'))
        with patch.object(nf, 'stream_run', side_effect=self.fake_run):
            run = nf.launch([self.task], self.root, args, 'test')
        summary = json.loads((run/'summary.json').read_text())
        self.assertTrue(Path(summary['work_dir']).exists())
        self.assertTrue(Path(summary['conda_cache']).exists())

    @patch.object(nf.shutil, 'which', return_value='/bin/nextflow')
    def test_failed_run_keeps_diagnostics_before_cleanup(self, _):
        def fail(*args):
            self.fake_run(*args)
            raise subprocess.CalledProcessError(1, ['nextflow'])
        with patch.object(nf, 'stream_run', side_effect=fail):
            with self.assertRaises(subprocess.CalledProcessError):
                nf.launch([self.task], self.root, self.args(max_cpus=8, max_memory='64 GB'), 'test')
        run = next((self.root/'logs/test').iterdir())
        self.assertEqual(json.loads((run/'summary.json').read_text())['status'], 'FAILED')
        self.assertEqual(len(list((run/'nextflow_tasks').rglob('.command.log'))), 1)


if __name__ == '__main__':
    unittest.main()
