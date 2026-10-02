"""Multi-sample discovery and DAG tests; no biological tools are executed."""
import csv
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from metapathways import analysis_workflow as wf, nextflow as nf, pipeline


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.inputs = self.root / 'inputs'
        for folder in ('assemblies', 'reads', 'mag_maps'):
            (self.inputs / folder).mkdir(parents=True)
        self.write('assemblies/Alpha.fa', '>ctg1\nACGT\n')
        self.write('assemblies/Beta.fna', '>ctg2\nACGT\n')
        self.write('reads/Alpha_R1.fq', '@r/1\nAC\n+\nII\n')
        self.write('reads/Alpha_R2.fq', '@r/2\nGT\n+\nII\n')
        self.write('reads/Beta_interleaved.fastq', '@b/1\nAC\n+\nII\n')
        self.write('mag_maps/Alpha.tsv', 'ctg1\tmb2.1\n')
        self.write('mag_maps/Beta.tsv', 'ctg2\tmb2.1\n')

    def write(self, name, text):
        (self.inputs / name).write_text(text)

    def args(self, *extra):
        return wf.parser().parse_args(['analysis_wf', '-i', str(self.inputs),
            '-o', str(self.root / 'out'), '-d', '/fake/db', *extra])

    def rows(self):
        return wf.validate(wf.discover(self.args()))

    def test_discovery_pairs_reads_and_normalizes_mag_ids(self):
        rows = self.rows()
        self.assertEqual([r['sample_id'] for r in rows], ['Alpha', 'Beta'])
        self.assertEqual([r['read_layout'] for r in rows], ['paired', 'interleaved'])
        self.assertEqual(rows[0]['entities'], ['mb2_1'])
        self.assertNotEqual(rows[0]['reads_1'], rows[1]['reads_1'])

    def test_missing_mate_points_to_documentation(self):
        (self.inputs / 'reads/Alpha_R2.fq').unlink()
        with self.assertRaisesRegex(ValueError, 'README.md#custom-analysis-manifest'):
            self.rows()

    def test_orphan_and_ambiguous_files_fail(self):
        for name, text in [('reads/unknown_single.fq', 'x'),
                           ('reads/Alpha_single.fq', 'x'),
                           ('assemblies/Alpha.fasta', '>ctg1\nAC\n'),
                           ('mag_maps/unknown.tsv', 'c\tm\n')]:
            with self.subTest(name=name):
                self.write(name, text)
                with self.assertRaises(ValueError):
                    self.rows()
                (self.inputs / name).unlink()

    def test_hardlinked_mates_fail(self):
        mate = self.inputs / 'reads/Alpha_R2.fq'
        mate.unlink()
        os.link(self.inputs / 'reads/Alpha_R1.fq', mate)
        with self.assertRaisesRegex(ValueError, 'same file'):
            self.rows()

    def test_bad_maps_fail_before_jobs(self):
        for text in ('ctg_missing\tMAG1\n', 'ctg1\tcommunity\n',
                     'ctg1\tmb2.1\nctg2\tmb2_1\n', 'ctg1\tMAG1\nctg1\tMAG2\n'):
            with self.subTest(text=text):
                self.write('mag_maps/Alpha.tsv', text)
                with self.assertRaises(ValueError):
                    self.rows()

    def test_manifest_relative_paths_and_mixed_optional_branches(self):
        manifest = self.inputs / 'custom.tsv'
        with manifest.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(wf.FIELDS)
            writer.writerow(['DifferentID', 'assemblies/Alpha.fa', 'none', '', '', ''])
            writer.writerow(['Beta', 'assemblies/Beta.fna', 'single', 'reads/Beta_interleaved.fastq', '', 'mag_maps/Beta.tsv'])
        rows = wf.validate(wf.read_manifest(manifest))
        self.assertEqual(rows[1]['read_layout'], 'none')
        self.assertEqual(rows[1]['entities'], [])
        self.assertTrue(Path(rows[0]['assembly']).is_absolute())

    def test_resolved_manifest_guards_resume_and_aliases(self):
        output = self.root / 'out'
        output.mkdir()
        rows = self.rows()
        staged = wf.save_inputs(rows, output)
        alias = Path(staged[0]['assembly'])
        self.assertTrue(alias.is_symlink())
        stamp = alias.lstat().st_mtime_ns
        wf.save_inputs(rows, output)
        self.assertEqual(alias.lstat().st_mtime_ns, stamp)
        with self.assertRaisesRegex(ValueError, 'different inputs'):
            wf.save_inputs(rows[:1], output)

    def test_dependencies_are_per_sample_and_pgdb_failures_optional(self):
        tasks = []
        for row in self.rows():
            annotation = [nf.task(row['sample_id'] + ':path', 'path', [],
                context={'name': 'PATHOLOGIC_INPUT'})]
            tasks += annotation + wf.downstream(row, self.root / 'out', annotation, self.args(), '/fake/image.sif')
        nf.ordered(tasks)
        indexed = {t['id']: t for t in tasks}
        self.assertEqual(indexed['Alpha:pgdb:community']['dependencies'], ['Alpha:path'])
        self.assertEqual(indexed['Alpha:pgdb:mb2_1']['dependencies'], ['Alpha:mag_split'])
        self.assertEqual(indexed['Beta:mag_split']['dependencies'], ['Beta:path'])
        for sample in ('Alpha', 'Beta'):
            split = indexed[f'{sample}:mag_split']
            self.assertEqual(split['cache_version'], 'authoritative-mag-coordinates-v1')
            self.assertFalse(split['adopt_existing'])
            self.assertIn(str(self.root / f'out/{sample}/results/annotation_table/{sample}.ptinput.tsv'), split['inputs'])
        self.assertTrue(indexed['Alpha:pgdb:mb2_1']['allow_failure'])
        self.assertFalse(indexed['Alpha:pgdb:community']['allow_failure'])
        self.assertEqual(indexed['Alpha:pgdb:mb2_1']['cpus'], 1)
        self.assertEqual(indexed['Alpha:pgdb:community']['memory'], '16 GB')
        for task in tasks:
            self.assertTrue(all(d.split(':')[0] == task['id'].split(':')[0] for d in task['dependencies']))

    def test_one_memory_request_applies_to_all_workflow_jobs_unless_overridden(self):
        row = self.rows()[0]
        annotation = [nf.task('Alpha:path', 'path', [], context={'name': 'PATHOLOGIC_INPUT'})]
        for executor in ('local', 'slurm'):
            args = self.args('--memory', '64 GB', '--executor', executor)
            tasks = wf.downstream(row, self.root / 'out', annotation, args, '/fake/image.sif')
            self.assertTrue(all(t['memory'] == '64 GB' for t in tasks))
            self.assertTrue(all(t['cpus'] == 1 for t in tasks))
            args.ptools_memory = '4 GB'
            tasks = wf.downstream(row, self.root / 'out', annotation, args, '/fake/image.sif')
            for t in tasks:
                self.assertEqual(t['memory'], '4 GB' if ':pgdb:' in t['id'] else '64 GB')

    def test_standalone_mag_split_invalidates_legacy_cache_and_tracks_coordinates(self):
        base = self.root / 'sample'
        (base / 'results/annotation_table').mkdir(parents=True)
        (base / 'preprocessed').mkdir()
        (base / 'results/annotation_table/sample.ORF_annotation_table.txt').touch()
        (base / 'preprocessed/sample.mapping.txt').touch()
        with patch('sys.argv', ['metapathways', 'mag_split', '-o', str(base),
                                '-m', str(self.inputs / 'mag_maps/Alpha.tsv')]), \
             patch.object(nf, 'launch') as launch:
            pipeline.mag_split()
        task = launch.call_args.args[0][0]
        self.assertEqual(task['cache_version'], 'authoritative-mag-coordinates-v1')
        self.assertFalse(task['adopt_existing'])
        self.assertIn(str(base / 'results/annotation_table/sample.ptinput.tsv'), task['inputs'])

    def test_one_launch_contains_all_samples_with_distinct_reads(self):
        planned = []
        def prepare(args, parser):
            planned.append(args)
            name = Path(args.input_file).stem
            return [nf.task(name + ':path', name, [], context={'name': 'PATHOLOGIC_INPUT'})], args.output_dir
        with patch.object(pipeline, 'prepare_annotation', side_effect=prepare), \
             patch.object(wf.shutil, 'which', return_value='/mock/tool'), \
             patch.object(nf, 'launch') as launch:
            wf.main(['-i', str(self.inputs), '-o', str(self.root / 'out'), '-d', '/fake/db',
                     '--skip_ptools', '--threads', '8', '--max_cpus', '32', '--dryrun'])
        launch.assert_called_once()
        self.assertEqual(len(launch.call_args.args[0]), 4)
        self.assertNotEqual(planned[0].fwd_fastq, planned[1].fwd_fastq)
        self.assertFalse(planned[0].interleaved)
        self.assertTrue(planned[1].interleaved)
        self.assertTrue(launch.call_args.kwargs['dryrun'])

    def test_invalid_input_never_launches(self):
        (self.inputs / 'reads/Alpha_R2.fq').unlink()
        with patch.object(nf, 'launch') as launch, self.assertRaises(ValueError):
            wf.main(['-i', str(self.inputs), '-o', str(self.root / 'out'), '-d', '/fake/db', '--skip_ptools'])
        launch.assert_not_called()

    def test_compact_cleanup_is_per_sample_terminal_and_resume_skips_complete(self):
        def prepare(args, parser):
            name = Path(args.input_file).stem
            return [nf.task(name+':path', name, [], sample=name, context={'name':'PATHOLOGIC_INPUT'}),
                    nf.task(name+':tpm', name, [], sample=name)], args.output_dir
        command = ['-i', str(self.inputs), '-o', str(self.root/'out'), '-d', '/fake/db',
                   '--skip_ptools', '--compact_results', '--dryrun']
        with patch.object(pipeline, 'prepare_annotation', side_effect=prepare), \
             patch.object(wf.shutil, 'which', return_value='/mock/tool'), \
             patch.object(nf, 'launch') as launch:
            wf.main(command)
            tasks = launch.call_args.args[0]
            for sample in ('Alpha','Beta'):
                cleanup = next(t for t in tasks if t['id'] == sample+':compact_results')
                self.assertEqual(set(cleanup['dependencies']), {sample+':path', sample+':tpm', sample+':mag_split'})
            with patch('metapathways.compact_results.marker_state', side_effect=['complete', None]):
                wf.main(command)
            self.assertTrue(all(t['sample']=='Beta' for t in launch.call_args.args[0]))

    def test_compact_rejects_shared_cache_and_retained_work(self):
        for flags in (['--keep_work'], ['--work_dir','/tmp/work'], ['--conda_cache','/tmp/cache']):
            with self.assertRaisesRegex(ValueError, 'cannot be combined'):
                wf.main(['--compact_results', *flags])

    def test_planners_do_not_share_parameter_state(self):
        from metapathways.jobscreator import ContextCreator
        a = ContextCreator({'g': {'key': 'a'}}, {'NUM_CPUS': 8})
        b = ContextCreator({'g': {'key': 'b'}}, {'NUM_CPUS': 3})
        self.assertEqual(a.params.get('g', 'key'), 'a')
        self.assertEqual(b.params.get('g', 'key'), 'b')
        self.assertEqual((a.configs.NUM_CPUS, b.configs.NUM_CPUS), (8, 3))

    def test_missing_mag_input_is_skipped_without_launching_tool(self):
        from metapathways.nf_worker import execute
        task = nf.task('mag', 'mag', ['must-not-run'], inputs=['absent'],
            outputs=['absent-output'], skip_if_missing=str(self.root / 'absent/0.pf'),
            receipt=str(self.root / 'cache.json'))
        previous = Path.cwd()
        try:
            os.chdir(self.root)
            with patch('metapathways.nf_worker.subprocess.run') as run:
                self.assertEqual(execute(task), 0)
                run.assert_not_called()
            self.assertEqual(json.loads(Path('receipt.json').read_text())['status'], 'SKIPPED')
        finally:
            os.chdir(previous)


if __name__ == '__main__':
    unittest.main()
