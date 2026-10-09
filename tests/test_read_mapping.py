"""Read-layout regressions, exercised without external bioinformatics tools."""
import inspect
import os
import shlex
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from metapathways import MetaPathways_tpm as tpm
from metapathways import jobscreator


class ReadMappingTests(unittest.TestCase):
    def test_distinct_paired_inputs(self):
        self.assertEqual(tpm.coverm_read_arguments('R1.fq', 'R2.fq'),
                         ['-1', 'R1.fq', '-2', 'R2.fq'])

    def test_interleaved_with_missing_reverse(self):
        for reverse in (None, '', 'None'):
            with self.subTest(reverse=reverse):
                self.assertEqual(tpm.coverm_read_arguments('reads.fq', reverse, True),
                                 ['--interleaved', 'reads.fq'])

    def test_single_with_missing_reverse(self):
        for reverse in (None, '', 'None'):
            with self.subTest(reverse=reverse):
                self.assertEqual(tpm.coverm_read_arguments('reads.fq', reverse),
                                 ['--single', 'reads.fq'])

    def test_reject_conflicting_or_missing_inputs(self):
        for args in [('reads.fq', 'reverse.fq', True), (None, 'reverse.fq', False),
                     ('None', None, True), ('', None, False), ('reads.fq', 'reads.fq', False)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                tpm.coverm_read_arguments(*args)

    def test_same_file_aliases_are_not_two_mates(self):
        with tempfile.TemporaryDirectory() as tmp:
            forward = Path(tmp) / 'reads.fq'
            forward.touch()
            for alias, link in [('symbolic.fq', os.symlink), ('hard.fq', os.link)]:
                reverse = Path(tmp) / alias
                link(forward, reverse)
                with self.subTest(alias=alias), self.assertRaises(ValueError):
                    tpm.coverm_read_arguments(str(forward), str(reverse))

    def test_main_dispatches_layout_and_stops_after_mapping_failure(self):
        # Include the archived '-2 None --interleaved' invocation, and a stale
        # BAM that must not be processed when CoverM fails.
        for mode in ('paired', 'interleaved', 'single'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                for name in ['contigs.fa', 'annotation.gff', 'forward reads.fq', 'reverse reads.fq', 'stale.bam']:
                    (root / name).touch()
                forward, reverse = str(root / 'forward reads.fq'), str(root / 'reverse reads.fq')
                args = ['-c', str(root/'contigs.fa'), '-g', str(root/'annotation.gff'),
                        '-o', str(root/'counts.tsv'), '--stats', str(root/'mapping.log'),
                        '--bwaFolder', tmp, '--sample_name', 'sample', '--rpkmExec', 'coverm',
                        '--bwaExec', 'bwa', '-1', forward, '--num_threads', '8']
                if mode == 'paired':
                    args += ['-2', reverse]
                elif mode == 'interleaved':
                    args += ['-2', 'None', '--interleaved']
                with patch.object(tpm.shutil, 'which', return_value='/mock/tool'), \
                     patch.object(tpm, 'runRPKMCommand', return_value=(1, 'mapping failed')) as run, \
                     patch.object(tpm, 'run_logged') as downstream, \
                     patch.object(tpm.gutils, 'eprintf'):
                    self.assertEqual(tpm.main(args), 1)
                    self.assertEqual(downstream.call_count, 1)
                    self.assertEqual(downstream.call_args.args[0][0], 'gff2gtf.py')
                cmd = shlex.split(run.call_args.kwargs['runcommand'])
                self.assertEqual(cmd[cmd.index('-t')+1], '8')
                if mode == 'paired':
                    self.assertEqual(cmd[cmd.index('-1')+1], forward)
                    self.assertEqual(cmd[cmd.index('-2')+1], reverse)
                else:
                    flag = '--interleaved' if mode == 'interleaved' else '--single'
                    self.assertEqual(cmd[cmd.index(flag)+1], forward)
                    self.assertNotIn('-2', cmd)
                self.assertIn('mapping failed', (root/'mapping.log').read_text())

    def test_counting_layout_and_sort_threads(self):
        for paired in (False, True):
            with self.subTest(paired=paired), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                for name in ('contigs.fa', 'annotation.gff', 'reads.fq', 'contigs.fa.reads.fq.bam',
                             'sample.sorted.bam', 'sample.old.bam'):
                    (root / name).touch()
                args = ['-c', str(root/'contigs.fa'), '-g', str(root/'annotation.gff'),
                    '-o', str(root/'counts.tsv'), '--stats', str(root/'mapping.log'),
                    '--bwaFolder', tmp, '--sample_name', 'sample', '--rpkmExec', 'coverm',
                    '--bwaExec', 'bwa', '-1', str(root/'reads.fq'), '--num_threads', '8']
                if paired:
                    args += ['--interleaved']
                with patch.object(tpm.shutil, 'which', return_value='/mock/tool'), \
                     patch.object(tpm, 'runRPKMCommand', return_value=(0, '')), \
                     patch.object(tpm, 'run_logged') as run:
                    self.assertEqual(tpm.main(args), 0)
                commands = [call.args[0] for call in run.call_args_list]
                sort = next(c for c in commands if c[0] == 'samtools')
                counts = next(c for c in commands if c[0] == 'featureCounts')
                self.assertEqual(sort[sort.index('-@')+1], '7')
                self.assertEqual(sort[-1], str(root/'contigs.fa.reads.fq.bam'))
                self.assertEqual(sort[sort.index('-o')+1], str(root/'sample.sorted.bam'))
                self.assertEqual(counts[counts.index('-T')+1], '8')
                self.assertEqual('-p' in counts, paired)
                # Existing sorted/stale BAMs cannot substitute for missing CoverM output.
                (root/'contigs.fa.reads.fq.bam').unlink()
                with patch.object(tpm.shutil, 'which', return_value='/mock/tool'), \
                     patch.object(tpm, 'runRPKMCommand', return_value=(0, '')), \
                     patch.object(tpm, 'run_logged') as run:
                    with self.assertRaisesRegex(RuntimeError, 'expected BAM'):
                        tpm.main(args)
                    self.assertEqual(run.call_count, 1)  # GFF preparation only

    def test_pipeline_does_not_serialize_missing_reverse(self):
        # Access the class wrapped by the project's Singleton decorator.
        cls = jobscreator.ContextCreator
        creator = object.__new__(cls)
        creator.configs = SimpleNamespace(RPKM_EXECUTABLE='coverm', BWA_EXECUTABLE='bwa',
                                         NUM_CPUS='8', RPKM_CALCULATION='MetaPathways_tpm')
        creator.params = SimpleNamespace(get=lambda *a: 'yes')
        sample = SimpleNamespace(fq_files=[['reads.fq', 'None'], True], bwa_folder='/tmp/bwa',
                                 genbank_dir='/tmp/genbank/', sample_name='sample',
                                 preprocessed_dir='/tmp/preprocessed', output_results_rpkm_dir='/tmp/rpkm')
        cmd = shlex.split(creator.create_rpkm_cmd(sample)[0].commands[0])
        self.assertEqual(cmd[cmd.index('--num_threads')+1], '8')
        self.assertIn('--interleaved', cmd)
        self.assertNotIn('-2', cmd)
        self.assertNotIn('None', cmd)
        sample.fq_files = [['R1.fq', 'R2.fq'], False]
        cmd = shlex.split(creator.create_rpkm_cmd(sample)[0].commands[0])
        self.assertEqual(cmd[cmd.index('-1')+1], 'R1.fq')
        self.assertEqual(cmd[cmd.index('-2')+1], 'R2.fq')
        self.assertNotIn('--interleaved', cmd)

    def test_orf_and_trna_commands_pass_configured_threads(self):
        cls = jobscreator.ContextCreator
        creator = object.__new__(cls)
        creator.params = SimpleNamespace(get=lambda *a: 'yes')
        sample = SimpleNamespace(sample_name='sample', preprocessed_dir='/tmp/input',
            orf_prediction_dir='/tmp/orfs/', output_results_tRNA_dir='/tmp/trna')
        for threads in (1, 3, 8):
            creator.configs = SimpleNamespace(NUM_CPUS=threads,
                ORF_PREDICTION='orf', PRODIGAL_EXECUTABLE='pprodigal',
                SCAN_tRNA='trna', SCAN_tRNA_EXECUTABLE='ptRNAscan.py')
            for method, flag in ((creator.create_orf_prediction_cmd, '--nthreads'),
                                 (creator.create_tRNA_scan_statistics, '-t')):
                with self.subTest(threads=threads, flag=flag):
                    cmd = shlex.split(method(sample)[0].commands[0])
                    self.assertEqual(cmd[cmd.index(flag)+1], str(threads))

    def test_barrnap_and_rrna_blast_use_the_same_budget(self):
        cls = jobscreator.ContextCreator
        creator = object.__new__(cls)
        creator.params = SimpleNamespace(get=lambda group, key, **kw:
            ['ssu'] if key == 'rRNA_refdbs' else 'yes')
        creator.configs = SimpleNamespace(NUM_CPUS=3, PARSE_FUNC_SEARCH='parse',
            SCAN_rRNA='rrna', REFDBS='/tmp/db')
        sample = SimpleNamespace(sample_name='sample', preprocessed_dir='/tmp/input',
            orf_prediction_dir='/tmp/orfs/', blast_results_dir='/tmp/blast',
            output_results_rRNA_dir='/tmp/rrna/')
        with patch.object(jobscreator.shutil, 'which', side_effect=lambda name: name):
            contexts = creator.create_scan_rRNA_seqs_cmd(sample)
        cmd = shlex.split(contexts[0].commands[0])
        self.assertEqual(cmd[cmd.index('--threads')+1], '3')
        cmd = shlex.split(contexts[1].commands[1])
        self.assertEqual(cmd[cmd.index('-num_threads')+1], '3')

    def test_trnascan_workers_do_not_create_nested_thread_pools(self):
        import importlib.util
        spec = importlib.util.spec_from_file_location('ptrna_test',
            Path(__file__).resolve().parents[1] / 'dev/ptRNAscan.py')
        module = importlib.util.module_from_spec(spec)
        # This worker test does not parse FASTA or require pyfastx.
        with patch.dict('sys.modules', pyfastx=SimpleNamespace()):
            spec.loader.exec_module(module)
        args = SimpleNamespace(**dict.fromkeys(('bacterial', 'archaeal', 'mito',
            'general', 'genomic', 'eukaryotic', 'infernal', 'max', 'legacy',
            'cove', 'nopseudo', 'quiet'), False))
        module.file_parts_queue.put(('chunk.fasta', {'-o': 'chunk.txt'}))
        with patch.object(module.subprocess, 'Popen') as run:
            module.run_tRNAscan(0, args)
        cmd = run.call_args.args[0]
        self.assertEqual(cmd[cmd.index('--thread')+1], '1')


if __name__ == '__main__':
    unittest.main()
