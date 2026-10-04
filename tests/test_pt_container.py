"""Container lifecycle tests; licensed installation is a separate integration check."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from metapathways import pt_container as pt


class ContainerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mp pt test ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.installer = self.root / 'pathway-tools-29.5-install'
        self.installer.write_bytes(b'fixture installer')
        self.image = self.root / 'pt.sif'
        patch_download = patch.object(pt, 'download_patches', return_value={'files': [], 'version': '29.5'})
        self.real_download_patches = pt.download_patches
        self.download_patches = patch_download.start()
        self.addCleanup(patch_download.stop)

    def fake_build(self, command, **kwargs):
        self.assertIn('--fakeroot', command)
        self.assertIn('-processors 2', command)
        Path(command[-2]).write_bytes(b'fixture SIF')
        return subprocess.CompletedProcess(command, 0)

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_publish_after_validation(self, _):
        with patch.object(pt.subprocess, 'run', side_effect=self.fake_build), patch.object(pt, 'validate', return_value='MP-PT-READY\n') as validate:
            pt.build(self.installer, self.image, pt.digest(self.installer), 2)
        validate.assert_called_once()
        metadata = json.loads(Path(str(self.image) + '.json').read_text())
        self.assertEqual(metadata['image_sha256'], pt.digest(self.image))
        self.assertEqual(metadata['installer_sha256'], pt.digest(self.installer))
        self.assertEqual(metadata['official_patches']['version'], '29.5')
        self.assertFalse(list(self.root.glob('.pt-build-*')))

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_failed_validation_preserves_previous_image(self, _):
        self.image.write_bytes(b'previous valid image')
        with patch.object(pt.subprocess, 'run', side_effect=self.fake_build), patch.object(pt, 'validate', side_effect=RuntimeError('startup failed')):
            with self.assertRaisesRegex(RuntimeError, 'startup failed'):
                pt.build(self.installer, self.image, pt.digest(self.installer), 2)
        self.assertEqual(self.image.read_bytes(), b'previous valid image')
        self.assertFalse(Path(str(self.image) + '.json').exists())
        self.assertFalse(list(self.root.glob('.pt-build-*')))

    def test_changed_installer_rejected_before_execution(self):
        with patch.object(pt.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'Installer changed'):
                pt.build(self.installer, self.image, 'wrong-checksum', 2)
        run.assert_not_called()

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_runtime_isolation(self, _):
        cmd = pt.exec_command(self.image, self.root / 'private state', ['sh', '-c', 'echo hello'])
        for flag in ['--containall', '--cleanenv']:
            self.assertIn(flag, cmd)
        self.assertIn(str(self.root / 'private state') + ':/data', cmd)
        self.assertEqual(cmd[cmd.index('--home') + 1], str(self.root / 'private state') + ':/data')
        self.assertEqual(cmd[cmd.index('--pwd') + 1], '/data')
        self.assertEqual(cmd[-3:], ['sh', '-c', 'echo hello'])

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_validation_requires_lisp_marker(self, _):
        with patch.object(pt.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, 'no startup')):
            with self.assertRaisesRegex(RuntimeError, 'validation marker'):
                pt.validate(self.image, self.root / 'state')

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_validation_preserves_startup_error(self, _):
        with patch.object(pt.subprocess, 'run', return_value=subprocess.CompletedProcess([], 1, 'missing shared library')):
            with self.assertRaisesRegex(RuntimeError, 'missing shared library'):
                pt.validate(self.image, self.root / 'state')

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_validation_marker_on_own_line_after_echo(self, _):
        result = subprocess.CompletedProcess([], 0, 'CMD: (format t "~%MP-PT-READY~%")\nMP-PT-READY\n')
        with patch.object(pt.subprocess, 'run', return_value=result) as run:
            pt.validate(self.image, self.root / 'state')
        self.assertIn('"~%MP-PT-READY~%"', run.call_args.args[0][-1])

    def test_registry_and_override(self):
        self.image.touch()
        with patch.dict(os.environ, {'XDG_CONFIG_HOME': str(self.root)}, clear=True):
            self.assertIsNone(pt.registered_image())
            pt.save_json(pt.registry_path(), {'image': str(self.image)})
            self.assertEqual(pt.registered_image(), str(self.image))
            with patch.dict(os.environ, {'METAPATHWAYS_PTOOLS_IMAGE': str(self.root / 'missing.sif')}):
                with self.assertRaisesRegex(ValueError, 'image is missing'):
                    pt.registered_image()

    def test_dryrun_without_runtimes_never_registers(self):
        with patch.dict(os.environ, {'XDG_CONFIG_HOME': str(self.root)}, clear=True), patch.object(pt.shutil, 'which', return_value=None), contextlib.redirect_stdout(io.StringIO()):
            pt.main(['-i', str(self.installer), '-o', str(self.root / 'images'), '--dryrun'])
            self.assertFalse(pt.registry_path().exists())
        plans = list(self.root.glob('images/logs/build_pt/*/tasks.json'))
        self.assertEqual(len(plans), 1)
        task = json.loads(plans[0].read_text())[0]
        self.assertEqual(task['cpus'], 2)
        self.assertFalse(task['adopt_existing'])

    def test_invalid_tag_rejected(self):
        with self.assertRaisesRegex(ValueError, 'PGDB tag'):
            pt.run_pgdb(self.image, self.root, self.root / 'out', "bad') (exit)")

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_missing_sequence_source_stops_before_pathway_tools_and_retains_diagnostics(self, _):
        inputs = self.root/'input'
        inputs.mkdir()
        (inputs/'0.pf').write_text('ID\tgene1\n//\n')
        output = self.root/'failed-input'
        with patch.object(pt.subprocess, 'run') as run:
            with self.assertRaisesRegex(ValueError, 'sequence-backed Pathway Tools input'):
                pt.run_pgdb(self.image, inputs, output, 'sample', sample_output=self.root/'sample')
        run.assert_not_called()
        record = json.loads(next(output.glob('diagnostics/*/execution.json')).read_text())
        self.assertEqual(record['stage'], 'input preparation')
        self.assertEqual(record['status'], 'FAILED')
        self.assertFalse(list(self.root.glob('.pt-run-*')))

    def test_build_with_mpdb_plans_export_after_image_without_running_it(self):
        with patch.object(pt.nextflow, 'launch') as launch:
            pt.main(['-i', str(self.installer), '-o', str(self.root/'images'), '-d', str(self.root/'MPDB'), '-a', 'blast', '--dryrun'])
        tasks = launch.call_args.args[0]
        self.assertEqual(len(tasks), 2)
        self.assertEqual(tasks[1]['dependencies'], [tasks[0]['id']])
        self.assertEqual(tasks[1]['inputs'], [tasks[0]['outputs'][0]])
        self.assertIn(str(self.root/'MPDB/functional/formatted/metacyc.pdb'), tasks[1]['outputs'])
        self.download_patches.assert_not_called()

    def test_patch_snapshot_only_accepts_official_listing_filenames(self):
        listing = b'<a href="p123.fasl">patch</a><a href="../evil.lisp">bad</a><a href="https://other/evil.fasl">bad</a><a href="p124.tar.gz">archive</a>'
        contents = [listing, b'official binary', b'official archive']
        with patch.object(pt, 'build_opener') as opener:
            opener.return_value.open.side_effect = [io.BytesIO(x) for x in contents]
            manifest = self.real_download_patches('29.5', self.root / 'patches')
        self.assertEqual([x['name'] for x in manifest['files']], ['p123.fasl', 'p124.tar.gz'])
        self.assertEqual(manifest['files'][0]['sha256'], pt.digest(self.root/'patches/files/p123.fasl'))
        self.assertTrue(all(c.args[0].startswith(manifest['source']) for c in opener.return_value.open.call_args_list))

    def test_failed_patch_download_does_not_build_or_publish(self):
        self.download_patches.side_effect = RuntimeError('official feed unavailable')
        with patch.object(pt.shutil, 'which', return_value='/bin/apptainer'), patch.object(pt.subprocess, 'run') as run:
            with self.assertRaisesRegex(RuntimeError, 'official feed unavailable'):
                pt.build(self.installer, self.image, pt.digest(self.installer), 2)
        run.assert_not_called()
        self.assertFalse(self.image.exists())

    def test_empty_vendor_listing_and_redirect_are_rejected(self):
        with patch.object(pt, 'build_opener') as opener:
            opener.return_value.open.return_value = io.BytesIO(b'<html>login required</html>')
            with self.assertRaisesRegex(RuntimeError, 'no recognized patch'):
                self.real_download_patches('29.5', self.root / 'empty')
        with self.assertRaisesRegex(ValueError, 'redirected'):
            pt._NoRedirect().redirect_request(None, None, 302, '', {}, 'https://other/patch')

    def test_release_detection_is_explicit(self):
        self.assertEqual(pt.installer_version(self.installer), '29.5')
        self.assertEqual(pt.installer_version('renamed-installer', '29.5'), '29.5')
        for name, version in [('renamed', None), (str(self.installer), '28.5'), ('renamed', '../29.5')]:
            with self.assertRaises(ValueError):
                pt.installer_version(name, version)

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_validation_rejects_blast_warning_even_with_startup_marker(self, _):
        result = subprocess.CompletedProcess([], 0, 'blastall or blastp could not be located\nMP-PT-READY\n')
        with patch.object(pt.subprocess, 'run', return_value=result):
            with self.assertRaisesRegex(RuntimeError, 'could not locate'):
                pt.validate(self.image, self.root / 'state')

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_pgdb_failure_preserves_private_error_log(self, _):
        inputs = self.root / 'input'
        inputs.mkdir()
        (inputs / '0.pf').write_text('ID\tgene1\n')
        output = self.root / 'failed'
        def fail(command, **kwargs):
            state = Path(command[command.index('--home')+1].rsplit(':', 1)[0])
            (state / 'stage.txt').write_text('build\n')
            (state / 'input/pathologic.log').write_text('Specific Lisp failure\n')
            kb = state / 'ptools-local/pgdbs/user/samplecyc/1.0/kb'
            kb.mkdir(parents=True)
            (kb / 'samplebase.ocelot').write_text('saved checkpoint')
            raise subprocess.CalledProcessError(255, command)
        with patch.object(pt.subprocess, 'run', side_effect=fail), contextlib.redirect_stdout(io.StringIO()) as console:
            with self.assertRaises(subprocess.CalledProcessError):
                pt.run_pgdb(self.image, inputs, output, 'sample')
        records = list(output.glob('diagnostics/*/execution.json'))
        self.assertEqual(len(records), 1)
        record = json.loads(records[0].read_text())
        self.assertEqual((record['status'], record['stage'], record['exit_code']), ('FAILED', 'build', 255))
        self.assertEqual((records[0].parent/'input/pathologic.log').read_text(), 'Specific Lisp failure\n')
        self.assertIn('Specific Lisp failure', console.getvalue())
        self.assertEqual((Path(record['retained_pgdbs'])/'samplecyc/1.0/kb/samplebase.ocelot').read_text(), 'saved checkpoint')
        self.assertTrue((records[0].parent/'input/0.pf').exists())
        self.assertFalse(list(self.root.glob('.pt-run-*')))
        self.assertFalse((inputs / 'pathologic.log').exists())

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_pgdb_success_preserves_log_and_publishes_archive(self, _):
        inputs = self.root / 'input'
        inputs.mkdir()
        output = self.root / 'success'
        def success(command, **kwargs):
            script = command[command.index('-ec') + 1]
            self.assertIn('mkdir -p /data/blastdb /tmp/.X11-unix', script)
            self.assertEqual(script.count('-nolisten local'), 2)
            self.assertIn('-e /data/build-xvfb.log', script)
            self.assertIn('-e /data/export-xvfb.log', script)
            state = Path(command[command.index('--home')+1].rsplit(':', 1)[0])
            (state / 'build-xvfb.log').write_text('build display diagnostics\n')
            (state / 'export-xvfb.log').write_text('export display diagnostics\n')
            (state / 'stage.txt').write_text('archive\n')
            (state / 'input/pathologic.log').write_text('Pathologic finished\n')
            (state / 'output/samplecyc.tar.bz2').write_bytes(b'fixture archive')
        with patch.object(pt.subprocess, 'run', side_effect=success):
            pt.run_pgdb(self.image, inputs, output, 'sample')
        record = next(output.glob('diagnostics/*/execution.json'))
        self.assertEqual(json.loads(record.read_text())['status'], 'SUCCESS')
        self.assertEqual((output/'samplecyc.tar.bz2').read_bytes(), b'fixture archive')
        self.assertTrue((record.parent/'input/pathologic.log').is_file())
        self.assertEqual((record.parent/'build-xvfb.log').read_text(), 'build display diagnostics\n')
        self.assertEqual((record.parent/'export-xvfb.log').read_text(), 'export display diagnostics\n')
        self.assertFalse(list(self.root.glob('.pt-run-*')))

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_missing_archive_is_failure_with_diagnostics(self, _):
        inputs = self.root / 'input'
        inputs.mkdir()
        output = self.root / 'no-archive'
        with patch.object(pt.subprocess, 'run'), self.assertRaisesRegex(RuntimeError, 'no PGDB archive'):
            pt.run_pgdb(self.image, inputs, output, 'sample')
        record = next(output.glob('diagnostics/*/execution.json'))
        self.assertEqual(json.loads(record.read_text())['status'], 'FAILED')

    def test_registered_image_plans_single_cpu_entities(self):
        from metapathways import pipeline
        self.image.touch()
        (self.root / 'magsplitter/results/MAG_1').mkdir(parents=True)
        (self.root / 'magsplitter/results/non_binned').mkdir()
        with patch.object(pipeline.sys, 'argv', ['metapathways', 'ptools', '-o', str(self.root)]), patch.object(pt, 'registered_image', return_value=str(self.image)), patch.object(pt.nextflow, 'launch') as launch:
            pipeline.ptools()
        tasks = launch.call_args.args[0]
        self.assertEqual(len(tasks), 2)
        self.assertEqual([t['cpus'] for t in tasks], [1, 1])
        self.assertEqual([t['allow_failure'] for t in tasks], [False, True])
        self.assertTrue(all('--image' in t['commands'][0] for t in tasks))

    def test_no_transport_and_single_entity_are_forwarded(self):
        from metapathways import pipeline
        self.image.touch()
        (self.root / 'magsplitter/results/MAG_1').mkdir(parents=True)
        with patch.object(pipeline.sys, 'argv', ['metapathways', 'ptools', '-o', str(self.root),
                '--entity', 'community', '--no_transport_inference', '--taxprune', '--taxon_id', '131567']), \
             patch.object(pt, 'registered_image', return_value=str(self.image)), \
             patch.object(pt.nextflow, 'launch') as launch:
            pipeline.ptools()
        tasks = launch.call_args.args[0]
        self.assertEqual(len(tasks), 1)
        self.assertIn('--no_transport_inference', tasks[0]['commands'][0])
        self.assertIn('--taxon_id 131567', tasks[0]['commands'][0])
        self.assertIn('--taxprune', tasks[0]['commands'][0])

    @patch.object(pt.shutil, 'which', return_value='/bin/apptainer')
    def test_transport_switch_controls_tip_argument(self, _):
        inputs = self.root/'switch-input'
        inputs.mkdir()
        for enabled in (True, False):
            def run(command, **kwargs):
                self.assertEqual('-tip' in command, enabled)
                state = Path(command[command.index('--home')+1].rsplit(':', 1)[0])
                (state/'output/samplecyc.tar.bz2').write_bytes(b'archive')
            with patch.object(pt.subprocess, 'run', side_effect=run):
                pt.run_pgdb(self.image, inputs, self.root/str(enabled), 'sample', transport_inference=enabled)

    def test_worker_rechecks_outputs_and_directory_inputs(self):
        from metapathways import nf_worker as worker
        source = self.root / 'source'
        source.mkdir()
        (source / 'input').write_text('first')
        output = self.root / 'result'
        t = pt.nextflow.task('fixture', 'fixture', ['mock-tool'], [str(source)], [str(output)], adopt_existing=False)
        t['receipt'] = str(self.root / 'durable.json')
        def run(*args, **kwargs):
            output.write_text((source / 'input').read_text())
        previous_cwd = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(worker.subprocess, 'run', side_effect=run) as invoke:
                self.assertEqual(worker.execute(t), 0)
                self.assertEqual(worker.execute(t), 0)
                self.assertEqual(invoke.call_count, 1)
                (source / 'input').write_text('changed input')
                self.assertEqual(worker.execute(t), 0)
                self.assertEqual(output.read_text(), 'changed input')
                self.assertEqual(invoke.call_count, 2)
                output.unlink()
                self.assertEqual(worker.execute(t), 0)
                self.assertEqual(invoke.call_count, 3)
                t['cache_version'] = 'private-fast-scratch-v1'
                self.assertEqual(worker.execute(t), 0)
                self.assertEqual(invoke.call_count, 4)
        finally:
            os.chdir(previous_cwd)

    def test_optional_failure_recorded_but_community_failure_fatal(self):
        from metapathways import nf_worker as worker
        t = pt.nextflow.task('fixture', 'fixture', ['mock-tool'], outputs=[str(self.root / 'missing')], adopt_existing=False)
        t['receipt'] = str(self.root / 'durable.json')
        previous_cwd = Path.cwd()
        try:
            os.chdir(self.root)
            with patch.object(worker.subprocess, 'run', side_effect=RuntimeError('PGDB failed')), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(worker.execute(t), 1)
                t['allow_failure'] = True
                self.assertEqual(worker.execute(t), 0)
            self.assertEqual(json.loads(Path(t['receipt']).read_text())['status'], 'FAILED')
        finally:
            os.chdir(previous_cwd)


if __name__ == '__main__':
    unittest.main()
