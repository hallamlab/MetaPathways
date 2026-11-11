"""Readiness failures must precede licensed downloads; public cache is recoverable."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from metapathways import pt_container as pt, pt_setup as setup


class SetupTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(prefix='mp readiness ')
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def test_build_command_supports_modern_and_older_runtimes(self):
        for help_text, supported in [('--fakeroot --mksquashfs-args', True), ('--fakeroot --sandbox', False)]:
            with self.subTest(supported=supported), patch.object(
                    setup.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, help_text)):
                command = setup.build_command('/site/apptainer', 'output.sif', 'input.def', 2)
            self.assertEqual('--mksquashfs-args' in command, supported)
            self.assertEqual(command[-2:], ['output.sif', 'input.def'])
            self.assertIn('--fakeroot', command)
            if supported:
                self.assertIn('-processors 2', command)

    def test_runtime_help_failure_is_not_treated_as_legacy_support(self):
        with patch.object(setup.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, ['apptainer'])):
            with self.assertRaises(subprocess.CalledProcessError):
                setup.build_command('/site/apptainer', 'output.sif', 'input.def', 2)

    def test_check_needs_neither_installer_nor_nextflow(self):
        with patch.object(pt.shutil, 'which', return_value='/env/bin/apptainer'), \
             patch.object(setup, 'select_runtime') as check, patch.object(pt.nextflow, 'launch') as run:
            pt.main(['--check', '-o', str(self.root)])
        check.assert_called_once_with(self.root/'logs/readiness', image=None)
        run.assert_not_called()

    def test_denied_host_stops_before_base_and_licensed_downloads(self):
        installer = self.root/'installer'
        installer.write_bytes(b'installer')
        with patch.object(pt.shutil, 'which', return_value='/env/bin/apptainer'), \
             patch.object(setup, 'select_runtime', side_effect=RuntimeError('host denied')), \
             patch.object(setup, 'prepared_base') as base, patch.object(pt, 'download_patches') as patches:
            with self.assertRaisesRegex(RuntimeError, 'host denied'):
                pt.build(installer, self.root/'image.sif', pt.digest(installer), 1)
        base.assert_not_called()
        patches.assert_not_called()

    def test_probe_reports_permission_failure_and_cleans_temporary_files(self):
        real_run = subprocess.run
        def run(command, **kwargs):
            if command[-1] == '--help':
                return subprocess.CompletedProcess(command, 0, '--mksquashfs-args')
            if command[0] == 'ldd':
                return real_run(command, **kwargs)
            kwargs['stdout'].write('Failed to create mount namespace: requires privileges\n')
            kwargs['stdout'].flush()
            raise subprocess.CalledProcessError(1, command)
        with patch.object(setup.subprocess, 'run', side_effect=run):
            with self.assertRaisesRegex(RuntimeError, 'host denied container setup'):
                setup.host_check('/env/bin/apptainer', self.root)
        self.assertIn('/env/bin/apptainer', (self.root/'apptainer-check.log').read_text())
        self.assertFalse(list(self.root.glob('.pt-check-*')))

    def test_probe_builds_then_executes_without_inherited_binds(self):
        real_run = subprocess.run
        calls = []
        def run(command, **kwargs):
            if command[-1] == '--help':
                return subprocess.CompletedProcess(command, 0, '--mksquashfs-args')
            if command[0] == 'ldd':
                return real_run(command, **kwargs)
            calls.append(command)
            self.assertNotIn('APPTAINER_BIND', kwargs['env'])
            self.assertNotIn('SINGULARITY_BINDPATH', kwargs['env'])
            return subprocess.CompletedProcess(command, 0)
        with patch.dict(os.environ, {'APPTAINER_BIND': '/bad', 'SINGULARITY_BINDPATH': '/bad'}), \
             patch.object(setup.subprocess, 'run', side_effect=run):
            setup.host_check('/env/bin/apptainer', self.root)
        self.assertEqual([call[1] for call in calls], ['build', 'exec'])
        self.assertIn('--fakeroot', calls[0])
        self.assertIn('--workdir', calls[1])

    def test_cache_reuses_valid_image_and_repairs_damage(self):
        def pull(command, **kwargs):
            Path(command[2]).write_bytes(b'public base')
        with patch.dict(os.environ, {'XDG_CACHE_HOME': str(self.root)}), \
             patch.object(setup.subprocess, 'run', side_effect=pull) as run:
            image = setup.prepared_base('/env/bin/apptainer')
            self.assertEqual(setup.prepared_base('/env/bin/apptainer'), image)
            self.assertEqual(run.call_count, 1)
            image.write_bytes(b'damaged')
            setup.prepared_base('/env/bin/apptainer')
            self.assertEqual(run.call_count, 2)
            image.with_suffix('.json').write_text('{incomplete')
            setup.prepared_base('/env/bin/apptainer')
            self.assertEqual(run.call_count, 3)
            self.assertEqual(image.read_bytes(), b'public base')

    def test_failed_pull_publishes_neither_image_nor_receipt(self):
        with patch.dict(os.environ, {'XDG_CACHE_HOME': str(self.root)}), \
             patch.object(setup.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, ['apptainer'])):
            with self.assertRaises(subprocess.CalledProcessError):
                setup.prepared_base('/env/bin/apptainer')
        self.assertFalse(list(self.root.rglob('*.sif')))
        self.assertFalse(list(self.root.rglob('*.json')))
        self.assertFalse(list(self.root.rglob('.download-*')))

    def test_explicit_base_does_not_download_and_recipe_skips_apt(self):
        base = self.root/'base.sif'
        base.write_bytes(b'public base')
        with patch.object(setup.subprocess, 'run') as run:
            self.assertEqual(setup.prepared_base('/env/bin/apptainer', base), base)
        run.assert_not_called()
        recipe, provenance = pt.build_recipe(base)
        self.assertIn('From: dependency-base.sif', recipe)
        self.assertNotIn('apt-get', recipe)
        self.assertEqual(provenance['sha256'], setup.checksum(base))


if __name__ == '__main__':
    unittest.main()
