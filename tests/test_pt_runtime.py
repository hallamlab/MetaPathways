"""Runtime selection must respect capabilities, host context and explicit choices."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from metapathways import pt_setup as setup


class RuntimeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.env = self.root/'conda/apptainer'
        self.site = self.root/'site/apptainer'
        self.image = self.root/'licensed.sif'
        self.image.write_bytes(b'fixture')
        for path in (self.env, self.site):
            path.parent.mkdir()
            path.write_text('#!/bin/sh\nexit 0\n')
            path.chmod(0o755)
        setup._SELECTED.clear()
        self.addCleanup(setup._SELECTED.clear)
        env = patch.dict(os.environ, {'XDG_CONFIG_HOME': str(self.root/'config')})
        env.start()
        self.addCleanup(env.stop)
        candidates = patch.object(setup, 'runtime_candidates', return_value=[str(self.env), str(self.site)])
        self.candidates = candidates.start()
        self.addCleanup(candidates.stop)

    def select(self, image=None):
        return setup.select_runtime(self.root/'logs', image=image)

    def test_environment_preferred_when_working(self):
        with patch.object(setup, 'host_check') as probe:
            self.assertEqual(self.select(), str(self.env))
        self.assertEqual(probe.call_count, 1)

    def test_fallback_is_saved_and_reprobed_in_a_new_process(self):
        with patch.object(setup, 'host_check', side_effect=[RuntimeError('namespace denied'), None]) as probe:
            self.assertEqual(self.select(), str(self.site))
        self.assertEqual(probe.call_count, 2)
        report = json.loads((self.root/'logs/runtime-selection.json').read_text())
        self.assertEqual([x['status'] for x in report['attempts']], ['failed', 'passed'])
        setup._SELECTED.clear()
        with patch.object(setup, 'host_check') as probe:
            self.assertEqual(self.select(), str(self.site))
        self.assertEqual(probe.call_args.args[0], str(self.site))

    def test_all_failures_are_retained(self):
        with patch.object(setup, 'host_check', side_effect=RuntimeError('denied')):
            with self.assertRaisesRegex(RuntimeError, 'Tested 2 runtime'):
                self.select()
        report = json.loads((self.root/'logs/runtime-selection.json').read_text())
        self.assertIsNone(report['selected'])
        self.assertEqual(len(report['attempts']), 2)
        self.assertFalse(setup.runtime_record().exists())

    def test_existing_image_does_not_require_build_permission(self):
        with patch.object(setup, 'host_check', side_effect=AssertionError('must not build')), \
             patch.object(setup, 'execution_check') as probe:
            self.assertEqual(self.select(self.image), str(self.env))
            self.assertEqual(self.select(self.image), str(self.env))
        self.assertEqual(probe.call_count, 1)

    def test_execution_success_never_certifies_build_support(self):
        with patch.object(setup, 'execution_check'):
            self.select(self.image)
        with patch.object(setup, 'host_check', side_effect=RuntimeError('building denied')) as probe:
            with self.assertRaisesRegex(RuntimeError, 'No usable Apptainer found for build'):
                self.select()
        self.assertEqual(probe.call_count, 2)

    def test_build_validation_uses_selected_executable(self):
        from metapathways import pt_container
        with patch.object(pt_container.subprocess, 'run', return_value=
                          subprocess.CompletedProcess([], 0, 'MP-PT-READY\n')) as run, \
             patch.object(setup, 'select_runtime', side_effect=AssertionError('must retain selected runtime')):
            pt_container.validate(self.image, self.root/'state', executable=str(self.site))
        self.assertEqual(run.call_args.args[0][0], str(self.site))

    def test_saved_builder_is_preferred_for_execution(self):
        with patch.object(setup, 'host_check', side_effect=[RuntimeError('denied'), None]):
            self.select()
        with patch.object(setup, 'execution_check') as probe:
            self.assertEqual(self.select(self.image), str(self.site))
        self.assertEqual(probe.call_args.args[0], str(self.site))

    def test_changed_image_is_rechecked(self):
        with patch.object(setup, 'execution_check') as probe:
            self.select(self.image)
            self.image.write_bytes(b'a different SIF')
            self.select(self.image)
        self.assertEqual(probe.call_count, 2)

    def test_compute_host_cannot_reuse_login_host_success(self):
        with patch.object(setup, 'host_check') as probe:
            with patch.object(setup, 'runtime_context', return_value={'host': 'login'}):
                self.select()
            with patch.object(setup, 'runtime_context', return_value={'host': 'compute'}):
                self.select()
        self.assertEqual(probe.call_count, 2)

    def test_stale_preference_falls_back(self):
        record = setup.runtime_record()
        record.parent.mkdir(parents=True)
        record.write_text(json.dumps({'build': str(self.site)}))
        with patch.object(setup, 'host_check', side_effect=[RuntimeError('site changed'), None]):
            self.assertEqual(self.select(), str(self.env))

    def test_explicit_candidate_discovery_does_not_fall_back(self):
        with patch.dict(os.environ, {'METAPATHWAYS_APPTAINER': str(self.site)}):
            self.assertEqual(REAL_CANDIDATES(), [str(self.site)])

    def test_missing_explicit_runtime_has_actionable_error(self):
        with patch.dict(os.environ, {'METAPATHWAYS_APPTAINER': str(self.root/'missing')}):
            with self.assertRaisesRegex(RuntimeError, 'METAPATHWAYS_APPTAINER is not executable'):
                REAL_CANDIDATES()

    def test_execution_probe_uses_containment_and_private_writable_storage(self):
        with patch.object(setup.subprocess, 'run') as run:
            setup.execution_check(str(self.site), self.image, self.root/'exec')
        command = run.call_args.args[0]
        self.assertIn('--containall', command)
        self.assertIn('--home', command)
        self.assertIn('--workdir', command)
        self.assertNotIn('--fakeroot', command)
        self.assertFalse(list((self.root/'exec').glob('.pt-exec-*')))


REAL_CANDIDATES = setup.runtime_candidates

if __name__ == '__main__':
    unittest.main()
