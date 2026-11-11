import subprocess
import sys
from pathlib import Path
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import retry_network


class RetryTests(unittest.TestCase):
    def test_network_errors_and_conflicts(self):
        for message in ('ServerError: HTTP 520', 'Unable to find the server at api.launchpad.net'):
            self.assertTrue(retry_network.transient(message))
        for message in ('401 Unauthorized', '403 Forbidden', '409 already exists',
                        '520 followed by already exists', 'invalid package', 'checksum mismatch'):
            self.assertFalse(retry_network.transient(message))

    def test_retry_is_bounded(self):
        with patch('retry_network.time.sleep') as sleep:
            with self.assertRaises(subprocess.CalledProcessError):
                retry_network.run([sys.executable, '-c', 'print("HTTP 520"); exit(1)'])
            self.assertEqual([c.args[0] for c in sleep.call_args_list], [15, 30])

    def test_auth_failure_does_not_retry(self):
        with patch('retry_network.time.sleep') as sleep:
            with self.assertRaises(subprocess.CalledProcessError):
                retry_network.run([sys.executable, '-c', 'print("401 Unauthorized"); exit(1)'])
            sleep.assert_not_called()

    def test_success_does_not_retry(self):
        with patch('retry_network.time.sleep') as sleep:
            retry_network.run([sys.executable, '-c', 'print("ok")'])
            sleep.assert_not_called()
