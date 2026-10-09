"""Exercise overlapping publishers without relying on shared filesystem locks."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch

from metapathways.nf_worker import atomic_json


class AtomicJsonTests(unittest.TestCase):
    def test_overlapping_writers_publish_complete_independent_records(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'cache.json'
            barrier = threading.Barrier(8)
            replace = Path.replace
            staged = []

            def simultaneous_replace(source, destination):
                staged.append(source)
                barrier.wait(timeout=10)
                return replace(source, destination)

            records = [{'writer': i, 'payload': str(i) * 10000} for i in range(8)]
            with patch.object(Path, 'replace', simultaneous_replace):
                with ThreadPoolExecutor(max_workers=8) as pool:
                    list(pool.map(lambda record: atomic_json(target, record), records))
            self.assertEqual(len(set(staged)), 8)
            self.assertIn(json.loads(target.read_text()), records)
            self.assertEqual(list(Path(directory).iterdir()), [target])

    def test_failed_publication_preserves_previous_record_and_cleans_temp(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / 'cache.json'
            atomic_json(target, {'old': True})
            with patch.object(Path, 'replace', side_effect=OSError('publication failed')):
                with self.assertRaisesRegex(OSError, 'publication failed'):
                    atomic_json(target, {'new': True})
            self.assertEqual(json.loads(target.read_text()), {'old': True})
            self.assertEqual(list(Path(directory).iterdir()), [target])


if __name__ == '__main__':
    unittest.main()
