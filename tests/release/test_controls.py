import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from check_release_controls import validate


class ControlTests(unittest.TestCase):
    def test_every_destination_combination_requires_manual_tagged_publication(self):
        for flags in itertools.product((False, True), repeat=3):
            selections = dict(zip(('anaconda', 'quay', 'github'), flags))
            with self.subTest(flags=flags):
                validate('workflow_dispatch', 'hallamlab/MetaPathways', 'tag', 'v3.5.2', **selections)
                validate('workflow_dispatch', 'hallamlab/MetaPathways', 'branch', 'dev', 'v3.5.2', **selections)
                if any(flags):
                    for event, repo, kind, ref in (
                        ('push', 'hallamlab/MetaPathways', 'tag', 'v3.5.2'),
                        ('workflow_dispatch', 'hallamlab/MetaPathways', 'branch', 'dev'),
                        ('workflow_dispatch', 'fork/MetaPathways', 'tag', 'v3.5.2'),
                        ('workflow_dispatch', 'hallamlab/MetaPathways', 'tag', 'invalid'),
                    ):
                        with self.assertRaises(ValueError):
                            validate(event, repo, kind, ref, **selections)
                else:
                    validate('push', 'fork/MetaPathways', 'branch', 'feature', **selections)

    def test_strings_cannot_enable_publication(self):
        with self.assertRaises(ValueError):
            validate('workflow_dispatch', 'hallamlab/MetaPathways', 'tag', 'v3.5.2', anaconda='false')
