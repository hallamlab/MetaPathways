import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch
from metapathways import MetaPathways_parse_blast as module

class BlastHitCountsTests(unittest.TestCase):
    def test_unique_counts_preserve_current_legacy_and_custom_identifiers(self):
        queries = ['C1-G1', 'C1-G1', 'C2-G1', 'sampleA_1_2', 'sampleB_1_2', 'custom-gene']
        fields = ['target', 'q_length', 'bitscore', 'bsr', 'expect', 'aln_length', 'identity', 'ec', 'product']
        records = [dict(query=q, **{f: 'value' for f in fields}) for q in queries]
        parser = MagicMock()
        parser.__iter__.return_value = iter(records)
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / 'parsed.tsv'
            opts = SimpleNamespace(ec_maps=None, taxonomy=False, parsed_output=str(output))
            with patch.object(module, 'BlastOutputParser', return_value=parser):
                hits, unique = module.process_blastoutput('test', 'unused', 'unused', 'unused', opts)
            self.assertEqual((hits, unique), (6, 5))
            self.assertEqual([line.split('\t')[0] for line in output.read_text().splitlines()[1:]], queries)

if __name__ == '__main__':
    unittest.main()
