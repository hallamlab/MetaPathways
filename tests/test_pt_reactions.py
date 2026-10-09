import json
from pathlib import Path
import tempfile
import unittest

from metapathways.pt_reactions import filter_reactions


class ReactionFilterTests(unittest.TestCase):
    def test_only_blocked_assignments_removed_and_audited(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = ('ID\tC137456-G1\nNAME\tC137456-G1\n'
                      'FUNCTION\tsublancin 168 maturation protease / ABC efflux transporter\n'
                      'METACYC\tTRANS-RXN8J2-121\nMETACYC\tSAFE-RXN\n'
                      'EC\t3.4.21.1\nPRODUCT-TYPE\tP\nSTARTBASE\t3\nENDBASE\t242\n//\n')
            original = root/'original.txt'
            original.write_text(source)
            for name in ('0.pf', 'contig_1.pf'):
                (root/name).write_text(source)
            audit = filter_reactions(root)
            for name in ('0.pf', 'contig_1.pf'):
                self.assertEqual((root/name).read_text(), source.replace('METACYC\tTRANS-RXN8J2-121\n', ''))
            self.assertEqual(original.read_text(), source)
            self.assertEqual(len(audit['removed']), 2)
            self.assertTrue(all(r['feature_id'] == 'C137456-G1' for r in audit['removed']))
            self.assertEqual(json.loads((root/'ptools-reaction-filter.json').read_text()), audit)

    def test_unlisted_reaction_ids_are_not_substring_matches(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = 'ID\tgene\nMETACYC\tTRANS-RXN8J2-1210\nMETACYC\tSAFE-RXN\n//\n'
            (root/'0.pf').write_text(source)
            self.assertEqual(filter_reactions(root)['removed'], [])
            self.assertEqual((root/'0.pf').read_text(), source)

    def test_matching_db_list_and_stale_mapping_fallback(self):
        from metapathways.pt_screen import digest
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            categories = root/'db/functional_categories'
            categories.mkdir(parents=True)
            mapping = categories/'MetaCyc-monomer-rxn-pairs.tsv'
            mapping.write_text('MC\tRXN\nM\tCUSTOM\n')
            image = root/'pt.sif'
            image.write_bytes(b'test image')
            inputs = root/'inputs'
            inputs.mkdir()
            (root/'metapathways_run_log.txt').write_text(f'Minimum Required Arguments:refdb_dir\t{root/"db"}\n')
            compatibility = categories/'ptools_reaction_compatibility.json'
            compatibility.write_text(json.dumps(dict(image_sha256=digest(image), mapping_sha256=digest(mapping),
                                                     reactions={'CUSTOM': {'reason': 'verified'}})))
            pf = inputs/'0.pf'
            pf.write_text('ID\tG1\nMETACYC\tCUSTOM\n//\n')
            audit = filter_reactions(inputs, image=image, sample_output=root)
            self.assertEqual(audit['blacklist_source'], str(compatibility))
            self.assertEqual(len(audit['removed']), 1)
            mapping.write_text('changed')
            pf.write_text('ID\tG1\nMETACYC\tCUSTOM\n//\n')
            self.assertEqual(filter_reactions(inputs, image=image, sample_output=root)['removed'], [])
            self.assertIn('METACYC\tCUSTOM', pf.read_text())


if __name__ == '__main__':
    unittest.main()
