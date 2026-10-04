from pathlib import Path
import tempfile
import unittest

from metapathways.pt_sequences import attach_sequences


class SequenceInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.base = self.root/'sample'
        self.inputs = self.root/'private-input'
        self.inputs.mkdir()
        (self.base/'preprocessed').mkdir(parents=True)
        (self.base/'preprocessed/sample.fasta').write_text('>contigA\nATGAAATAG\n>contigB\nTTATTTCAT\n')
        table = self.base/'results/annotation_table/sample.ptinput.tsv'
        table.parent.mkdir(parents=True)
        table.write_text('orf_id\tseqname\tstart\tend\tstrand\nA-G1\tcontigA\t1\t9\t+\nB-G1\tcontigB\t1\t9\t-\n')
        self.records = ['ID\tA-G1\nSTARTBASE\t1\nENDBASE\t9\nFUNCTION\tx\nPRODUCT-TYPE\tP\n//\n',
                        'ID\tB-G1\nSTARTBASE\t9\nENDBASE\t1\nFUNCTION\ty\nPRODUCT-TYPE\tP\n//\n']
        (self.inputs/'0.pf').write_text(''.join(self.records))
        (self.inputs/'genetic-elements.dat').write_text('old compact manifest')

    def test_real_contigs_keep_feature_ids_coordinates_and_strands(self):
        attach_sequences(self.inputs, self.base)
        manifest = (self.inputs/'genetic-elements.dat').read_text()
        self.assertEqual(manifest.count('SEQ-FILE\t'), 2)
        self.assertEqual(set((self.inputs/'contig_1.pf').read_text().splitlines()), set(self.records[0].splitlines()))
        self.assertEqual(set((self.inputs/'contig_2.pf').read_text().splitlines()), set(self.records[1].splitlines()))
        self.assertEqual((self.inputs/'contig_1.fasta').read_text(), '>contig_1\nATGAAATAG\n')
        self.assertEqual((self.inputs/'contig_2.fasta').read_text(), '>contig_2\nTTATTTCAT\n')
        self.assertEqual((self.inputs/'0.pf').read_text(), ''.join(self.records))

    def test_taxon_override_replaces_only_taxon_and_rejects_invalid_ids(self):
        from metapathways.pt_sequences import set_organism_taxon
        params = self.inputs/'organism-params.dat'
        params.write_text('ID\tsample\nNCBI-TAXON-ID\t12908\nSTORAGE\tFILE\n')
        set_organism_taxon(self.inputs, 131567)
        self.assertEqual(params.read_text(), 'ID\tsample\nSTORAGE\tFILE\nNCBI-TAXON-ID\t131567\n')
        for value in (0, -1, '131567', True):
            with self.assertRaises(ValueError):
                set_organism_taxon(self.inputs, value)

    def test_trna_name_formatting_preserves_ids_and_unknown_anticodons(self):
        import json
        from metapathways.pt_sequences import normalize_trna_name
        for label in ('GluTTC', 'fMetCAT', 'LeuTAA', 'SeCTCA'):
            name = 'C1094.tRNA2-' + label
            expected = name[:-3] + '-' + name[-3:]
            self.assertEqual(normalize_trna_name(name), expected)
            self.assertEqual(normalize_trna_name(expected), expected)
        for name in ('C1.tRNA1-UndetNNN', 'C1.tRNA1-GluNNC', 'customGluTTC'):
            self.assertEqual(normalize_trna_name(name), name)
        text = self.records[0].replace('FUNCTION\tx', 'NAME\tC1094.tRNA2-GluTTC').replace('PRODUCT-TYPE\tP', 'PRODUCT-TYPE\tTRNA')
        (self.inputs/'0.pf').write_text(text)
        attach_sequences(self.inputs, self.base)
        written = (self.inputs/'contig_1.pf').read_text()
        self.assertIn('ID\tA-G1\n', written)
        self.assertIn('NAME\tC1094.tRNA2-Glu-TTC\n', written)
        self.assertEqual((self.inputs/'0.pf').read_text(), text)
        audit = json.loads((self.inputs/'sequence-input.json').read_text())
        self.assertEqual(audit['normalized_trna_names']['A-G1'],
                         {'original': 'C1094.tRNA2-GluTTC', 'staged': 'C1094.tRNA2-Glu-TTC'})

    def test_prodigal_genetic_codes_are_preserved_per_contig(self):
        folder = self.base/'orf_prediction'
        folder.mkdir()
        gff = folder/'sample.cds.gff'
        gff.write_text('# Sequence Data: seqhdr="contigA"\n# Model Data: transl_table=4;uses_sd=0\n'
                       '# Sequence Data: seqhdr="contigB"\n# Model Data: transl_table=11;uses_sd=1\n')
        attach_sequences(self.inputs, self.base)
        manifest = (self.inputs/'genetic-elements.dat').read_text()
        self.assertIn('NAME\tcontigA\nTYPE\t:CONTIG\nCODON-TABLE\t4', manifest)
        self.assertIn('NAME\tcontigB\nTYPE\t:CONTIG\nCODON-TABLE\t11', manifest)
        gff.write_text(gff.read_text() + '# Sequence Data: seqhdr="contigA"\n# Model Data: transl_table=11;\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            attach_sequences(self.inputs, self.base)

    def test_new_pf_writer_normalizes_string_and_list_ec_values(self):
        from io import StringIO
        from metapathways.MetaPathways_create_genbank_ptinput import write_to_pf_file
        for value in ('1.2.3.4,2.3.4.5', ['1.2.3.4,2.3.4.5', '1.2.3.4']):
            handle = StringIO()
            write_to_pf_file(str(self.inputs), 'A-G1',
                             dict(strand='+', start=1, end=9, feature='CDS', ec=value), handle, True)
            self.assertEqual([line for line in handle.getvalue().splitlines() if line.startswith('EC\t')],
                             ['EC\t1.2.3.4', 'EC\t2.3.4.5'])

    def test_legacy_ec_lists_are_split_without_losing_provisional_values(self):
        import json
        text = self.records[0].replace('//', 'EC\t1.2.3.4, 2.3.4.5\nEC\t1.2.3.4|3.6.5.n1\n//')
        (self.inputs/'0.pf').write_text(text)
        attach_sequences(self.inputs, self.base)
        written = (self.inputs/'contig_1.pf').read_text()
        self.assertEqual([line for line in written.splitlines() if line.startswith('EC\t')],
                         ['EC\t1.2.3.4', 'EC\t2.3.4.5', 'EC\t3.6.5.n1'])
        self.assertEqual((self.inputs/'0.pf').read_text(), text)
        self.assertEqual(json.loads((self.inputs/'sequence-input.json').read_text())['normalized_ec_records'], 1)

    def test_mag_gets_only_the_contigs_its_records_reference(self):
        (self.inputs/'0.pf').write_text(self.records[1])
        attach_sequences(self.inputs, self.base)
        self.assertEqual((self.inputs/'genetic-elements.dat').read_text().count('SEQ-FILE\t'), 1)
        self.assertIn('NAME\tcontigB', (self.inputs/'genetic-elements.dat').read_text())
        self.assertIn('TTATTTCAT', (self.inputs/'contig_1.fasta').read_text())

    def test_missing_sequence_fails_before_replacing_manifest(self):
        (self.base/'preprocessed/sample.fasta').write_text('>contigA\nATGAAATAG\n')
        with self.assertRaisesRegex(ValueError, 'Missing Pathway Tools contig sequences: contigB'):
            attach_sequences(self.inputs, self.base)
        self.assertEqual((self.inputs/'genetic-elements.dat').read_text(), 'old compact manifest')

    def test_invalid_coordinates_fail_before_replacing_manifest(self):
        table = self.base/'results/annotation_table/sample.ptinput.tsv'
        table.write_text(table.read_text().replace('contigA\t1\t9', 'contigA\t1\t10'))
        with self.assertRaisesRegex(ValueError, 'outside contigA'):
            attach_sequences(self.inputs, self.base)
        self.assertEqual((self.inputs/'genetic-elements.dat').read_text(), 'old compact manifest')

    def test_mag_representative_coordinates_are_restored_to_actual_member(self):
        (self.inputs/'0.pf').write_text(self.records[1].replace('STARTBASE\t9', 'STARTBASE\t900'))
        attach_sequences(self.inputs, self.base)
        restored = (self.inputs/'contig_1.pf').read_text()
        self.assertIn('STARTBASE\t9\n', restored)
        self.assertIn('ENDBASE\t1\n', restored)
        self.assertIn('"restored_coordinates": 1', (self.inputs/'sequence-input.json').read_text())

    def test_unannotated_contig_is_not_required_or_exported(self):
        with (self.base/'preprocessed/sample.fasta').open('a') as stream:
            stream.write('>unannotated\nATGC\n')
        attach_sequences(self.inputs, self.base)
        self.assertNotIn('unannotated', (self.inputs/'genetic-elements.dat').read_text())

    def test_unmapped_or_duplicate_features_fail(self):
        for records in [self.records[0].replace('A-G1', 'unknown'), self.records[0]*2]:
            (self.inputs/'0.pf').write_text(records)
            with self.assertRaisesRegex(ValueError, 'Duplicate or unmapped'):
                attach_sequences(self.inputs, self.base)


if __name__ == '__main__':
    unittest.main()
