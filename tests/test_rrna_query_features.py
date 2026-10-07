import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from metapathways.MetaPathways_rRNA_stats_calculator import process_blastout_file, append_taxonomic_information, main

class RrnaFeatureTests(unittest.TestCase):
    def test_coordinate_and_legacy_headers_resolve_subunit(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            queries = ['ctg:0-250(+)', 'ctg:299-599(-)', '16S_rRNA::old:0-250(+)', 'old23_23S_rRNA']
            (d/'query.fna').write_text(''.join('>'+q+'\n'+'A'*300+'\n' for q in queries))
            (d/'query.gff').write_text('ctg\tbarrnap\trRNA\t1\t250\t.\t+\t.\tName=16S_rRNA\nctg\tbarrnap\trRNA\t300\t599\t.\t-\t.\tName=23S_rRNA\n')
            (d/'blast.tsv').write_text(''.join(q+'\tref\t99\t240\t0\t0\t1\t240\t1\t240\t1e-40\t200\n' for q in queries))
            for subunit, expected in [('16S',queries[:1]+queries[2:3]), ('23S',[queries[1],queries[3]])]:
                table = {}
                process_blastout_file(d/'blast.tsv','test',table,subunit,d/'query.fna',query_gff=d/'query.gff')
                self.assertEqual(set(table),set(expected))

    def test_inclusive_alignment_length_and_reverse_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            db = Path(directory) / 'silva.fasta'
            db.write_text('>ref Bacteria;Test species\nAAAA\n')
            table = {'forward': [99, 1e-40, 200, 'ref', 1, 180],
                     'reverse': [99, 1e-40, 200, 'ref', 180, 1],
                     'short': [99, 1e-40, 200, 'ref', 1, 179]}
            append_taxonomic_information(str(db), table, dict(length=180, similarity=20, evalue=1e-6, bitscore=50))
            self.assertEqual(table['forward'][6], 'Bacteria;Test species')
            self.assertEqual(table['reverse'][6], 'Bacteria;Test species')
            self.assertEqual(table['short'][6], '-')

    def test_cli_writes_taxonomy_for_coordinate_only_query(self):
        with tempfile.TemporaryDirectory() as directory:
            d = Path(directory)
            (d/'query.fna').write_text('>ctg:0-250(+)\n'+'A'*250+'\n')
            (d/'query.gff').write_text('ctg\tbarrnap\trRNA\t1\t250\t.\t+\t.\tName=16S_rRNA\n')
            (d/'blast.tsv').write_text('ctg:0-250(+)\tref\t99\t240\t0\t0\t1\t240\t1\t240\t1e-40\t200\n')
            (d/'silva.fasta').write_text('>ref Bacteria;Test species\n'+'A'*250+'\n')
            args=['-q',str(d/'query.fna'),'--query-gff',str(d/'query.gff'),'-r','16S','-i',str(d/'blast.tsv'),'-d',str(d/'silva.fasta'),'-o',str(d/'stats.tsv')]
            with patch('metapathways.MetaPathways_rRNA_stats_calculator.runBlastCommandrRNA', return_value=0):
                main(args)
            text=(d/'stats.tsv').read_text()
            self.assertIn('#Number of rRNA sequences detected:\t1', text)
            self.assertIn('Bacteria;Test species', text)

if __name__ == '__main__': unittest.main()
