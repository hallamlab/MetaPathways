import sqlite3
import tempfile
import unittest
from pathlib import Path
from metapathways.reporting import build_report


class RNAReportTests(unittest.TestCase):
    def test_rna_links_preserve_samples_loci_and_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for sample in ('alpha', 'beta'):
                def write(path, text):
                    target = root/sample/path
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_text(text)
                write('preprocessed/test.mapping.txt', 'C1\toriginal\t1000\n')
                write('genbank/test.annot.gff',
                      'C1\tbarrnap\trRNA\t11\t100\t.\t+\t.\tID=R1;product=16S ribosomal RNA\n'
                      'C1\tbarrnap\trRNA\t201\t300\t.\t-\t.\tID=R2;product=23S ribosomal RNA\n'
                      'C1\ttRNAscan\ttRNA\t401\t480\t50\t-\t.\tID=T1;isotype=Ala;anticodon=TGC\n')
                write('results/rRNA/test.SSU.rRNA.stats.txt',
                      'preamble\nsequence\tstart\tend\tsimilarity\tevalue\tbitscore\ttaxonomy\n'
                      '16S_rRNA::C1:10-100(+)\t2\t80\t99\t1e-20\t100\tSpecies A\n'
                      'missing\t1\t2\t95\t1e-10\t40\tUnlinked evidence\n')
                write('results/rRNA/test.LSU.rRNA.stats.txt',
                      'sequence\tstart\tend\tsimilarity\tevalue\tbitscore\ttaxonomy\n'
                      'C1:200-300(-)\t1\t90\t97\t1e-15\t90\tSpecies B\n')
                write('magsplitter/contig_to_mag.tsv', 'original\tMAG.1\n')
                write('results/rpkm/test.orf_counts.tsv',
                      'Gene_ID\tCount\tRPKM\tTPM\tfeature\tseqname\n'
                      'G1\t10\t1\t2\tCDS\tC1\nT1\t3\t4\t5\ttRNA\tC1\n')
            report = build_report(root)
            with sqlite3.connect(report/'results.sqlite') as db:
                self.assertEqual(db.execute('PRAGMA foreign_key_check').fetchall(), [])
                self.assertEqual(db.execute('SELECT count(*) FROM rna_features').fetchone()[0], 6)
                self.assertEqual(db.execute('SELECT count(*) FROM mag_rna_explorer').fetchone()[0], 6)
                self.assertEqual(db.execute('SELECT count(*) FROM orf_rna_explorer').fetchone()[0], 6)
                self.assertEqual(db.execute('SELECT count(*) FROM rna_abundance_explorer').fetchone()[0], 2)
                self.assertEqual(db.execute("SELECT isotype,anticodon FROM rna_features WHERE sample_id='alpha' AND rna_id='T1'").fetchone(), ('Ala','TGC'))
                self.assertEqual(db.execute("SELECT rna_id,start,hit_start,taxonomy FROM rna_taxonomy_explorer WHERE sample_id='alpha' ORDER BY evidence_id").fetchall(), [('R2',201,1,'Species B'), ('R1',11,2,'Species A'), (None,None,1,'Unlinked evidence')])
                self.assertEqual(db.execute('SELECT count(*) FROM issues WHERE message LIKE ?',( '%no unique RNA locus%',)).fetchone()[0], 2)
                self.assertEqual(db.execute("SELECT count(*) FROM rna_taxonomy WHERE rna_id='T1'").fetchone()[0],0)

    def test_missing_taxonomy_and_raw_gff_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'preprocessed').mkdir()
            (root/'orf_prediction').mkdir()
            (root/'orf_prediction/test.rRNA.gff').write_text(
                'C1\tbarrnap\trRNA\t11\t100\t.\t+\t.\tName=16S_rRNA\n')
            report = build_report(root)
            with sqlite3.connect(report/'results.sqlite') as db:
                self.assertEqual(db.execute('SELECT subunit FROM rna_features').fetchall(), [('16S',)])
                self.assertEqual(db.execute('SELECT count(*) FROM rna_taxonomy').fetchone()[0],0)
