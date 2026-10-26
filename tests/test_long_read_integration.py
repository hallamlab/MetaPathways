"""Opt-in known-count mapping test using installed CoverM/minimap2/Subread tools."""
import csv
import os
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from metapathways import MetaPathways_tpm as tpm


@unittest.skipUnless(os.environ.get('MP_TEST_BIO'), 'set MP_TEST_BIO with the MP environment activated')
class LongReadRuntimeTests(unittest.TestCase):
    def test_known_counts_for_hifi_pacbio_and_ont(self):
        for tool in ('coverm', 'minimap2', 'samtools', 'featureCounts', 'gff2gtf.py', 'abund_calc.py'):
            self.assertIsNotNone(shutil.which(tool), tool+' is required')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rng = random.Random(724)
            sequence = ''.join(rng.choice('ACGT') for _ in range(12000))
            (root/'contigs.fa').write_text('>contig1\n'+sequence+'\n')
            (root/'reads.fq').write_text(''.join('@read'+str(i)+'\n'+sequence[start:start+4000]+
                '\n+\n'+'I'*4000+'\n' for i,start in enumerate((1000,1500,2500))))
            (root/'annotation.gff').write_text('contig1\ttest\tCDS\t2001\t4000\t.\t+\t0\tID=gene1\n')
            for preset in ('map-hifi', 'map-pb', 'map-ont'):
                with self.subTest(preset=preset):
                    output = root/preset
                    output.mkdir()
                    args = ['-c',str(root/'contigs.fa'),'-g',str(root/'annotation.gff'),
                        '-1',str(root/'reads.fq'),'-o',str(output/'sample.contig_counts.tsv'),
                        '--stats',str(output/'mapping.log'),'--bwaFolder',str(output),
                        '--sample_name','sample','--rpkmExec','coverm','--num_threads','2',
                        '--long_read_preset',preset]
                    self.assertEqual(tpm.main(args), 0)
                    with (output/'sample.annot.featurecounts.txt.summary').open() as stream:
                        summary = list(csv.reader(stream, delimiter='\t'))
                    self.assertEqual(dict(summary[1:])['Assigned'], '3')
                    with (output/'sample.orf_counts.tsv').open() as stream:
                        records = list(csv.DictReader(stream, delimiter='\t'))
                    self.assertEqual(len(records), 1)
                    self.assertEqual(records[0]['Gene_ID'], 'gene1')
                    self.assertEqual(float(records[0]['Count']), 3)
                    self.assertEqual(float(records[0]['TPM']), 1000000)
