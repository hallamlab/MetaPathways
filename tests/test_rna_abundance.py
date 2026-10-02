"""RNA locus identity, annotation buffer boundaries and abundance regressions."""
import importlib.util
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


def helper(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parents[1] / 'dev' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RnaAbundanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def put(self, name, text):
        p = self.root / name
        p.write_text(text)
        return p

    def test_rna_emitted_once_across_buffers_and_on_rna_only_contigs(self):
        from metapathways import MetaPathways_annotate_fast as annot
        cds = self.put('cds.gff', ''.join(f's-C1\tprodigal\tCDS\t{i}\t{i+2}\t1\t+\t0\tID=s-C1-G{i}\n' for i in (1, 4)))
        rrna = self.put('rrna.gff', ''.join(f'{contig}\tbarrnap\trRNA\t{start}\t{start+9}\t1\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n' for contig, start in [('s-C1', 10), ('s-C1', 30), ('s-C2', 10)]))
        trna = self.put('trna.gff', 's-C1\ttRNAscan\ttRNA\t50\t59\t1\t+\t.\tID=s-C1.trna1;Name=s-C1.tRNA1-Ala;isotype=Ala;anticodon=CGC\n')
        original = annot.GffFileParser.__init__
        def small_buffer(obj, *args, **kwargs):
            original(obj, *args, **kwargs)
            obj.Size = 1
        output = self.root / 'annot.gff'
        with patch.object(annot.GffFileParser, '__init__', small_buffer):
            annot.create_annotation({}, {}, str(cds), ['stats'], str(rrna), ['stats'], str(trna),
                str(output), str(self.root/'comparison'), {'s-C1': 100, 's-C2': 100}, [],
                '>s-C1\nACGT\n>s-C2\nACGT\n', 's')
        rows = [line.split('\t') for line in output.read_text().splitlines() if '\t' in line]
        self.assertEqual(len(rows), 4)
        ids = [row[8].split(';')[0] for row in rows]
        self.assertEqual(len(set(ids)), 4)
        self.assertEqual(sum(row[0] == 's-C2' for row in rows), 1)

    def fixtures(self, duplicate=False, zero=False):
        counts = self.put('counts.tsv', 'Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\nG1\tC1\t1\t100\t+\t100\t' + ('0' if zero else '10') + '\nG2\tC1\t201\t400\t+\t200\t' + ('0' if zero else '10') + '\n')
        gtf = 'C1\tx\tCDS\t1\t100\t.\t+\t0\tgene_id "G1";\nC1\tx\trRNA\t201\t400\t.\t+\t.\tgene_id "G2";\n'
        if duplicate:
            gtf += 'C1\tx\trRNA\t501\t700\t.\t+\t.\tgene_id "G2";\n'
        return counts, self.put('input.gtf', gtf)

    def test_abundance_uses_featurecounts_lengths_and_one_to_one_join(self):
        abund = helper('abund_calc')
        counts, gtf = self.fixtures()
        legacy = self.put('legacy.tsv', 'G1\t999\nG2\t999\n')
        result = abund.abundance_table(counts, gtf, legacy)
        self.assertEqual(result.Length.tolist(), [100, 200])
        self.assertEqual(result.Count.tolist(), [10, 10])
        self.assertAlmostEqual(result.TPM.sum(), 1e6)
        self.assertAlmostEqual(result.TPM.iloc[0] / result.TPM.iloc[1], 2)
        self.assertAlmostEqual(result.RPKM.iloc[0], 5000000.)
        self.assertAlmostEqual(result.RPKM.iloc[1], 2500000.)
        self.assertFalse(result.Gene_ID.duplicated().any())

    def test_duplicate_loci_rejected_instead_of_silently_pooled(self):
        counts, gtf = self.fixtures(duplicate=True)
        with self.assertRaisesRegex(ValueError, 'duplicate gene IDs'):
            helper('abund_calc').abundance_table(counts, gtf)

    def test_zero_counts_are_zero_not_nan(self):
        counts, gtf = self.fixtures(zero=True)
        result = helper('abund_calc').abundance_table(counts, gtf)
        self.assertEqual(result.TPM.tolist(), [0, 0])
        self.assertEqual(result.RPKM.tolist(), [0, 0])

    def test_gff_conversion_rejects_duplicate_ids_before_mapping(self):
        gff = self.put('bad.gff', 'C1\tx\trRNA\t1\t100\t.\t+\t.\tID=RNA1\nC1\tx\trRNA\t201\t300\t.\t+\t.\tID=RNA1\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate gene IDs'):
            helper('gff2gtf').gff_to_gtf(gff, self.root/'bad.gtf', ['rRNA'])

    def test_child_stderr_is_logged_and_failure_propagates(self):
        from metapathways.MetaPathways_tpm import run_logged
        log = io.StringIO()
        with self.assertRaises(subprocess.CalledProcessError) as error:
            run_logged([sys.executable, '-c', 'import sys; print("specific failure", file=sys.stderr); sys.exit(7)'], log)
        self.assertEqual(error.exception.returncode, 7)
        self.assertIn('specific failure', log.getvalue())


if __name__ == '__main__':
    unittest.main()
