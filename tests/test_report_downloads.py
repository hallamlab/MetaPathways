import csv
import gzip
import io
import json
import sqlite3
import tempfile
import unittest
import zipfile
from pathlib import Path
import test_reporting
from metapathways.reporting import build_report, VIEWS
from metapathways.report_server import query
from metapathways.report_downloads import write_download, preview, selections
from metapathways.report_tables import TABLE_GROUPS
from metapathways.report_sequences import records


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        fixture=test_reporting.ReportTests();fixture.root=self.root
        for sample in ('alpha','beta'):
            fixture.fixture(sample)
            fixture.write(sample,'preprocessed/test.mapping.txt','C1\toriginal1\t300\nC2\toriginal2\t100\n')
            fixture.write(sample,'preprocessed/test.fasta','>C1 description\n'+'ACGT'*75+'\n>C2\nAAAA\n')
            fixture.write(sample,'orf_prediction/test.faa','>G1 protein\nMALW\n>G2\nMFG\n')
            fixture.write(sample,'orf_prediction/test.fna','>G1\nATGGCC\n>G2\nATGTTT\n')
            fixture.write(sample,'reads/test.fastq','@C1\nACGT\n+\n!!!!\n')
            fixture.write(sample,'genbank/test.annot.gff','C1\tx\trRNA\t11\t20\t.\t+\t.\tID=R1;product=16S ribosomal RNA\nC2\tx\trRNA\t1\t4\t.\t+\t.\tID=R2;product=23S ribosomal RNA\n')
            fixture.write(sample,'orf_prediction/test.rRNA.fna','>16S_rRNA::C1:10-20(+)\nACGTACGTAC\n')
            fixture.write(sample,'genbank/test.gbk','LOCUS       C1 300 bp DNA\nFEATURES             Location/Qualifiers\n     CDS             1..90\n                     /product="kinase"\nORIGIN\n        1 acgt\n//\nLOCUS       C2 100 bp DNA\nORIGIN\n        1 aaaa\n//\n')
            fixture.write(sample,'bwa/sample.bam','shared sample BAM')
        self.reports=build_report(self.root)
        self.db=sqlite3.connect(self.reports/'results.sqlite');self.db.row_factory=sqlite3.Row;self.addCleanup(self.db.close)

    def bundle(self,selections,mode='complete'):
        out=io.BytesIO();write_download(self.db,self.root,{'mode':mode,'selections':selections},out)
        return zipfile.ZipFile(out)

    def test_mixed_samples_mags_deduplicate_without_leaking_contigs(self):
        with self.bundle([{'sample_id':'alpha'}, {'sample_id':'alpha','entity_id':'MAG1'}, {'sample_id':'beta','entity_id':'MAG1'}]) as archive:
            manifest=json.loads(archive.read('manifest.json'))
            self.assertEqual(len(manifest['selections']),2)
            names=archive.namelist()
            self.assertIn('sample-alpha/complete-sample/files/alpha/bwa/sample.bam',names)
            self.assertFalse(any('/files/beta/bwa/' in n for n in names))
            prefix='sample-beta/mag-MAG1'
            gbk=archive.read(prefix+'/annotations/test.gbk').decode()
            self.assertIn('/product="kinase"',gbk);self.assertNotIn('LOCUS       C2',gbk)
            gff=archive.read(prefix+'/annotations/test.annot.gff').decode()
            self.assertIn('ID=R1',gff);self.assertNotIn('ID=R2',gff)
            contigs=archive.read(prefix+'/sequences/contig.fasta').decode()
            self.assertIn('C1 description',contigs);self.assertNotIn('C2',contigs)
            self.assertIn('MALW',archive.read(prefix+'/sequences/protein.fasta').decode())
            self.assertIn('ATGGCC',archive.read(prefix+'/sequences/nucleotide.fasta').decode())
            self.assertIn('ACGTACGTAC',archive.read(prefix+'/sequences/RNA.fasta').decode())
            rows=list(csv.DictReader(io.StringIO(archive.read(prefix+'/tables/annotation_explorer.csv').decode())))
            self.assertEqual(len(rows),2);self.assertEqual({r['sample_id'] for r in rows},{'beta'})
            seqs=list(csv.DictReader(io.StringIO(archive.read(prefix+'/tables/sequence_explorer.csv').decode())))
            self.assertNotIn('read',{r['sequence_type'] for r in seqs})
            self.assertTrue(all(r['contig_id']=='C1' for r in seqs))

    def test_table_only_and_all_samples(self):
        chosen=selections(self.db,{'mode':'tables','all_samples':True})
        self.assertEqual(chosen,[{'sample_id':'alpha'},{'sample_id':'beta'}])
        with self.bundle(chosen,'tables') as archive:
            self.assertFalse(any('/files/' in p or '/annotations/' in p for p in archive.namelist()))
            self.assertIn('sample-alpha/complete-sample/tables/sequence_explorer.csv',archive.namelist())
        self.assertEqual(len(preview(self.db,{'mode':'complete','all_samples':True})['items']),2)

    def test_review_counts_match_zip_members_in_each_mode(self):
        choices=[{'sample_id':'alpha'}, {'sample_id':'alpha','entity_id':'MAG1'},
                 {'sample_id':'beta','entity_id':'MAG1'}]
        for mode in ('tables','complete'):
            spec={'mode':mode,'selections':choices}
            review=preview(self.db,spec)
            with self.bundle(choices,mode) as archive:
                names=archive.namelist()
                self.assertEqual(review['total_files'],len(names))
                manifest=json.loads(archive.read('manifest.json'))
                self.assertEqual(len(review['items']),2)
                for item,entry in zip(review['items'],manifest['selections']):
                    prefix=entry['directory']+'/'
                    files=[n for n in names if n.startswith(prefix)]
                    self.assertEqual(item['total_files'],len(files))
                    self.assertEqual(item['file_counts']['table_csv'],sum('/tables/' in n for n in files))
                    self.assertEqual(item['file_counts']['sequence_fasta'],sum('/sequences/' in n for n in files))
                    self.assertEqual(item['file_counts']['annotation_subset'],sum('/annotations/' in n for n in files))
                    copied=[n for n in files if n.startswith(prefix+'files/')]
                    self.assertEqual(item['retained_files'],len(copied))
                    self.assertEqual(item['retained_bytes'],sum(archive.getinfo(n).file_size for n in copied))
                    if mode=='tables':
                        self.assertEqual(item['retained_files'],0)
                        self.assertEqual(item['retained_bytes'],0)
        # Empty FASTA subsets are still ZIP members and must be counted.
        self.db.execute("DELETE FROM sequences WHERE sample_id='beta'")
        self.db.execute("DELETE FROM files WHERE path LIKE 'beta/genbank/%'")
        spec={'mode':'complete','selections':[{'sample_id':'beta','entity_id':'MAG1'}]}
        with self.bundle(spec['selections']) as archive:
            self.assertEqual(preview(self.db,spec)['total_files'],len(archive.namelist()))

    def test_sequence_inspection_is_one_record_and_bulk_tables_stay_private(self):
        from metapathways.report_fasta import options,record
        for spec in ({'table':'sequence_explorer'}, {'table':'sequence_explorer','scope':{'sample_id':'alpha'}}):
            with self.assertRaisesRegex(ValueError,'Unknown table'):query(self.db,spec)
        data=options(self.db,{'sample_id':'alpha','contig_id':'C1'})
        self.assertEqual(len(data['records']),1)
        selected={'sample_id':'alpha','sequence_record_id':data['records'][0]['sequence_record_id']}
        self.assertEqual(record(self.db,selected)['length'],300)
        self.assertFalse(record(self.db,selected)['truncated'])
        self.db.execute('UPDATE sequences SET sequence=?,length=? WHERE sequence_record_id=?',('A'*20000,20000,selected['sequence_record_id']))
        self.assertTrue(record(self.db,selected)['truncated'])
        self.assertEqual(record(self.db,selected,full=True)['fasta'].count('A'),20000)
        feature=options(self.db,{'sample_id':'alpha','feature_id':'G1'})
        self.assertEqual({r['sequence_type'] for r in feature['records']},{'protein','nucleotide'})
        with self.assertRaises(ValueError):options(self.db,{'sample_id':'alpha'})
        with self.assertRaises(ValueError):record(self.db,{**selected,'sample_id':'beta'})
        self.assertIsNone(self.db.execute("SELECT contig_id FROM sequences WHERE sample_id='alpha' AND sequence_type='read'").fetchone()[0])

    def test_changed_files_and_symlinks_are_rejected(self):
        path=self.root/'alpha/bwa/sample.bam';path.write_text('changed')
        with self.assertRaisesRegex(ValueError,'changed since'): self.bundle([{'sample_id':'alpha'}])
        path.unlink();path.symlink_to('/etc/passwd')
        with self.assertRaises(ValueError): self.bundle([{'sample_id':'alpha'}])

    def test_missing_genbank_is_reported_and_invalid_scope_rejected(self):
        self.db.execute("DELETE FROM files WHERE path LIKE 'alpha/genbank/%.gbk'")
        with self.bundle([{'sample_id':'alpha','entity_id':'MAG1'}]) as archive:
            item=json.loads(archive.read('manifest.json'))['selections'][0]
            self.assertIn('No retained genbank_records available for this selection',item['availability_notes'])
        for items in ([],[{'sample_id':'missing'}],[{'entity_id':'MAG1'}]):
            with self.assertRaises(ValueError): selections(self.db,{'mode':'complete','selections':items})
        flattened=[t for group in TABLE_GROUPS.values() for t in group]
        self.assertEqual(set(flattened),set(VIEWS)-{'abundance'})
        self.assertEqual(len(flattened),len(set(flattened)))

    def test_gzip_and_wrapped_fastq(self):
        path=self.root/'wrapped.fastq.gz'
        with gzip.open(path,'wt') as f:f.write('@read1\nAC\nGT\n+\n!!\n##\n')
        self.assertEqual(list(records(path)),[('read1','ACGT','!!##')])
