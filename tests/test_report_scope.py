import csv
import io
import json
import sqlite3
import zipfile
import unittest
import test_reporting
from metapathways.reporting import build_report
from metapathways.report_server import query, write_mag_bundle
from metapathways.report_scope import BUNDLE_TABLES


class MAGScopeTests(unittest.TestCase):
    fixture = test_reporting.ReportTests.fixture
    write = test_reporting.ReportTests.write
    def setUp(self):
        test_reporting.ReportTests.setUp(self)
        self.write('alpha','genbank/alpha.annot.gff', 'C1\tx\ttRNA\t20\t50\t.\t+\t.\tID=T1;isotype=Ala;anticodon=TGC\n')
        build_report(self.root)
        self.db.close()
        self.db=sqlite3.connect(self.reports/'results.sqlite')
        self.db.row_factory=sqlite3.Row
        self.addCleanup(self.db.close)
        self.scope={'sample_id':'alpha','entity_id':'MAG1'}

    def test_scoped_tables_and_no_cross_sample_or_community_rows(self):
        for table in BUNDLE_TABLES:
            result=query(self.db,{'table':table,'scope':self.scope},allow_sequence_scan=True)
            for row in result['rows']:
                if 'sample_id' in row: self.assertEqual(row['sample_id'],'alpha')
                if 'entity_id' in row: self.assertEqual(row['entity_id'],'MAG1')
        # Both mapped genes belong, even though only G1 is a Pathway Tools input.
        self.assertEqual(query(self.db,{'table':'orf_explorer','scope':self.scope})['total'],2)
        self.assertEqual(query(self.db,{'table':'feature_explorer','scope':self.scope})['total'],3)
        self.assertEqual(query(self.db,{'table':'pathway_gene_explorer','scope':self.scope})['total'],1)
        self.assertEqual(query(self.db,{'table':'annotation_explorer','scope':self.scope})['total'],2)
        # An unrelated contig must not leak into the selection.
        self.db.execute("INSERT INTO contigs VALUES('alpha','C2','original2',100)")
        self.assertEqual(query(self.db,{'table':'contigs','scope':self.scope})['total'],1)

    def test_scope_errors_and_metadata_not_silently_unfiltered(self):
        for scope in ({'entity_id':'MAG1'}, {'sample_id':'beta','entity_id':'unknown'},
                      {'sample_id':'alpha','entity_id':'community'},
                      {'sample_id':'alpha','entity_id':'MAG_failed'}):
            with self.assertRaises(ValueError): query(self.db,{'table':'contigs','scope':scope})
        with self.assertRaises(ValueError): query(self.db,{'table':'files','scope':self.scope})

    def test_zip_has_all_rows_separate_tables_and_provenance(self):
        output=io.BytesIO()
        write_mag_bundle(self.db,self.scope,output)
        with zipfile.ZipFile(output) as archive:
            manifest=json.loads(archive.read('manifest.json'))
            self.assertEqual(manifest['scope'],self.scope)
            self.assertEqual(manifest['rows']['feature_explorer'],3)
            self.assertEqual(set(manifest['rows']),set(BUNDLE_TABLES))
            for table in BUNDLE_TABLES:
                data=list(csv.DictReader(io.StringIO(archive.read(table+'.csv').decode())))
                self.assertEqual(len(data),manifest['rows'][table])
                for row in data:
                    if 'sample_id' in row: self.assertEqual(row['sample_id'],'alpha')
            self.assertIn('sources.csv',archive.namelist())
            self.assertIn('README.txt',archive.namelist())

    def test_shared_contig_does_not_duplicate_scoped_results(self):
        self.db.execute("INSERT INTO entities VALUES('alpha','MAG2','MAG','unavailable',NULL)")
        self.db.execute("INSERT INTO contig_mags VALUES('alpha','C1','MAG2','MAG2',NULL)")
        self.assertEqual(query(self.db,{'table':'annotation_explorer','scope':self.scope})['total'],2)
        self.assertEqual(query(self.db,{'table':'contigs','scope':{'sample_id':'alpha','entity_id':'MAG2'}})['total'],1)
