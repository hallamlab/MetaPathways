import csv
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import threading
import unittest
from urllib.request import urlopen, Request
from urllib.error import HTTPError
from urllib.parse import urlencode

from metapathways.reporting import build_report
from metapathways.report_server import ReportServer, query, csv_value


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        for sample in ('alpha','beta'):
            self.fixture(sample)
        self.reports=build_report(self.root)
        self.db=sqlite3.connect(self.reports/'results.sqlite')
        self.db.row_factory=sqlite3.Row
        self.addCleanup(self.db.close)

    def test_taxonomy_is_joined_by_sample_orf_database_and_target(self):
        self.write('alpha', 'results/annotation_table/test.annotation_taxonomy.tsv',
            'orf_id\treference_db\ttarget\ttaxid\ttaxonomy\tlca_taxonomy\n'
            'G1\tswissprot\tACC1\t3\tSpecies A\tBacteria\n'
            'G1\tmetacyc\tWRONG_TARGET\t5\tArchaea\tArchaea\n')
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            actual = db.execute("SELECT reference_db,taxonomy,lca_taxonomy FROM annotation_explorer WHERE sample_id='alpha' ORDER BY reference_db").fetchall()
            self.assertEqual(actual, [('metacyc','Not computed','Not computed'), ('swissprot','Species A','Bacteria')])
            self.assertEqual(db.execute("SELECT DISTINCT taxonomy FROM annotation_explorer WHERE sample_id='beta'").fetchall(), [('Not computed',)])

    def test_wide_abundance_preserves_values_nulls_and_sample_identity(self):
        for sample in ('alpha','beta'):
            prefix = f'{sample}.fasta/{sample}.fastq.gz '
            headers = ['Length','Read Count','Mean','Variance','Trimmed Mean','RPKM','TPM','Extra']
            self.write(sample, 'results/rpkm/test.contig_counts.tsv',
                'Contig\t'+'\t'.join(prefix+h for h in headers)+'\nC1\t200\t0\t1.5\t2.5\t1.2\t3.4\t4.5\t9\n')
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            db.row_factory = sqlite3.Row
            self.assertEqual(db.execute('SELECT COUNT(*) FROM abundance_explorer').fetchone()[0],4)
            row = db.execute("SELECT * FROM abundance_explorer WHERE sample_id='alpha' AND feature_type='contig'").fetchone()
            self.assertEqual((row['length_bp'],row['count'],row['mean_coverage'],row['coverage_variance'],row['trimmed_mean_coverage'],row['rpkm'],row['tpm']), (200,0,1.5,2.5,1.2,3.4,4.5))
            orf = db.execute("SELECT * FROM abundance_explorer WHERE sample_id='alpha' AND feature_type='orf'").fetchone()
            self.assertEqual((orf['orf_id'],orf['contig_id'],orf['count']),('G1','C1',10))
            self.assertIsNone(orf['mean_coverage'])
            self.assertIsNone(orf['length_bp'])  # Absent in this legacy fixture; do not invent a length.
            self.assertEqual(db.execute("SELECT value FROM abundance WHERE sample_id='alpha' AND measurement LIKE '% Extra'").fetchone()[0],9)
            spec = {'table':'abundance_explorer','filters':[{'column':'count','op':'ge','value':10}], 'columns':['sample_id','feature_id','count']}
            columns, result = query(db, spec, export=True)
            self.assertEqual(columns, spec['columns'])
            self.assertEqual(len(list(result)),2)

    def test_report_footer_run_details_and_sample_default(self):
        self.write('alpha', 'logs/run/example/summary.json', json.dumps({
            'mp_version':'3.5.1', 'status':'SUCCESS', 'resources':{'executor':'slurm'}, 'tasks':[]}))
        build_report(self.root)
        meta = json.loads((self.reports/'schema.json').read_text())
        self.assertEqual(meta['run_details'], {'mp_version':'3.5.1','status':'SUCCESS','executor':'slurm','command':'run'})
        self.assertIn('report_mp_version', meta)
        page = (self.reports/'EDA_portal.html').read_text()
        self.assertIn('issues/new/choose', page)
        self.assertNotIn('Source abundance values are not recomputed', page)
        self.assertIn('id="runDetails"', page)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            db.row_factory = sqlite3.Row
            self.assertEqual(query(db, {})['columns'], ['sample_id','output_path'])
        js = (self.reports/'portal.js').read_text()
        self.assertIn("params.get('table')||'samples'", js)

    def write(self, sample, relative, value):
        path=self.root/sample/relative
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(value)

    def fixture(self,sample):
        self.write(sample,'preprocessed/test.mapping.txt','C1\toriginal1\t200\n')
        self.write(sample,'results/annotation_table/test.functional_and_taxonomic_table.txt',
            'ORF_ID\tORF_length\tstart\tend\tContig_Name\tContig_length\tstrand\ttarget\tproduct\ttaxonomy\n'
            f'G1\t90\t1\t90\tC1\t200\t+\tACC1\t{sample} kinase\tBacteria; Test\n'
            'G2\t90\t100\t190\tC1\t200\t-\tACC2\t</script><script>alert(1)</script>\troot\n')
        self.write(sample,'results/annotation_table/test.1.txt',
            'orf_id\tref dbname\ttarget\tproduct\tvalue\nG1\tswissprot\tACC1\tkinase\t12\n\tmetacyc\tACC2\tkinase\t5\n')
        self.write(sample,'results/annotation_table/test.EC_RXN_map.tsv',
            'orf_id\tref dbname\ttarget\tproduct\tvalue\tEC\tRXN\nG1\tswissprot\tACC1\tkinase\t12\t2.7.7.2|2.7.1.26\tRXN-1\nG1\tmetacyc\tACC2\tkinase\t5\t2.7.7.2\tRXN-2\n')
        self.write(sample,'ptools/orf_map.txt','G1\tG2\n')
        self.write(sample,'magsplitter/contig_to_mag.tsv','original1\tMAG1\n')
        self.write(sample,'magsplitter/results/MAG1/0.pf','ID\tG1\nNAME\tG1\n//\n')
        self.write(sample,'magsplitter/results/MAG_failed/0.pf','ID\tG2\n//\n')
        header='SAMPLE\tPWY_NAME\tPWY_COMMON_NAME\tPWY_SCORE\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\tORFS\n'
        self.write(sample,'results/pgdb/community/tag_pwy.tsv',header+'tag\tPWY-1\tTest pathway\t0.5\t4\t2\t2\tG1,G2\n')
        self.write(sample,'results/pgdb/MAGs/MAG1/tag_pwy.tsv',header+'MAG1\tPWY-1\tTest pathway\t0.5\t4\t2\t1\tG1\n')
        self.write(sample,'results/rpkm/test.orf_counts.tsv','Gene_ID\tCount\tRPKM\tTPM\nG1\t10\t1.2\t3.4\n')
        self.write(sample,'logs/ptools/run1/summary.json',json.dumps({'tasks':[{'task':'pgdb-MAG_failed','label':'MAG_failed','status':'FAILED','error':'Expected optional failure'}]}))

    def test_join_keys_and_no_annotation_cartesian_product(self):
        self.assertFalse(self.db.execute('PRAGMA foreign_key_check').fetchall())
        self.assertEqual(query(self.db,{'table':'orf_explorer'})['total'],4)
        self.assertEqual(query(self.db,{'table':'annotation_explorer'})['total'],4)
        self.assertEqual(query(self.db,{'table':'pathway_gene_explorer'})['total'],6)
        result=query(self.db,{'table':'pathway_gene_explorer','filters':[{'column':'sample_id','op':'eq','value':'alpha'},{'column':'entity_id','op':'eq','value':'MAG1'}]})
        self.assertEqual(result['total'],1)
        self.assertEqual(result['rows'][0]['orf_id'],'G1')
        self.assertEqual(query(self.db,{'table':'mag_gene_explorer'})['total'],4)
        self.assertEqual(query(self.db,{'table':'mag_orf_explorer'})['total'],4)
        failed=query(self.db,{'table':'entities','filters':[{'column':'entity_id','op':'eq','value':'MAG_failed'}]})
        self.assertEqual(failed['rows'][0]['last_task_status'],'FAILED')
        self.assertEqual(failed['rows'][0]['pathway_status'],'unavailable')

    def test_combined_workflow_statuses_preserve_sample_identity(self):
        self.write('', 'logs/analysis_wf/run2/summary.json', json.dumps({'tasks': [
            {'task': 'alpha:pgdb:MAG1', 'label': 'alpha:pgdb:MAG1', 'sample': 'alpha',
             'entity': 'MAG1', 'status': 'SKIPPED'},
            {'task': 'beta:pgdb:MAG1', 'label': 'beta:pgdb:MAG1', 'sample': 'beta',
             'entity': 'MAG1', 'status': 'SUCCESS'}]}))
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            self.assertEqual(db.execute("SELECT sample_id,last_task_status,pathway_status FROM entities WHERE entity_id='MAG1' ORDER BY sample_id").fetchall(),
                             [('alpha', 'SKIPPED', 'unavailable'), ('beta', 'SUCCESS', 'available')])
            self.assertEqual(db.execute("SELECT sample_id FROM pathways WHERE entity_id='MAG1'").fetchall(), [('beta',)])

    def test_related_filters_preserve_rows_and_same_annotation(self):
        spec={'table':'orf_explorer','related':{'pathway_id':'PWY-1','entity_id':'MAG1','ec':'2.7.7.2','reference_db':'swissprot'}}
        self.assertEqual(query(self.db,spec)['total'],2)
        spec['related']['reaction']='RXN-2'
        self.assertEqual(query(self.db,spec)['total'],0)

    def test_projection_numeric_filter_and_literal_search(self):
        result=query(self.db,{'table':'annotation_explorer','search':'kinase','columns':['orf_id','score'],
            'filters':[{'column':'score','op':'ge','value':'10'}]})
        self.assertEqual(result['total'],2)
        self.assertEqual(set(result['rows'][0]),{'orf_id','score'})
        self.assertEqual(query(self.db,{'table':'orf_explorer','search':'%'})['total'],0)
        self.assertEqual(query(self.db,{'table':'orf_explorer','search':'</script>'})['total'],2)

    def test_query_validation_and_read_only_csv(self):
        for spec in ({'table':'orfs; DROP TABLE samples'}, {'columns':['";DROP TABLE samples']},
                     {'filters':[{'column':'orf_id','op':'raw','value':'1=1'}]}, {'limit':0}):
            with self.assertRaises(ValueError): query(self.db,spec)
        columns, cursor=query(self.db,{'table':'orf_explorer','columns':['sample_id','orf_id'],'limit':1},export=True)
        self.assertEqual(columns,['sample_id','orf_id'])
        self.assertEqual(len(list(cursor)),4)
        self.assertEqual(csv_value('=SUM(A1)'),"'=SUM(A1)")
        self.assertEqual(csv_value(-3.4),-3.4)

    def test_missing_annotations_retained_and_report_refresh_does_not_touch_sources(self):
        source=self.root/'alpha/results/annotation_table/test.functional_and_taxonomic_table.txt'
        source.unlink()
        extra=self.reports/'user-notes.txt'; extra.write_text('retain')
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM orfs WHERE sample_id='alpha' AND annotation_present=0").fetchone()[0],2)
            self.assertGreater(db.execute('SELECT COUNT(*) FROM issues').fetchone()[0],0)
        self.assertEqual(extra.read_text(),'retain')
        self.assertIn('EDA_portal.html',(self.reports/'MP_run_report.html').read_text())

    def test_mag_period_normalization_matches_splitter_directory(self):
        self.write('alpha','magsplitter/contig_to_mag.tsv','original1\tMAG.1\n')
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            row=db.execute("SELECT entity_id,original_mag_id FROM contig_mags WHERE sample_id='alpha'").fetchone()
            self.assertEqual(row,('MAG_1','MAG.1'))
            self.assertEqual(db.execute("SELECT COUNT(*) FROM mag_orf_explorer WHERE sample_id='alpha' AND entity_id='MAG_1'").fetchone()[0],2)

    def test_compact_hit_fallback_fills_blank_orf_identifiers(self):
        (self.root/'alpha/results/annotation_table/test.EC_RXN_map.tsv').unlink()
        build_report(self.root)
        with sqlite3.connect(self.reports/'results.sqlite') as db:
            self.assertEqual(db.execute("SELECT COUNT(*) FROM annotations WHERE sample_id='alpha' AND orf_id='G1'").fetchone()[0],2)

    def test_http_query_export_and_path_restrictions(self):
        server=ReportServer(self.root)
        self.addCleanup(server.server_close)
        worker=threading.Thread(target=server.serve_forever,daemon=True);worker.start()
        self.addCleanup(server.shutdown)
        prefix=server.url.rsplit('/reports/',1)[0]
        spec={'table':'orf_explorer','columns':['sample_id','orf_id'],'filters':[{'column':'sample_id','op':'eq','value':'alpha'}],'limit':1}
        with urlopen(prefix+'/api/query?'+urlencode({'spec':json.dumps(spec)})) as response:
            data=json.load(response)
        self.assertEqual(data['total'],2)
        self.assertEqual(len(data['rows']),1)
        with urlopen(prefix+'/api/export?'+urlencode({'spec':json.dumps(spec)})) as response:
            records=list(csv.reader(io.StringIO(response.read().decode())))
        self.assertEqual(len(records),3)
        for path in ('/%2e%2e/etc/passwd','/.secret','/unknown'):
            with self.assertRaises(HTTPError):urlopen(prefix+path)
        with self.assertRaises(HTTPError):
            urlopen(Request(prefix+'/api/meta',headers={'Origin':'https://example.com'}))
        with urlopen(server.url) as response:
            self.assertIn(b'Explore and export',response.read())


if __name__=='__main__':unittest.main()
