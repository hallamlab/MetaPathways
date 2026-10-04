import csv
import io
import shlex
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from metapathways.LCAComputation import LCAComputation
from metapathways.MetaPathways_create_reports_fast import create_annotation, print_orf_table
from metapathways.protein_taxonomy import hit_taxid, hit_taxonomy, raw_lca


def tree():
    lca = LCAComputation([])
    lca.id_to_name = {'1': 'root', '2': 'Bacteria', '3': 'Species A', '4': 'Species B', '5': 'Archaea'}
    lca.name_to_id = {v: k for k, v in lca.id_to_name.items()}
    lca.taxid_to_ptaxid = {i: [p, 0, 0] for i, p in [('1','1'), ('2','1'), ('3','2'), ('4','2'), ('5','1')]}
    return lca


class TaxonomyStatusTests(unittest.TestCase):
    def test_rna_only_and_empty_ptinput_keep_coordinate_schema(self):
        from metapathways.MetaPathways_create_genbank_ptinput import ptinput_dataframe
        rna = {'rna1': {'id': 'rna1', 'seqname': 'contig1',
                       'start': 2, 'end': 9, 'strand': '+'}}
        frame = ptinput_dataframe(rna, {'contig1': 'A' * 20})
        self.assertEqual(frame.loc['rna1', 'contig_length'], 20)
        empty = ptinput_dataframe({}, {})
        self.assertTrue(empty.empty)
        self.assertTrue({'id', 'seqname', 'start', 'end', 'strand', 'contig_length'} <= set(empty))

    def test_orf_map_preserves_unannotated_cds_across_batches(self):
        output = io.StringIO()
        with tempfile.TemporaryDirectory() as directory:
            print_orf_table({'swissprot_test': {}}, {'C1-G1': 'sample-C1'}, directory, output)
            print_orf_table({'swissprot_test': {'C2-G1': [
                {'query': 'C2-G1', 'product': 'enzyme'}]}},
                {'C2-G1': 'sample-C2', 'C2-G2': 'sample-C2'}, directory, output)
        rows = list(csv.reader(io.StringIO(output.getvalue()), delimiter='\t'))
        self.assertEqual(rows, [['# ORF_ID', 'CONTIG_ID', 'swissprot_test'],
            ['C1-G1', 'sample-C1', ''], ['C2-G1', 'sample-C2', 'enzyme'],
            ['C2-G2', 'sample-C2', '']])

    def test_empty_orf_map_has_header(self):
        output = io.StringIO()
        with tempfile.TemporaryDirectory() as directory:
            print_orf_table({'swissprot_test': {}}, {}, directory, output)
        self.assertEqual(output.getvalue(), '# ORF_ID\tCONTIG_ID\tswissprot_test\n')

    def test_report_task_lists_selected_inputs_and_checkpoints_taxonomy(self):
        from metapathways.jobscreator import ContextCreator
        creator = ContextCreator.__new__(ContextCreator)
        creator.configs = SimpleNamespace(REFDBS='/db', CREATE_ANNOT_REPORTS='report-script')
        creator.params = MagicMock()
        creator.params.get.return_value = 'yes'
        creator.get_dbs = lambda: ['swissprot', 'metacyc']
        sample = SimpleNamespace(genbank_dir='/out/genbank', sample_name='sample',
            output_results_annotation_table_dir='/out/annotations', blast_results_dir='/out/blast', algorithm='FAST')
        context = creator.create_report_files_cmd(sample)[0]
        command = shlex.split(context.commands[0])
        self.assertNotIn('-D', command)
        self.assertEqual([command[i+1] for i,x in enumerate(command) if x == '-d'], ['swissprot','metacyc'])
        self.assertEqual(context.outputs['annotation_taxonomy'], '/out/annotations/sample.annotation_taxonomy.tsv')

    def test_taxid_header_forms(self):
        for text in ('Protein OS=Name OX=3 GN=g', 'Protein OS Name OX 3 GN g'):
            self.assertEqual(hit_taxid({'product': text}, 'swissprot'), '3')
        self.assertEqual(hit_taxid({'product':'protein', 'comment':'OS=Name OX=3'}, 'swissprot_test'), '3')
        self.assertEqual(hit_taxid({'product':'protein TaxID=4'}, 'uniref50'), '4')
        self.assertEqual(hit_taxid({'target':'5.protein'}, 'eggnog'), '5')
        self.assertIsNone(hit_taxid({'product':'OX=bad'}, 'swissprot'))
        self.assertIsNone(hit_taxid({'product':'OX=3bad'}, 'swissprot'))

    def test_missing_unknown_and_disabled_are_distinct(self):
        lca = tree()
        self.assertEqual(hit_taxonomy({'product':'OX=999'}, 'swissprot', lca), ('999','Unclassified'))
        self.assertEqual(hit_taxonomy({}, 'swissprot', lca), ('','Unclassified'))
        self.assertEqual(hit_taxonomy({'product':'OX=3'}, 'metacyc', lca), ('','Not computed'))
        self.assertEqual(hit_taxonomy({'product':'OX=1'}, 'swissprot', lca), ('1','root'))

    def test_lca_thresholds_and_real_root(self):
        lca = tree()
        hits = [{'product':'OX=3', 'bitscore':100}, {'product':'OX=4', 'bitscore':95},
                {'product':'OX=5', 'bitscore':20}, {'product':'OX=999', 'bitscore':100}]
        self.assertEqual(raw_lca(hits, 'swissprot', lca), '2')
        self.assertEqual(raw_lca([{'product':'OX=3','bitscore':100}, {'product':'OX=5','bitscore':100}], 'swissprot', lca), '1')
        self.assertIsNone(raw_lca([{'product':'OX=3','bitscore':20}], 'swissprot', lca))
        self.assertIsNone(raw_lca([{'product':'OX=999','bitscore':100}], 'swissprot', lca))
        self.assertEqual(raw_lca(hits, 'swissprot', lca), '2')  # LCA scratch counters were cleared.

    def test_report_generation_is_independent_of_database_order(self):
        from metapathways.MetaPathways_create_reports_fast import main
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ncbi = root/'tree.tsv'
            ncbi.write_text('root\t1\t1\nBacteria\t2\t1\nSpecies A\t3\t2\nSpecies B\t4\t2\nArchaea\t5\t1\n')
            gff = root/'annot.gff'
            gff.write_text(''.join(
                f'C1\tprodigal\tCDS\t{start}\t{start+299}\t.\t+\t0\tID=G{i};orf_length=300;contig_length=1000;sourcedb=swissprot;target=S{i};product=protein\n'
                for i,start in [(1,1),(2,400)]))
            header = '#query\ttarget\tq_length\tbitscore\tbsr\texpect\taln_length\tidentity\tec\tproduct\n'
            swiss = root/'swiss.parsed.txt'
            swiss.write_text(header+''.join(f'G{i}\tS{i}\t300\t100\t1\t0\t300\t100\t\tProtein OS=Name OX={i+2} GN=g\n' for i in (1,2)))
            uniref = root/'uniref.parsed.txt'
            uniref.write_text(header+''.join(f'G{i}\tU{i}\t300\t100\t1\t0\t300\t100\t\tProtein n=1 Tax=Archaea TaxID=5 RepID=foo\n' for i in (1,2)))
            eggnog = root/'eggnog.parsed.txt'
            eggnog.write_text(header+'G1\t3.protein\t300\t100\t1\t0\t300\t100\t\tprotein\n')
            expected = None
            for iteration, databases in enumerate(([('swissprot',swiss),('uniref50',uniref),('eggnog',eggnog)], [('eggnog',eggnog),('uniref50',uniref),('swissprot',swiss)])):
                output = root/str(iteration)
                args = ['--input-annotated-gff',str(gff),'--output-dir',str(output),'--ncbi-taxonomy-map',str(ncbi),'-s','sample']
                for db, filename in databases:
                    args += ['-d',db,'-b',str(filename)]
                main(args)
                with (output/'sample.annotation_taxonomy.tsv').open() as stream:
                    rows = sorted((r['orf_id'],r['reference_db'],r['taxonomy'],r['lca_taxonomy']) for r in csv.DictReader(stream, delimiter='\t'))
                self.assertIn(('G1','swissprot','Species A','Bacteria'), rows)
                self.assertIn(('G1','uniref50','Archaea','Archaea'), rows)
                self.assertIn(('G1','eggnog','Species A','root'), rows)  # Support must not leak from SwissProt.
                if expected is not None:
                    self.assertEqual(rows, expected)
                expected = rows

    def test_primary_annotation_uses_its_own_database_and_target(self):
        lca = tree()
        reader = MagicMock()
        reader.__iter__.return_value = iter(['contig1'])
        reader.orf_dictionary = {'contig1': [{
            'id': 'orf1', 'seqname': 'contig1', 'feature': 'CDS', 'sourcedb':'swissprot',
            'product': 'protein', 'orf_length':90, 'start':1, 'end':90,
            'contig_length':100, 'strand':'+', 'target':'hit1',
        }]}
        hits = {'swissprot': {'orf1':[{'target':'hit1','product':'OX=3'}]},
                'uniref50': {'orf1':[{'target':'hit1','product':'TaxID=5'}]}}
        taxons = {'swissprot':{'orf1':'Bacteria'}, 'uniref50':{'orf1':'Archaea'}}
        with tempfile.TemporaryDirectory() as directory:
            with patch('metapathways.MetaPathways_create_reports_fast.mputils.GffFileParser', return_value=reader):
                create_annotation(hits, list(hits), 'unused.gff', directory, taxons,
                                  {'orf1': True}, {}, lca, sample_name='sample')
            row = (Path(directory)/'sample.functional_and_taxonomic_table.txt').read_text().strip().split('\t')
        self.assertEqual(row[-3:], ['Species A','swissprot','Bacteria'])


if __name__ == '__main__':
    unittest.main()
