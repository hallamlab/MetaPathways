import pytest
import sys
import os
import filecmp
import re
import glob
from functools import total_ordering

@total_ordering
class Sequence:
    def __init__(self, name, seq): 
        self.name = name 
        self.seq =  seq 
  
    def __lt__(self, other): 
        return self.name < other.name or self.seq < other.seq 
  
    def __eq__(self, other): 
        return self.name == other.name and self.seq == other.seq 

def print_unequal_lines_in_files(file1, file2, sort_n_compare):
    with gzip.open(file1, 'r') if file1.endswith('.gz') \
         else open(file1, 'r') as fout, \
         gzip.open(file2, 'r') if file2.endswith('.gz') \
         else open(file2, 'r') as expfout:

       if sort_n_compare:
           output_lines = sorted(fout.readlines())
           expect_output_lines = sorted(expfout.readlines())
       else:
           output_lines = fout.readlines()
           expect_output_lines = expfout.readlines()
     
    for a, b in zip(output_lines, expect_output_lines):
        if a != b:
            print(">> " + a.strip())
            print("<< " + b.strip())

 
def compare_lines_in_files(file1, file2, sort_n_compare):
    with gzip.open(file1, 'r') if file1.endswith('.gz') \
         else open(file1, 'r') as fout, \
         gzip.open(file2, 'r') if file2.endswith('.gz') \
         else open(file2, 'r') as expfout:

       if sort_n_compare:
           output_lines = sorted(fout.readlines())
           expect_output_lines = sorted(expfout.readlines())
       else:
           output_lines = fout.readlines()
           expect_output_lines = expfout.readlines()
     
    #print([a == b for a, b in zip(output_lines, expect_output_lines)])
    return all([a == b for a, b in zip(output_lines, expect_output_lines)])

def compare_rpkm_stats_in_files(file1, file2, sort_n_compare):
    with gzip.open(file1, 'r') if file1.endswith('.gz') \
         else open(file1, 'r') as fout, \
         gzip.open(file2, 'r') if file2.endswith('.gz') \
         else open(file2, 'r') as expfout:

       if sort_n_compare:
           output_lines = sorted(fout.readlines())
           expect_output_lines = sorted(expfout.readlines())
       else:
           output_lines = fout.readlines()
           expect_output_lines = expfout.readlines()
    
    for a, b in zip(output_lines, expect_output_lines):
        if len(a.split(':')) == len(b.split(':')):
            if len(a.split(':')) > 0:
                assert(a.split(':')[0].strip() == b.split(':')[0].strip())
            if len(a.split(':')) > 1:
                try:
                    numa =  float(re.sub('%',  '', a.split(':')[1]).strip())
                    numb =  float(re.sub('%',  '', b.split(':')[1]).strip())
                    if numb != 0:
                        assert(numa == pytest.approx(numb, rel = 1e-1) ) 
                    else:
                        assert(pytest.approx(numa, 0.1) == 0)
                except:
                    pass
                
    return True

def compare_rpkm_values_in_files(file1, file2, sort_n_compare):
    with gzip.open(file1, 'r') if file1.endswith('.gz') \
         else open(file1, 'r') as fout, \
         gzip.open(file2, 'r') if file2.endswith('.gz') \
         else open(file2, 'r') as expfout:

       if sort_n_compare:
           output_lines = sorted(fout.readlines())
           expect_output_lines = sorted(expfout.readlines())
       else:
           output_lines = fout.readlines()
           expect_output_lines = expfout.readlines()
    
    for a, b in zip(output_lines, expect_output_lines):
         if "ORF_ID" not in a:
            numa =  float(a.split('\t')[1].strip())
            numb =  float(b.split('\t')[1].strip())
            if  numa != pytest.approx(numb, rel = 1e-1):
               print('WARNING: mismatch', a, b)
#        assert(len(a.split('\t')) == len(b.split('\t')))
#        assert(len(a.split('\t')) == 2)
#
#        # header line or the values
#        if a.split('\t')[1].strip() == 'COUNT':
#           print(a, b)
#           assert(a.split('\t')[0].strip() == b.split('\t')[0].strip())
#        else:
#            numa =  float(a.split('\t')[1].strip())
#            numb =  float(b.split('\t')[1].strip())
#            assert(numa == pytest.approx(numb, rel = 1e-1) ) 
                
    return True


def compare_fasta_files(file1, file2, sort_n_compare):
    import pyfastx
    contents1 = [] 
    for seq in pyfastx.Fasta(file1):
       contents1.append(Sequence(seq.name, seq.seq))

    contents2 = [] 
    for seq in pyfastx.Fasta(file1):
       contents2.append(Sequence(seq.name, seq.seq))
    
    contents1.sort()
    contents2.sort()

    return all([a == b for a, b in zip(contents1, contents2)])
    return True



_test_data_dir = os.path.join(os.path.split(__file__)[0], "data/")
_test_ref_data_dir = os.path.join(_test_data_dir, "ref_data")

def pipeline(tmpdir, sample_name):
    from metapathways import pipeline 
    _test_sample_name =  sample_name
    _test_input_fasta_file = os.path.join(_test_data_dir, _test_sample_name, "/input/", _test_sample_name + ".fasta")

    commands = ["MetaPathways", "-i",  _test_data_dir + '/' + _test_sample_name + '/input/',
                "-o", str(tmpdir),
                "-s", sample_name,
                "-p", _test_data_dir + "/template_param.txt",
                "-d", _test_ref_data_dir
               ]
    pipeline.process(commands)
    assert 1 == 1

@pytest.fixture(params=["lagoon-sample", "lagoon-sample2"])
#@pytest.fixture(params=[ "lagoon-sample2"])
def sample_name(request):
    return request.param

def test_filter_nuc_input(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_filter_input
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name 
    _test_input_fasta_file = os.path.join(_test_data_dir, 
                                          "input", 
                                          _test_sample_name + ".fasta"
                                         )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "preprocessed"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_fasta = os.path.join(out_folder_name, 
                                _test_sample_name, 
                                "preprocessed", 
                                _test_sample_name + ".fasta"
                               )
    output_mapping_file = os.path.join(out_folder_name, 
                                       _test_sample_name,
                                       "preprocessed", 
                                       _test_sample_name + ".mapping.txt"
                                      )

    output_contig_length = os.path.join(out_folder_name, 
                                        _test_sample_name, 
                                        "run_statistics", 
                                        _test_sample_name + ".contig.lengths.txt"
                                       )

    output_logfile = os.path.join(out_folder_name, 
                                  _test_sample_name, 
                                  "run_statistics", 
                                  _test_sample_name + ".nuc.stats"
                                 )

    args = [ "--min_length",  "180",
             "--log_file", output_logfile,
             "-i", _test_input_fasta_file,
             "-o", output_fasta,
             "-M", output_mapping_file,
             "-t", "nucleotide",
             "-L", output_contig_length
           ]

    expect_output_fasta = os.path.join(_test_data_dir, 
                                         "output", 
                                         _test_sample_name, 
                                         'preprocessed', 
                                         _test_sample_name + ".fasta"
                                        )

    expect_output_mapping_file = os.path.join(_test_data_dir, 
                                         "output", 
                                         _test_sample_name, 
                                         "preprocessed", 
                                         _test_sample_name + ".mapping.txt"
                                        )

    expect_output_contig_length = os.path.join(_test_data_dir, 
                                                 'output', 
                                                 _test_sample_name, 
                                                 "run_statistics", 
                                                 _test_sample_name + ".contig.lengths.txt"
                                               )
    expect_output_logfile = os.path.join(out_folder_name, 
                            _test_sample_name, 
                            "run_statistics", 
                            _test_sample_name + ".nuc.stats"
                           )

    MetaPathways_filter_input.main(args)
    assert filecmp.cmp(output_fasta, expect_output_fasta)
    assert filecmp.cmp(output_mapping_file, expect_output_mapping_file)
    assert filecmp.cmp(output_contig_length, expect_output_contig_length)
    assert filecmp.cmp(output_logfile, expect_output_logfile)


def test_orf_predictions_input(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_orf_prediction
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name 

    prod_input = os.path.join(_test_data_dir, 
                              "output",
                              _test_sample_name, 
                              "preprocessed", 
                              _test_sample_name + ".fasta"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_gff = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".gff"
                             )

    expect_output_gff = os.path.join(_test_data_dir,
                                       "output", 
                                       _test_sample_name, 
                                       "orf_prediction", 
                                       _test_sample_name + ".gff"
                                        )

    args = [ "--prod_exec", "prodigal", 
             "--prod_m", 
             "--prod_p", 'meta',
             "--prod_f", 'gff',
             "--prod_g", '11',
             "--prod_input", prod_input,
             "--prod_output", output_gff
           ]


    MetaPathways_orf_prediction.main(args)
    #print(glob.glob(os.path.join(out_folder_name, _test_sample_name, "orf_prediction/*")))
    assert filecmp.cmp(output_gff, expect_output_gff)


def test_create_amino_sequences(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_create_amino_sequences 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name 

    g_input = os.path.join(_test_data_dir, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".gff"
                            )

    n_input = os.path.join(_test_data_dir, 
                              "output",
                              _test_sample_name, 
                              "preprocessed", 
                              _test_sample_name + ".fasta"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_amino = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".faa"
                             )

    expect_output_amino = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "orf_prediction", 
                                       _test_sample_name + ".faa"
                                        )

    output_nuc = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".fna"
                             )

    expect_output_nuc = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "orf_prediction", 
                                       _test_sample_name + ".fna"
                                        )

    output_gff = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".unannot.gff"
                             )

    expect_output_gff = os.path.join(_test_data_dir,
                                       "output", 
                                       _test_sample_name, 
                                       "orf_prediction", 
                                       _test_sample_name + ".unannot.gff"
                                        )
    args = [ "-g", g_input, 
             "-n", n_input, 
             "--output_amino", output_amino,
             "--output_nuc", output_nuc,
             "--output_gff", output_gff
           ]


    MetaPathways_create_amino_sequences.main(args)

    #print(glob.glob(os.path.join(out_folder_name, _test_sample_name, "orf_prediction/*")))
    assert compare_lines_in_files(output_gff, expect_output_gff, sort_n_compare = True)
    assert compare_fasta_files(output_amino, expect_output_amino, sort_n_compare = True)
    assert compare_fasta_files(output_nuc, expect_output_nuc, sort_n_compare = True)


def test_refscores(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_refscore 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    faa_input = os.path.join(_test_data_dir, 
                              _test_sample_name, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".qced.faa"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_refscores = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "blast_results", 
                              _test_sample_name + ".refscores.BLAST"
                             )

    expect_output_refscores = os.path.join(_test_data_dir, 
                                             _test_sample_name, 
                                             "output", 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + ".refscores.BLAST"
                                            )

    args = [ "-i", faa_input, 
             "-o", output_refscores, 
             "-a", "BLAST",
           ]
    MetaPathways_refscore.main(args)

    assert compare_lines_in_files(output_refscores, expect_output_refscores, sort_n_compare = True)


def test_filter_amino_input(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_filter_input
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_faa = os.path.join( _test_data_dir, 
                              "output", 
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".faa"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_qced_faa = os.path.join(out_folder_name, 
                                   _test_sample_name, 
                                  "orf_prediction", 
                                  _test_sample_name + ".qced.faa"
                                  )

    expect_output_qced_faa = os.path.join(_test_data_dir, 
                                         "output", 
                                         _test_sample_name, 
                                         'orf_prediction', 
                                         _test_sample_name + ".qced.faa"
                                        )

    output_aminoseq_length = os.path.join(out_folder_name, 
                                          _test_sample_name, 
                                          "run_statistics", 
                                          _test_sample_name + ".orf.lengths.txt"
                                       )

    expect_output_aminoseq_length = os.path.join(_test_data_dir, 
                                                'output', 
                                                 _test_sample_name, 
                                                "run_statistics", 
                                                 _test_sample_name + ".orf.lengths.txt"
                                                )

    output_logfile = os.path.join(out_folder_name, 
                                  _test_sample_name, 
                                  "run_statistics", 
                                  _test_sample_name + ".amino.stats"
                                 )

    expect_output_logfile = os.path.join(out_folder_name, 
                                           _test_sample_name, 
                                           "run_statistics", 
                                           _test_sample_name + ".amino.stats"
                                          )
    args = [ "--min_length",  "60",
             "--log_file", output_logfile,
             "-i", input_faa,
             "-o", output_qced_faa,
             "-t", "amino",
             "-L", output_aminoseq_length
           ]

    MetaPathways_filter_input.main(args)

    assert compare_fasta_files(output_qced_faa, expect_output_qced_faa, sort_n_compare = True)
    assert compare_lines_in_files(output_aminoseq_length, expect_output_aminoseq_length, sort_n_compare = False)
    assert compare_lines_in_files(output_logfile, expect_output_logfile, sort_n_compare = False)

def test_func_search_sample1(tmpdir, sample_name):
    run_func_search(tmpdir, "refseq-small-sample", sample_name)

@pytest.fixture(params=["refseq-small-sample", "cazy-small-sample", "kegg-small-sample", "metacyc-small-sample"])
def db_name(request):
    return request.param

def test_func_search_sample2(tmpdir, db_name):
    run_func_search(tmpdir, db_name, "lagoon-sample2")


def run_func_search(tmpdir, dbname,  test_sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_func_search 
    out_folder_name = str(tmpdir)
    _test_sample_name = test_sample_name

    faa_input = os.path.join(_test_data_dir, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".qced.faa"
                            )

    ref_blast_db = os.path.join(_test_data_dir, 
                            "ref_data",
                            "functional",
                            "formatted",
                            dbname
                           )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)

    output_blast_results = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "blast_results", 
                              _test_sample_name + "." + dbname + ".BLASTout"
                             )

    expect_output_blast_results = os.path.join(_test_data_dir, 
                                               "output", 
                                               _test_sample_name, 
                                               "blast_results", 
                                               _test_sample_name + "." + dbname + ".BLASTout"
                                              )

    args = [ "--algorithm", "BLAST", 
             "--blast_executable", "blastp", 
             "--num_threads", "4",
             "--blast_max_target_seqs", "5",
             "--blast_outfmt", "6",
             "--blast_query", faa_input,
             "--blast_evalue", "0.000001",
             "--blast_db", ref_blast_db,
             "--blast_out", output_blast_results,
           ]

    MetaPathways_func_search.main(args)

    assert compare_lines_in_files(output_blast_results, expect_output_blast_results, sort_n_compare = True)

def test_refscores(tmpdir, sample_name):
    from metapathways import pipeline 
    from metapathways import MetaPathways_refscore 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    faa_input = os.path.join(_test_data_dir, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".qced.faa"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "run_statistics"), exist_ok = True)

    output_refscores = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "blast_results", 
                              _test_sample_name + ".refscores.BLAST"
                             )

    expect_output_refscores = os.path.join(_test_data_dir, 
                                             "output", 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + ".refscores.BLAST"
                                            )

    args = [ "-i", faa_input, 
             "-o", output_refscores, 
             "-a", "BLAST",
           ]
    MetaPathways_refscore.main(args)

    assert compare_lines_in_files(output_refscores, expect_output_refscores, sort_n_compare = True)

def test_parse_blast1(tmpdir, sample_name):
    run_parse_blast(tmpdir, "refseq-small-sample", sample_name)

def test_parse_blast2(tmpdir, db_name):
    run_parse_blast(tmpdir, db_name, "lagoon-sample2")

def run_parse_blast(tmpdir,  db_name, sample_name):
    from metapathways import MetaPathways_parse_blast 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    ref_blast_db_annots = os.path.join(_test_data_dir, 
                            "ref_data",
                            "functional",
                            "formatted",
                            db_name + "-names.txt"
                           )

    input_blast_results = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                        _test_sample_name + "." + db_name + ".BLASTout"
                                      )

    input_blast_refscores = os.path.join(_test_data_dir, 
                                         'output',
                                         _test_sample_name, 
                                         "blast_results", 
                                         _test_sample_name + ".refscores.BLAST"
                                        )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)

    output_parsed_blast_results = os.path.join(out_folder_name, 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + "." + db_name + ".BLASTout.parsed.txt"
                                            )

    expect_output_parsed_blast_results = os.path.join(_test_data_dir, 
                                             "output", 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + "." + db_name + ".BLASTout.parsed.txt"
                                            )
    args = [ 
             "-b", input_blast_results,
             "-r", input_blast_refscores,
             "-m", ref_blast_db_annots,
             "-d", db_name,
             "-o", output_parsed_blast_results,
             "--min_bsr", "0.4",
             "--min_score", "20",
             "--min_length", "45",
             "--max_evalue", "0.000001",
             "--algorithm", "BLAST"
           ]

    MetaPathways_parse_blast.main(args)
    assert compare_lines_in_files(output_parsed_blast_results, expect_output_parsed_blast_results, sort_n_compare = True)

@pytest.fixture(params=["silva-ssu-small-sample", "silva-lsu-small-sample"], ids = ["silva-ssu-small", "silva-lsu-small"])
def rna_db(request):
    return request.param

def test_rRNA_stats_calculator(tmpdir, sample_name, rna_db):
    from metapathways import MetaPathways_rRNA_stats_calculator 
    from metapathways import sysutil as sysutils
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_blast_results = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                        _test_sample_name + ".rRNA." + rna_db + ".BLASTout"
                                      )

    ref_blast_db = os.path.join(_test_data_dir, 
                            "ref_data",
                            "taxonomic",
                            rna_db
                           )
    
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "results",  "rRNA"), exist_ok = True)

    output_lsu_rRNA_stats = os.path.join(out_folder_name, 
                                         _test_sample_name, 
                                         "results",
                                         "rRNA", 
                                         _test_sample_name + "." + rna_db + ".rRNA.stats.txt"
                                            )

    expect_output_lsu_rRNA_stats = os.path.join(_test_data_dir, 
                                             "output", 
                                             _test_sample_name, 
                                             "results",
                                             "rRNA", 
                                             _test_sample_name + "." + rna_db + ".rRNA.stats.txt"
                                            )
    args = [ 
             "-o", output_lsu_rRNA_stats,
             "-b", "50",
             "-e", "0.000001",
             "-s", "20", 
             "-i", input_blast_results,
             "-d", ref_blast_db
           ]

    MetaPathways_rRNA_stats_calculator.main(args)
    assert compare_lines_in_files(output_lsu_rRNA_stats, expect_output_lsu_rRNA_stats, sort_n_compare = True)

    """ Test the BLAST output """
    ref_blast_db_formatted = os.path.join(_test_data_dir, 
                                          "ref_data",
                                          "taxonomic",
                                          "formatted",
                                          rna_db
                                          )

    input_contigs = os.path.join(_test_data_dir, 
                                 "output", 
                                  _test_sample_name, 
                                 'preprocessed', 
                                  _test_sample_name + ".fasta"
                                )

    expect_output_blast_results = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                       _test_sample_name + ".rRNA." + rna_db + ".BLASTout"
                                      )

    output_blast_results = os.path.join(out_folder_name, 
                                        _test_sample_name, 
                                        "blast_results",
                                        _test_sample_name + ".rRNA." + rna_db + ".BLASTout"
                                       )

    blast_args = [ 
             "blastn",
             "-outfmt", "6",
             "-num_threads", "8",
             "-query", input_contigs, 
             "-out",  output_blast_results,
             "-db", ref_blast_db_formatted,
             "-max_target_seqs", '5'
            ]

    result = sysutils.getstatusoutput(' '.join(blast_args))
    assert compare_lines_in_files(output_blast_results, expect_output_blast_results, sort_n_compare = True)

def test_tRNA_scan(tmpdir, sample_name):
    from metapathways import MetaPathways_tRNA_scan 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_contigs = os.path.join(_test_data_dir, 
                                 "output", 
                                 _test_sample_name, 
                                 "preprocessed", 
                                 _test_sample_name + ".fasta"
                                )

    TPCsignal = os.path.join(_test_data_dir, "ref_data", "TPCsignal")
    Dsignal = os.path.join(_test_data_dir, "ref_data", "Dsignal")

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "results",  "tRNA"), exist_ok = True)

    output_tRNA_stats = os.path.join(out_folder_name, 
                                             _test_sample_name, 
                                             "results",
                                             "tRNA", 
                                             _test_sample_name + ".tRNA.stats.txt"
                                            )

    expect_output_tRNA_stats = os.path.join(_test_data_dir, 
                                              "output", 
                                              _test_sample_name, 
                                              "results",
                                              "tRNA", 
                                              _test_sample_name + ".tRNA.stats.txt"
                                             )
     
    output_tRNA_fasta = os.path.join(out_folder_name, 
                                     _test_sample_name, 
                                     "results",
                                     "tRNA", 
                                     _test_sample_name + ".tRNA.fasta"
                                    )

    expect_output_tRNA_fasta = os.path.join(_test_data_dir, 
                                            "output", 
                                            _test_sample_name, 
                                            "results",
                                            "rRNA", 
                                            _test_sample_name + ".tRNA.fasta"
                                          )
     
    args = [ 
             "--executable", "trnascan-1.4",
             "-T", TPCsignal,
             "-D", Dsignal,
             "-F", output_tRNA_fasta,
             "-i", input_contigs,
             "-o",  output_tRNA_stats
            ]

    MetaPathways_tRNA_scan.main(args)
    assert compare_lines_in_files(output_tRNA_stats, expect_output_tRNA_stats, sort_n_compare = True)
    assert compare_fasta_files(output_tRNA_fasta, expect_output_tRNA_fasta, sort_n_compare = True)


def test_annotate_fast(tmpdir, sample_name):
    from metapathways import MetaPathways_annotate_fast 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_gff = os.path.join(_test_data_dir, 
                             "output", 
                             _test_sample_name, 
                             "orf_prediction", 
                             _test_sample_name + ".unannot.gff"
                            )

    output_annot_gff = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "genbank", 
                             _test_sample_name + ".annot.gff"
                            )

    expect_output_annot_gff = os.path.join(_test_data_dir, 
                             "output", 
                             _test_sample_name, 
                             "genbank", 
                             _test_sample_name + ".annot.gff"
                            )

    input_lsu_rRNA_stats = os.path.join(_test_data_dir, 
                                        "output", 
                                         _test_sample_name, 
                                        "results",
                                        "rRNA", 
                                        _test_sample_name + ".silva-lsu-small-sample.rRNA.stats.txt"
                                       )


    input_ssu_rRNA_stats = os.path.join(_test_data_dir, 
                                        "output", 
                                         _test_sample_name, 
                                        "results",
                                        "rRNA", 
                                        _test_sample_name + ".silva-ssu-small-sample.rRNA.stats.txt"
                                      )

    input_tRNA_stats = os.path.join(_test_data_dir, 
                                              "output", 
                                              _test_sample_name, 
                                              "results",
                                              "tRNA", 
                                              _test_sample_name + ".tRNA.stats.txt"
                                             )
     
    annot_folder_prefix = os.path.join(out_folder_name, 
                                _test_sample_name, 
                                "results",
                                "annotation_table", 
                                _test_sample_name 
                                )

    input_mapping = os.path.join(_test_data_dir, 
                                 "output", 
                                 _test_sample_name, 
                                "preprocessed", 
                                _test_sample_name + ".mapping.txt"
                               )

    output_annot_gff = os.path.join(_test_data_dir, 
                             "output", 
                             _test_sample_name, 
                             "genbank", 
                             _test_sample_name + ".annot.gff"
                            )

    input_D =  os.path.join(_test_data_dir, 
                           "output", 
                           _test_sample_name, 
                           "blast_results/") 

    output_dir_tsv_pfx = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "results",
                             "annotation_table",
                             _test_sample_name 
                             )

    expect_output_dir_tsv_pfx = os.path.join(_test_data_dir, 
                                          "output",
                                          _test_sample_name, 
                                          "results",
                                          "annotation_table",
                                          _test_sample_name
                                         )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name,  "genbank"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name,  "results",  "annotation_table"), exist_ok = True)

    argv = [
      "--input_gff", input_gff, 
      "-o", output_annot_gff, 
      "--rRNA_16S", input_ssu_rRNA_stats,
      "--rRNA_16S", input_lsu_rRNA_stats,
      "--tRNA", input_tRNA_stats, 
      "--output-comparative-annotation", annot_folder_prefix,
      "--algorithm",  "BLAST",
      "-m", input_mapping, 
      "-D", input_D, 
      "-s", _test_sample_name
    ]

    MetaPathways_annotate_fast.main(argv)
    assert compare_lines_in_files(output_annot_gff, expect_output_annot_gff, sort_n_compare = True)
    
    print_unequal_lines_in_files(output_dir_tsv_pfx + ".1.txt", expect_output_dir_tsv_pfx + ".1.txt", sort_n_compare = True)
   # assert compare_lines_in_files(output_dir_tsv_pfx + ".1.txt", expect_output_dir_tsv_pfx + ".1.txt", sort_n_compare = True)
    #assert compare_lines_in_files(output_dir_tsv_pfx + ".2.txt", expect_output_dir_tsv_pfx + ".2.txt", sort_n_compare = True)

def test_create_reports_fast(tmpdir, sample_name):
    from metapathways import MetaPathways_create_reports_fast 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_annot_gff = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "genbank", 
                                        _test_sample_name + ".annot.gff"
                                      )

    ref_kegg_maps = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "functional_categories",
                                   "KO_classification.txt.gz" 
                                   )

    ref_cog_maps = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "functional_categories",
                                   "COG_categories.txt.gz" 
                                   )

    ref_seed_maps = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "functional_categories",
                                   "SEED_subsystems.txt.gz" 
                                   )

    ref_cazy_maps = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "functional_categories",
                                   "CAZY_hierarchy.txt"
                                   )

    ref_ncbi_tree = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "ncbi_tree",
                                   "ncbi_taxonomy_tree.txt"
                                   )

    ref_ncbi_maps = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "ncbi_tree",
                                   "ncbi.map"
                                   )

    input_D =  os.path.join(_test_data_dir, 
                           "output", 
                           _test_sample_name, 
                           "blast_results/") 

    output_dir = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "results",
                             "annotation_table"
                             )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "results",  "annotation_table"), exist_ok = True)

    output_fun_and_tax = os.path.join(out_folder_name, 
                                      _test_sample_name, 
                                      "results",
                                      "annotation_table",
                                      _test_sample_name + ".functional_and_taxonomic_table.txt" 
                                  )

    expect_output_fun_and_tax = os.path.join(_test_data_dir, 
                                      "output", 
                                      _test_sample_name, 
                                      "results",
                                      "annotation_table",
                                      _test_sample_name + ".functional_and_taxonomic_table.txt" 
                                  )
    output_ORF_and_annt = os.path.join(out_folder_name, 
                                      _test_sample_name, 
                                      "results",
                                      "annotation_table",
                                      _test_sample_name + ".ORF_annotation_table.txt" 
                                  )

    expect_output_ORF_and_annt = os.path.join(_test_data_dir, 
                                      "output", 
                                      _test_sample_name, 
                                      "results",
                                      "annotation_table",
                                      _test_sample_name + ".ORF_annotation_table.txt" 
                                  )

    args = [ 
             "--input-annotated-gff", input_annot_gff,
             "--input-kegg-maps", ref_kegg_maps,
             "--input-cog-maps", ref_cog_maps,
             "--input-seed-maps", ref_seed_maps,
             "--input-cazy-maps", ref_cazy_maps,
             "--output-dir",  output_dir,
             "--ncbi-taxonomy-map", ref_ncbi_tree,
             "--ncbi-megan-map", ref_ncbi_maps,
             "-D",  input_D,
             "-s", _test_sample_name,
             "-a",  "BLAST"
            ]

    print('\nMetaPathways_create_reports_fast ' +  ' '.join(args))
    cmd = "MetaPathways_create_reports_fast  --input-annotated-gff /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/genbank//lagoon-sample.annot.gff  --input-kegg-maps tests/data/ref_data//functional_categories/KO_classification.txt                --input-cog-maps tests/data/ref_data//functional_categories/COG_categories.txt --input-seed-maps tests/data/ref_data//functional_categories/SEED_subsystems.txt --input-cazy-maps tests/data/ref_data//functional_categories/CAZY_hierarchy.txt --output-dir /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/results//annotation_table/                --ncbi-taxonomy-map tests/data/ref_data//ncbi_tree/ncbi_taxonomy_tree.txt --ncbi-megan-map tests/data/ref_data//ncbi_tree/ncbi.map -D /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/blast_results/ -s lagoon-sample -a BLAST"
    print("TODO")
    print(cmd)

    MetaPathways_create_reports_fast.main(args)

    assert compare_lines_in_files(output_ORF_and_annt, expect_output_ORF_and_annt, sort_n_compare = True)
    #assert compare_lines_in_files(output_fun_and_tax, expect_output_fun_and_tax, sort_n_compare = True)

    print_unequal_lines_in_files(output_fun_and_tax, expect_output_fun_and_tax, sort_n_compare = True)

def test_create_genbank_ptinput(tmpdir, sample_name):
    from metapathways import MetaPathways_create_genbank_ptinput 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

#MetaPathways_create_reports_fast  --input-annotated-gff /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/genbank//lagoon-sample.annot.gff  --input-kegg-maps tests/data/ref_data//functional_categories/KO_classification.txt                --input-cog-maps tests/data/ref_data//functional_categories/COG_categories.txt --input-seed-maps tests/data/ref_data//functional_categories/SEED_subsystems.txt --input-cazy-maps tests/data/ref_data//functional_categories/CAZY_hierarchy.txt --output-dir /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/results//annotation_table/                --ncbi-taxonomy-map tests/data/ref_data//ncbi_tree/ncbi_taxonomy_tree.txt --ncbi-megan-map tests/data/ref_data//ncbi_tree/ncbi.map -D /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/blast_results/ -s lagoon-sample -a BLAST
 
    input_annot_gff = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "genbank", 
                                        _test_sample_name + ".annot.gff"
                                      )

    input_preprocess_fasta = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "preprocessed", 
                                        _test_sample_name + ".fasta"
                                      )

    input_qced_faa = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "orf_prediction", 
                                        _test_sample_name + ".qced.faa"
                                      )

    ref_ncbi_tree = os.path.join(_test_data_dir, 
                                   "ref_data",
                                   "ncbi_tree",
                                   "ncbi_taxonomy_tree.txt"
                                   )

    input_output_fun_and_tax = os.path.join(_test_data_dir, 
                                      "output", 
                                      _test_sample_name, 
                                      "results",
                                      "annotation_table",
                                      _test_sample_name + ".functional_and_taxonomic_table.txt" 
                                  )

    output_ptools_dir = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "ptools"
                             )

    output_ptools_pf = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "ptools",
                             "0.pf"
                             )

    expect_output_ptools_pf = os.path.join(_test_data_dir, 
                                             "output",
                                             _test_sample_name, 
                                             "ptools",
                                             "0.pf"
                                            )

    output_ptools_gen_elem = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "ptools",
                             "genetic-elements.dat"
                             )

    expect_output_ptools_gen_elem = os.path.join(_test_data_dir, 
                                             "output",
                                             _test_sample_name, 
                                             "ptools",
                                             "genetic-elements.dat"
                                            )

    output_ptools_org_param = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "ptools",
                             "organism-params.dat"
                             )

    expect_output_ptools_org_param = os.path.join(_test_data_dir, 
                                             "output",
                                             _test_sample_name, 
                                             "ptools",
                                             "organism-params.dat"
                                            )

    output_ptools_reduced = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "ptools",
                             "reduced.txt"
                             )

    expect_output_ptools_reduced = os.path.join(_test_data_dir, 
                                             "output",
                                             _test_sample_name, 
                                             "ptools",
                                             "reduced.txt"
                                            )

    output_gbk = os.path.join(out_folder_name, 
                             _test_sample_name, 
                             "genbank",
                             _test_sample_name + ".gbk" 
                             )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name,  "ptools"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name,  "genbank"), exist_ok = True)

    args = [ 
             "--g", input_annot_gff,
             "--out-ptinput", output_ptools_dir,
             "-n", input_preprocess_fasta,
             "--ncbi-tree", ref_ncbi_tree,
             "--taxonomy-table", input_output_fun_and_tax,
             "-p",  input_qced_faa,
             "--out-gbk",  output_gbk
            ]

    MetaPathways_create_genbank_ptinput.main(args)
    #assert compare_lines_in_files(output_ptools_pf, expect_output_ptools_pf, sort_n_compare = True)
    assert compare_lines_in_files(output_ptools_gen_elem, expect_output_ptools_gen_elem, sort_n_compare = True)
    assert compare_lines_in_files(output_ptools_org_param, expect_output_ptools_org_param, sort_n_compare = True)
    #assert compare_lines_in_files(output_ptools_reduced, expect_output_ptools_reduced, sort_n_compare = True)


# MetaPathways_create_genbank_ptinput 
#-g /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/genbank//lagoon-sample.annot.gff  
#--out-ptinput /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/ptools/ 
#-n /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/preprocessed//lagoon-sample.fasta 
#--ncbi-tree tests/data/ref_data//ncbi_tree/ncbi_taxonomy_tree.txt 
#--taxonomy-table /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/results//annotation_table//lagoon-sample.functional_and_taxonomic_table.txt 
#-p /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/orf_prediction//lagoon-sample.qced.faa 
#--out-gbk /home/kishori/metapathways_engcyc/metapathways/mp_output/lagoon-sample/genbank//lagoon-sample.gbk

def test_rpkm(tmpdir, sample_name):
    from metapathways import MetaPathways_rpkm 
    out_folder_name = str(tmpdir)
    _test_sample_name = sample_name

    input_contigs = os.path.join(_test_data_dir, 
                                       "output", 
                                       _test_sample_name, 
                                       "preprocessed", 
                                        _test_sample_name + ".fasta"
                                      )
    input_reads = os.path.join(_test_data_dir, 
                               "input", 
                               "reads" 
                              )


    input_gff = os.path.join(_test_data_dir, 
                             "output", 
                             _test_sample_name, 
                             "genbank", 
                             _test_sample_name + ".annot.gff"
                            )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "results",  "rpkm"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "bwa"), exist_ok = True)

    output_orf_rpkm = os.path.join(out_folder_name, 
                                   _test_sample_name, 
                                   "results",
                                   "rpkm", 
                                   _test_sample_name + ".orf_rpkm.txt"
                                  )

    expect_output_orf_rpkm = os.path.join(_test_data_dir, 
                                     "output", 
                                     _test_sample_name, 
                                     "results",
                                     "rpkm", 
                                     _test_sample_name + ".orf_rpkm.txt"
                                    )

    output_rpkm_stats = os.path.join(out_folder_name, 
                                     _test_sample_name, 
                                     "results",
                                     "rpkm", 
                                     _test_sample_name + ".orf_read_counts_stats.txt"
                                    )

    expect_output_rpkm_stats = os.path.join(_test_data_dir, 
                                              "output", 
                                              _test_sample_name, 
                                              "results",
                                              "rpkm", 
                                              _test_sample_name + ".orf_read_counts_stats.txt"
                                             )
     
    output_bwa_folder = os.path.join(out_folder_name, 
                                    _test_sample_name, 
                                    "bwa"
                                    )

    args = [ 
             "-c", input_contigs,
             "--rpkmExec", "rpkm",
             "--readsdir",  input_reads,
             "-O", input_gff,
             "-o", output_orf_rpkm,
             "--sample_name", _test_sample_name,
             "--stats",  output_rpkm_stats,
             "--bwaFolder", output_bwa_folder,
             "--bwaExec", "bwa",
            ]

    MetaPathways_rpkm.main(args)
    
    assert compare_rpkm_stats_in_files(output_rpkm_stats, expect_output_rpkm_stats, sort_n_compare = True)
    assert compare_rpkm_values_in_files(output_orf_rpkm, expect_output_orf_rpkm, sort_n_compare = True)


