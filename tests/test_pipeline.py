import pytest
import sys
import os
import filecmp
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
  

def compare_lines_in_files(file1, file2, sort_n_compare):
    with open(file1, 'r') as fout, open(file2, 'r') as expfout:
       if sort_n_compare:
           output_lines = sorted(fout.readlines())
           expected_output_lines = sorted(expfout.readlines())
       else:
           output_lines = fout.readlines()
           expected_output_lines = expfout.readlines()
     
    #print([a == b for a, b in zip(output_lines, expected_output_lines)])
    return all([a == b for a, b in zip(output_lines, expected_output_lines)])


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
_test_input_fasta_file = os.path.join(_test_data_dir, "lagoon-sample/input/lagoon-sample.fasta")
_test_ref_data_dir = os.path.join(_test_data_dir, "ref_data")

def xest_pipeline(tmpdir):
    from metapathways import pipeline 
    _test_sample_name = "lagoon-sample"
    _test_input_fasta_file = os.path.join(_test_data_dir, _test_sample_name, "/input/", _test_sample_name + ".fasta")

    commands = ["MetaPathways", "-i",  _test_data_dir + '/lagoon-sample/input/',
                "-o", str(tmpdir),
                "-s", "lagoon-sample",
                "-p", _test_data_dir + "/template_param.txt",
                "-d", _test_ref_data_dir
               ]
    #print( pipeline.__dict__['__file__'])
    pipeline.process(commands)
    assert 1 == 1

                           
def test_filter_nuc_input(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_filter_input
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"
    _test_input_fasta_file = os.path.join(_test_data_dir, 
                                          _test_sample_name, 
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

    expected_output_fasta = os.path.join(_test_data_dir, 
                                         _test_sample_name, 
                                         "output", 
                                         _test_sample_name, 
                                         'preprocessed', 
                                         _test_sample_name + ".fasta"
                                        )

    expected_output_mapping_file = os.path.join(_test_data_dir, 
                                         _test_sample_name, 
                                         "output", _test_sample_name, 
                                         "preprocessed", 
                                         _test_sample_name + ".mapping.txt"
                                        )

    expected_output_contig_length = os.path.join(_test_data_dir, 
                                                 _test_sample_name, 
                                                 'output', 
                                                 _test_sample_name, 
                                                 "run_statistics", 
                                                 _test_sample_name + ".contig.lengths.txt"
                                               )
    expected_output_logfile = os.path.join(out_folder_name, 
                            _test_sample_name, 
                            "run_statistics", 
                            _test_sample_name + ".nuc.stats"
                           )

    MetaPathways_filter_input.main(args)
    assert filecmp.cmp(output_fasta, expected_output_fasta)
    assert filecmp.cmp(output_mapping_file, expected_output_mapping_file)
    assert filecmp.cmp(output_contig_length, expected_output_contig_length)
    assert filecmp.cmp(output_logfile, expected_output_logfile)


def test_orf_predictions_input(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_orf_prediction
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    prod_input = os.path.join(_test_data_dir, 
                              _test_sample_name, 
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

    expected_output_gff = os.path.join(_test_data_dir,
                                       _test_sample_name, 
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
    assert filecmp.cmp(output_gff, expected_output_gff)


def test_create_amino_sequences(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_create_amino_sequences 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    g_input = os.path.join(_test_data_dir, 
                              _test_sample_name, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".gff"
                            )

    n_input = os.path.join(_test_data_dir, 
                              _test_sample_name, 
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

    expected_output_amino = os.path.join(_test_data_dir, 
                                       _test_sample_name, 
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

    expected_output_nuc = os.path.join(_test_data_dir, 
                                       _test_sample_name, 
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

    expected_output_gff = os.path.join(_test_data_dir,
                                       _test_sample_name, 
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
    assert compare_lines_in_files(output_gff, expected_output_gff, sort_n_compare = True)
    assert compare_fasta_files(output_amino, expected_output_amino, sort_n_compare = True)
    assert compare_fasta_files(output_nuc, expected_output_nuc, sort_n_compare = True)


def test_refscores(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_refscore 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

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

    expected_output_refscores = os.path.join(_test_data_dir, 
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

    assert compare_lines_in_files(output_refscores, expected_output_refscores, sort_n_compare = True)


def test_filter_amino_input(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_filter_input
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    input_faa = os.path.join( _test_data_dir, 
                              _test_sample_name, 
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

    expected_output_qced_faa = os.path.join(_test_data_dir, 
                                         _test_sample_name, 
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

    expected_output_aminoseq_length = os.path.join(_test_data_dir, 
                                                   _test_sample_name, 
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

    expected_output_logfile = os.path.join(out_folder_name, 
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

    assert compare_fasta_files(output_qced_faa, expected_output_qced_faa, sort_n_compare = True)
    assert compare_lines_in_files(output_aminoseq_length, expected_output_aminoseq_length, sort_n_compare = False)
    assert compare_lines_in_files(output_logfile, expected_output_logfile, sort_n_compare = False)

def test_func_search(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_func_search 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    faa_input = os.path.join(_test_data_dir, 
                              _test_sample_name, 
                              "output",
                              _test_sample_name, 
                              "orf_prediction", 
                              _test_sample_name + ".qced.faa"
                            )

    ref_blast_db = os.path.join(_test_data_dir, 
                            "ref_data",
                            "functional",
                            "formatted",
                            "refseq-small-sample"
                           )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "orf_prediction"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)

    output_blast_results = os.path.join(out_folder_name, 
                              _test_sample_name, 
                              "blast_results", 
                              _test_sample_name + ".refseq-small-sample.BLASTout"
                             )

    expected_output_blast_results = os.path.join(_test_data_dir, 
                                             _test_sample_name, 
                                             "output", 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + ".refseq-small-sample.BLASTout"
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

    assert compare_lines_in_files(output_blast_results, expected_output_blast_results, sort_n_compare = True)

def test_refscores(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_refscore 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

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

    expected_output_refscores = os.path.join(_test_data_dir, 
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

    assert compare_lines_in_files(output_refscores, expected_output_refscores, sort_n_compare = True)

def test_parse_blast(tmpdir):
    from metapathways import MetaPathways_parse_blast 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    ref_blast_db_annots = os.path.join(_test_data_dir, 
                            "ref_data",
                            "functional",
                            "formatted",
                            "refseq-small-sample-names.txt"
                           )

    input_blast_results = os.path.join(_test_data_dir, 
                                        _test_sample_name, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                        _test_sample_name + ".refseq-small-sample.BLASTout"
                                      )

    input_blast_refscores = os.path.join(_test_data_dir, 
                                         _test_sample_name, 
                                         'output',
                                         _test_sample_name, 
                                         "blast_results", 
                                         _test_sample_name + ".refscores.BLAST"
                                        )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)

    output_parsed_blast_results = os.path.join(out_folder_name, 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + ".refseq-small-sample.BLASTout.parsed.txt"
                                            )

    expected_output_parsed_blast_results = os.path.join(_test_data_dir, 
                                             _test_sample_name, 
                                             "output", 
                                             _test_sample_name, 
                                             "blast_results", 
                                             _test_sample_name + ".refseq-small-sample.BLASTout.parsed.txt"
                                            )
    args = [ 
             "-b", input_blast_results,
             "-r", input_blast_refscores,
             "-m", ref_blast_db_annots,
             "-d", "refseq-small-sample",
             "-o", output_parsed_blast_results,
             "--min_bsr", "0.4",
             "--min_score", "20",
             "--min_length", "45",
             "--max_evalue", "0.000001",
             "--algorithm", "BLAST"
           ]

    MetaPathways_parse_blast.main(args)
    assert compare_lines_in_files(output_parsed_blast_results, expected_output_parsed_blast_results, sort_n_compare = True)

def test_rRNA_stats_calculator(tmpdir):
    from metapathways import MetaPathways_rRNA_stats_calculator 
    from metapathways import sysutil as sysutils
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    input_blast_results = os.path.join(_test_data_dir, 
                                        _test_sample_name, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                        _test_sample_name + ".rRNA.silva-lsu-small-sample.BLASTout"
                                      )

    ref_blast_db = os.path.join(_test_data_dir, 
                            "ref_data",
                            "taxonomic",
                            "silva-lsu-small-sample"
                           )

    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "blast_results"), exist_ok = True)
    os.makedirs(os.path.join(out_folder_name, _test_sample_name, "results",  "rRNA"), exist_ok = True)

    output_lsu_rRNA_stats = os.path.join(out_folder_name, 
                                             _test_sample_name, 
                                             "results",
                                             "rRNA", 
                                             _test_sample_name + ".silva-lsu-small-sample.rRNA.stats.txt"
                                            )

    expected_output_lsu_rRNA_stats = os.path.join(_test_data_dir, 
                                             _test_sample_name, 
                                             "output", 
                                             _test_sample_name, 
                                             "results",
                                             "rRNA", 
                                             _test_sample_name + ".silva-lsu-small-sample.rRNA.stats.txt"
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
    assert compare_lines_in_files(output_lsu_rRNA_stats, expected_output_lsu_rRNA_stats, sort_n_compare = True)

    """ Test the BLAST output """
    ref_blast_db_formatted = os.path.join(_test_data_dir, 
                                          "ref_data",
                                          "taxonomic",
                                          "formatted",
                                          "silva-lsu-small-sample"
                                          )

    input_contigs = os.path.join(_test_data_dir, 
                                         _test_sample_name, 
                                         "output", 
                                         _test_sample_name, 
                                         'preprocessed', 
                                         _test_sample_name + ".fasta"
                                        )


    expected_output_blast_results = os.path.join(_test_data_dir, 
                                        _test_sample_name, 
                                       "output", 
                                       _test_sample_name, 
                                       "blast_results", 
                                       _test_sample_name + ".rRNA.silva-lsu-small-sample.BLASTout"
                                      )

    output_blast_results = os.path.join(out_folder_name, 
                                        _test_sample_name, 
                                        "blast_results",
                                         _test_sample_name + ".rRNA.silva-lsu-small-sample.BLASTout"
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
    assert compare_lines_in_files(output_blast_results, expected_output_blast_results, sort_n_compare = True)

def test_tRNA_scan(tmpdir):
    from metapathways import MetaPathways_tRNA_scan 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    input_contigs = os.path.join(_test_data_dir, 
                                        _test_sample_name, 
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

    expected_output_tRNA_stats = os.path.join(_test_data_dir, 
                                              _test_sample_name, 
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

    expected_output_tRNA_fasta = os.path.join(_test_data_dir, 
                                             _test_sample_name, 
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
    assert compare_lines_in_files(output_tRNA_stats, expected_output_tRNA_stats, sort_n_compare = True)
    
    assert compare_fasta_files(output_tRNA_fasta, expected_output_tRNA_fasta, sort_n_compare = True)

def test_rpkm(tmpdir):
    from metapathways import MetaPathways_rpkm 
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"

    input_contigs = os.path.join(_test_data_dir, 
                                        _test_sample_name, 
                                       "output", 
                                       _test_sample_name, 
                                       "preprocessed", 
                                        _test_sample_name + ".fasta"
                                      )
    input_reads = os.path.join(_test_data_dir, 
                               _test_sample_name, 
                               "input", 
                               "reads" 
                              )


    input_gff = os.path.join(_test_data_dir, 
                             _test_sample_name, 
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

    expected_output_orf_rpkm = os.path.join(_test_data_dir, 
                                     _test_sample_name, 
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
                                     _test_sample_name + ".orf_rpkm.txt"
                                    )

    expected_output_rpkm_stats = os.path.join(_test_data_dir, 
                                              _test_sample_name, 
                                              "output", 
                                              _test_sample_name, 
                                              "results",
                                              "rpkm", 
                                              _test_sample_name + ".orf_rpkm.txt"
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
             "--sample_name", 'lagoon-sample',
             "--stats",  output_rpkm_stats,
             "--bwaFolder", output_bwa_folder,
             "--bwaExec", "bwa",
            ]

    MetaPathways_rpkm.main(args)
   
    assert compare_lines_in_files(output_rpkm_stats, expected_output_rpkm_stats, sort_n_compare = True)
    assert compare_lines_in_files(output_rpkm_stats, expected_output_rpkm_stats, sort_n_compare = True)
