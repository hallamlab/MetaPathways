import pytest
import sys
import os
import filecmp

_test_data_dir = os.path.join(os.path.split(__file__)[0], "data/")
_test_input_fasta_file = os.path.join(_test_data_dir, "lagoon-sample/input/lagoon-sample.fasta")
_test_ref_data_dir = os.path.join(_test_data_dir, "ref_data")

def test_pipeline(tmpdir):
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

                           
def test_filter_input(tmpdir):
    from metapathways import pipeline 
    from metapathways import MetaPathways_filter_input
    out_folder_name = str(tmpdir)
    _test_sample_name = "lagoon-sample"
    _test_input_fasta_file = '/'.join([_test_data_dir, _test_sample_name, "input", _test_sample_name + ".fasta"])

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

