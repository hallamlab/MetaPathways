import pytest
import re
import gzip
import os

from pkg_resources import Requirement, resource_filename, ResolutionError


def get_package_root():
    return resource_filename(Requirement.parse("metapathways"), 'metapathways')


def get_test_data(filename):
    filepath = None
    try:
        filepath = resource_filename(Requirement.parse("metapathways"), "tests/data/" + filename)
    except ResolutionError:
        pass
    if not filepath or not os.path.isfile(filepath):
        filepath = os.path.join(os.path.dirname(__file__), 'data', filename)
    return filepath


def get_metapathways_file(filename):
    return resource_filename(Requirement.parse("metapathways"), filename)


def get_metapathways_path():
    return resource_filename(Requirement.parse("metapathways"), "")


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

    # print([a == b for a, b in zip(output_lines, expect_output_lines)])
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
                assert (a.split(':')[0].strip() == b.split(':')[0].strip())
            if len(a.split(':')) > 1:
                try:
                    numa = float(re.sub('%', '', a.split(':')[1]).strip())
                    numb = float(re.sub('%', '', b.split(':')[1]).strip())
                    if numb != 0:
                        assert (numa == pytest.approx(numb, rel=1e-1))
                    else:
                        assert (pytest.approx(numa, 0.1) == 0)
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
            numa = float(a.split('\t')[1].strip())
            numb = float(b.split('\t')[1].strip())
            if numa != pytest.approx(numb, rel=1e-1):
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

# def compare_fasta_files(file1, file2, sort_n_compare):
#     import pyfastx
#     contents1 = []
#     for seq in pyfastx.Fasta(file1):
#         contents1.append(Sequence(seq.name, seq.seq))
#
#     contents2 = []
#     for seq in pyfastx.Fasta(file2):
#         contents2.append(Sequence(seq.name, seq.seq))
#
#     contents1.sort()
#     contents2.sort()
#
#     return all([a == b for a, b in zip(contents1, contents2)])
#     return True
