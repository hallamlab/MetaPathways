#!/usr/bin/env python
# file parse.py: parsers for map file, distance matrix file, env file

__author__ = "Kishori Mohan Konwar"
__copyright__ = "MetaPathways"
__credits__ = ["Kishori M Konwar"]
__license__ = "L"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

from collections import defaultdict
from copy import deepcopy
import os, re

from metapathways.utils.metapathways_utils import (
    parse_command_line_parameters,
    eprintf,
    halt_process,
)
from metapathways.utils.utils import *


class MetaPathwaysError(Exception):
    pass


def parse_mapping_file(lines, strip_quotes=True, suppress_stripping=False):
    """Parser for map file that relates samples to metadata.

    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    Result: list of lists of fields, incl. headers.
    """
    if hasattr(lines, "upper"):
        # Try opening if a string was passed
        try:
            lines = open(lines, "U")
        except IOError:
            raise (
                "A string was passed that doesn't refer " "to an accessible filepath."
            )

    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"', "")
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"', "").strip()
    else:
        if suppress_stripping:
            # don't remove quotes or spaces
            strip_f = lambda x: x
        else:
            # remove spaces but not quotes
            strip_f = lambda x: x.strip()

    # Create lists to store the results
    mapping_data = []
    header = []
    comments = []

    # Begin iterating over lines
    for line in lines:
        line = strip_f(line)
        if not line or (suppress_stripping and not line.strip()):
            # skip blank lines when not stripping lines
            continue

        if line.startswith("#"):
            line = line[1:]
            if not header:
                header = line.strip().split("\t")
            else:
                comments.append(line)
        else:
            mapping_data.append(map(strip_f, line.split("\t")))
    if not header:
        raise ("No header line was found in mapping file.")
    if not mapping_data:
        raise ("No data found in mapping file.")

    return mapping_data, header, comments


def parse_mapping_file_to_dict(*args, **kwargs):
    """Parser for map file that relates samples to metadata.

    input format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    calls parse_mapping_file, then processes the result into a 2d dict, assuming
    the first field is the sample id
    e.g.: {'sample1':{'age':'3','sex':'male'},'sample2':...

    returns the dict, and a list of comment lines"""
    mapping_data, header, comments = parse_mapping_file(*args, **kwargs)
    return mapping_file_to_dict(mapping_data, header), comments


def mapping_file_to_dict(mapping_data, header):
    """processes mapping data in list of lists format into a 2 deep dict"""
    map_dict = {}
    for i in range(len(mapping_data)):
        sam = mapping_data[i]
        map_dict[sam[0]] = {}
        for j in range(len(header)):
            if j == 0:
                continue  # sampleID field
            map_dict[sam[0]][header[j]] = sam[j]
    return Dict2D(map_dict)


def parse_prefs_file(prefs_string):
    """Returns prefs dict evaluated from prefs_string.

    prefs_string: read buffer from prefs file or string containing prefs
        dict.  Must be able to evauluated as a dict using eval.
    """
    try:
        prefs = dict(eval(prefs_string))
    except TypeError:
        raise ("Invalid prefs file. Prefs file must contain a valid prefs dictionary.")
    return prefs


def group_by_field(table, name):
    """Returns dict of field_state:[row_headers] from table.

    Use to extract info from table based on a single field.
    """
    try:
        col_index = table[0].index(name)
    except ValueError:
        raise (alueError, "Couldn't find name %s in headers: %s" % (name, table[0]))
    result = defaultdict(list)
    for row in table[1:]:
        header, state = row[0], row[col_index]
        result[state].append(header)
    return result


def group_by_fields(table, names):
    """Returns dict of (field_states):[row_headers] from table.

    Use to extract info from table based on combinations of fields.
    """
    col_indices = map(table[0].index, names)
    result = defaultdict(list)
    for row in table[1:]:
        header = row[0]
        states = tuple([row[i] for i in col_indices])
        result[states].append(header)
    return result


def parse_distmat_to_dict(table):
    """Parse a dist matrix into an 2d dict indexed by sample ids.

    table: table as lines
    """

    col_headers, row_headers, data = parse_matrix(table)
    assert col_headers == row_headers

    result = defaultdict(dict)
    for (sample_id_x, row) in zip(col_headers, data):
        for (sample_id_y, value) in zip(row_headers, row):
            result[sample_id_x][sample_id_y] = value
    return result


def parse_bootstrap_support(lines):
    """Parser for a bootstrap/jackknife support in tab delimited text"""
    bootstraps = {}
    for line in lines:
        if line[0] == "#":
            continue
        wordlist = line.strip().split()
        bootstraps[wordlist[0]] = float(wordlist[1])

    return bootstraps


def fields_to_dict(lines, delim="\t"):
    """makes a dict where first field is key, rest are vals."""
    result = {}
    for line in lines:
        # skip empty lines
        if strip_f:
            fields = map(strip_f, line.split(delim))
        else:
            fields = line.split(delim)
        if not fields[0]:  # empty string in first field implies problem
            continue
        result[fields[0]] = fields[1:]
    return result


def parse_metapaths_parameters(filename):
    """Return 2D dict of params (and values, if applicable) which should be on"""
    # The qiime_config object is a default dict: if keys are not
    # present, {} is returned
    def return_empty_dict():
        return dict()

    with open(filename, "r") as filep:
        result = {}
    
        lines = filep.readlines()
        for line in lines:
            line = line.strip()
            if line and not line.startswith("#"):
                fields = line.split()
                try:
                    script_id, parameter_id = fields[0].split(":")
                    value = ",".join([x.strip() for x in fields[1:]])
                    value = re.sub(",,", ",", value)
                    # if value.upper() == 'FALSE' or value.upper() == 'NONE':
                    #    continue
                    # elif value.upper() == 'TRUE':
                    #    value = None
                    # else:
                    #    pass
                    if script_id not in result:
                        result[script_id] = {}
                    result[script_id][parameter_id] = value
                except KeyError:
                    result[script_id] = {parameter_id: value}
    # result['filename'] = filename
    return result


def parse_parameter_file(filename):
    """Return 2D dict of params (and values, if applicable) which should be on"""
    # The qiime_config object is a default dict: if keys are not
    # present, {} is returned
    def return_empty_dict():
        return dict()

    result = defaultdict(return_empty_dict)
    file = open(filename, "r")
    lines = file.readlines()
    file.close()

    for line in lines:
        line = line.strip()
        if line and not line.startswith("#"):
            fields = line.split()
            script_id, parameter_id = fields[0].split(":")
            try:
                value = " ".join([x.strip() for x in fields[1:]])
            except IndexError:
                continue

            if value.upper() == "FALSE" or value.upper() == "NONE":
                continue
            elif value.upper() == "TRUE":
                value = None
            else:
                pass

            try:
                result[script_id][parameter_id] = value
            except KeyError:
                result[script_id] = {parameter_id: value}
    return result
