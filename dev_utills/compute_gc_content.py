#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import os, re
     from os import sys
     from sys import path
     from optparse import OptionParser

except:
     print """ Could not load some user defined  module functions"""
     sys.exit(3)


usage= sys.argv[0] + """ -i <input file> -o <output file>"""

parser = None
def createParser():
    global parser
    epilog = """This script computes the GC content of individual sequences  in the input file and writes 
the results into an output file.
"""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-i", "--input_file", dest="input_fasta",
                      help='the input fasta filename [REQUIRED]')
    parser.add_option("-o", "--output_file", dest="output",
                      help='the output file name [REQUIRED]')
    

def valid_arguments(opts, args):
    state = True
    if opts.input_fasta == None :
        print 'ERROR: Missing input fasta file'
        state = False

    if opts.output == None :
        print 'ERROR: Missing output file name'
        state = False
    return state


class FastaRecord(object):
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

#    return FastaRecord(title, sequence)

def read_fasta_records(input_file):
    records = []
    sequence=""
    name=""
    while 1:
         line = input_file.readline()
         if line == "": 
            if sequence!="" and name!="":
               records.append(FastaRecord(name, sequence))
            return  records

         if line=='\n':
            continue

         line = line.rstrip()
         if  line.startswith(">") :
            if sequence!="" and name!="":
               records.append(FastaRecord( re.sub('>', '', name), sequence))

            name = line.rstrip()
            sequence =""
         else:
            sequence = sequence + line.rstrip()
    return records

        
"""file name:gc_content.py"""
def validate_base_sequence(base_sequence, Nucltype='DNA'):
    """ return true if no other letters"""
    seq=base_sequence.upper()
    return len(seq)==(seq.count('T' if Nucltype=='DNA' else 'U')+seq.count('A')+seq.count('C')+seq.count('G')+seq.count('N'))

def gc_content(base_sequence,Nucltype='DNA'):
    """return the GC percentage in base_sequence"""
    assert validate_base_sequence(base_sequence,Nucltype), 'sequence has invalid characters'
    seq= base_sequence.upper()
    return ((seq.count('G')+seq.count('C'))/len(seq))

#print(gc_content(seq,Nucltype))

def fprintf(file, fmt, *args):
    file.write(fmt % args)

# the main function
def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    # read the fasta records
    with open(opts.input_fasta, 'r') as input_file:
      with open(opts.output, 'w') as output_file:
         records = read_fasta_records(input_file)
         for record in records:
            fprintf(output_file, "%s\t%.4f\n", record.name, gc_content(record.sequence, Nucltype='DNA') )
       


# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

