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
     from os import makedirs, sys, remove, rename
     from sys import path
     from optparse import OptionParser

     from libs.python_modules.utils.metapathways_utils  import parse_command_line_parameters, fprintf , printf
     from libs.python_modules.utils.sysutil import getstatusoutput, pathDelim
     from libs.python_modules.parsers.fastareader  import FastaReader
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()



usage= sys.argv[0] + """ -s <sample_name>  -f <output_folder> -i <input_folder> """

parser = None
def createParser():
    global parser
    epilog = """
This script generates run stats of a MetaPathways processed sample."""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-s", "--sample_name", dest="sample_name",
                      help='the sample name [REQUIRED]')

    parser.add_option("-f", "--output_folder", dest="output_folder",
                      help='the output folder where the sample is [REQUIRED]')

    parser.add_option("-i", "--input_folder", dest="input_folder",
                      help='the input folder where the input is [REQUIRED]')


def valid_arguments(opts, args):
    state = True

    if opts.sample_name == None :
        print 'ERROR: Missing sample name'
        state = False

    if opts.output_folder == None :
        print 'ERROR: Missing output folder for sample'
        state = False

    if opts.input_folder == None :
        print 'WARNING: Missing input folder for sample'

    return state




def isAminoAcidSequence(sequence):
    if sequence:
        count = 0 

        list= {
           'A':  0,  'R':  0,  'N':  0,  'D':  0,  'C':  0,  'Q':  0,  'E':  0,  'G':  0,  
           'H':  0,  'I':  0,  'L':  0,  'K':  0,  'M':  0,  'F':  0,  'P':  0, 'S':  0,  
           'T':  0, 'W':  0, 'Y':  0, 'V':  0,  'B':  0,  'J':  0,  'Z':  0,  }

        for x in sequence:
            if x.upper() in list:
               list[x.upper()]=1

        count = 0 
        for x in list:
           count += list[x]

        if count > 10: 
            return True
        else:
             return False
    return True
       

def filter_sequence(sequence):
   if isAminoAcidSequence(sequence):
       return sequence

   sequence = re.sub(r'[^atcgATCG]','-', sequence.strip())
   subsequences =  sequence.split('-')

   max_length = 0;
   longest_sequence = ""; 
   for seq  in subsequences: 
      if len(seq) > max_length :
          longest_sequence = seq
          max_length = len(seq)

   return  longest_sequence



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
               records.append(FastaRecord(name, sequence))

            name = line.rstrip()
            sequence =""
         else:
            sequence = sequence + line.rstrip()
    return records

        

# the main function
SIZE = 1000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    sample_name = opts.sample_name
    output_folder = opts.output_folder
    input_folder = opts.input_folder


    _MAX = 1000000000000

    seq_count = 0
    allNames= dict()
    outputStr = ""

    if  opts.input_folder!=None:
       formats= "%s\tNumber of sequences in input file BEFORE QC (nucleotide)\t%s\n" 
       compute_fasta_stats(formats, input_folder  + PATHDELIM + sample_name + ".fasta", 'nucleotide', 995)

    formats= "%s\tNumber of sequences AFTER QC (nucleotide)\t%s\n"
    compute_fasta_stats(formats, output_folder  + PATHDELIM + sample_name + PATHDELIM + 'preprocessed' + PATHDELIM +  sample_name + ".fasta", 'nucleotide', 1000)

    formats= "%s\tNumber of translated ORFs BEFORE QC (amino)\t%s\n" 
    compute_fasta_stats(formats, output_folder  + PATHDELIM + sample_name + PATHDELIM + 'orf_prediction' + PATHDELIM +  sample_name + ".faa", 'amino', 1995)

    formats= "%s\tNumber of translated ORFs AFTER QC (amino)\t%s\n" 
    compute_fasta_stats(formats, output_folder  + PATHDELIM + sample_name + PATHDELIM + 'orf_prediction' + PATHDELIM +  sample_name + ".qced.faa", 'amino', 2000)



def compute_fasta_stats(formats, input_file, seqtype, priority):
    MIN_LENGTH='MIN_LENGTH'
    MAX_LENGTH='MAX_LENGTH'
    NUMSEQ='NUMSEQ' 
    TOTAL_LENGTH='TOTAL_LENGTH'
    AVG_LENGTH='AVG_LENGTH' 

    stats = { 
              MIN_LENGTH: 0,  
              MAX_LENGTH: 0,  
              NUMSEQ : 0,    
              TOTAL_LENGTH: 0,
              AVG_LENGTH : 0 
    }  

    """ min length """
    _MAX = 1000000000000
    stats[MAX_LENGTH] = -(_MAX)
    stats[MIN_LENGTH]= _MAX

    fastareader= FastaReader(input_file)

    """ process one fasta sequence at a time """
    lengths_str=""
    for record in fastareader:
        seqname = record.name
        seq = record.sequence
        length = len(seq)
        
        stats[NUMSEQ] += 1
        
        stats[AVG_LENGTH]  =  stats[AVG_LENGTH] + length

        if stats[MIN_LENGTH] > length:
           stats[MIN_LENGTH] = length

        if stats[MAX_LENGTH] < length:
           stats[MAX_LENGTH] = length



    
    if stats[NUMSEQ] > 0 :
      stats[AVG_LENGTH]  = stats[AVG_LENGTH]/stats[NUMSEQ]
    else:
      stats[AVG_LENGTH]  = 0



    #     printf("%s\tNumber of sequences in input file BEFORE QC (%s)\t%s\n" %(str(priority), opts.seqtype,  str(stats[NUMSEQ][BEFORE])) )

    #     printf("%s\tNumber of sequences AFTER QC (%s)\t%s\n" %(str(priority + 5), opts.seqtype, str(stats[NUMSEQ][AFTER])))
    printf(formats  %(str(priority + 5), str(stats[NUMSEQ])))
    printf("%s\t-min length\t%s\n" %(str(priority + 6), str(stats[MIN_LENGTH])) )
    printf("%s\t-avg length\t%s\n" %( str(priority + 7), str(int(stats[AVG_LENGTH]))))
    printf("%s\t-max length\t%s\n" %( str(priority + 8), str(stats[MAX_LENGTH])) )
    printf("%s\t-total base pairs (bp)\t%s\n" %( str(priority + 9), str(int(stats[AVG_LENGTH]* stats[NUMSEQ])) ))
    #     printf("%s\tNumber of translated ORFs BEFORE QC (%s)\t%s\n" %(str(priority), opts.seqtype,  str(stats[NUMSEQ][BEFORE])) )

#         printf("%s\tNumber of tranlated ORFs AFTER QC (%s)\t%s\n" %(str(priority + 5), opts.seqtype, str(stats[NUMSEQ][AFTER])))



def MetaPathways_filter_input(argv, errorlogger = None, runstatslogger = None):
    createParser()
    main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger) 
    return (0,'')

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

