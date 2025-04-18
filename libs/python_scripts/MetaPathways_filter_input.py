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
     import os
     import re
     from os import makedirs, sys, remove
     from sys import path
     
     from optparse import OptionParser
     from python_modules.metapaths_utils  import parse_command_line_parameters, fprintf
     from python_modules.sysutil import getstatusoutput, pathDelim
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)

PATHDELIM = pathDelim()



usage= """ Usage: ./MetaPathway_filter_input.py  -i file.fna  --min_length N --log_file logfile.log """ +\
              """ -o outfasta  [ -M map file ]"""

parser = OptionParser(usage)
parser.add_option("-i", "--input_file", dest="input_fasta",
                  help='the input fasta file [REQUIRED]')
parser.add_option("-o", "--output_file", dest="output_fasta",
                  help='the output fasta file [REQUIRED]')
parser.add_option("-m", "--min_length", type='int', dest="min_length",
                  help='minmum length of the sequence allowed')
parser.add_option("-l", "--log_file", dest="log_file",  
                  help='file name to write the statsitics and log into')
parser.add_option("-M", "--map_file", dest="map_file", type='str',  
                  help='file name to store the sequence  name maps')


def valid_arguments(opts, args):
    state = True
    if opts.input_fasta == None :
        print 'Must have an input fasta file'
        state = False

    if opts.output_fasta == None :
        print 'Must have an output fasta file'
        state = False

    if opts.min_length == None :
        print 'Must have a minimum sequence length'
        state = False

    if opts.log_file == None :
        print 'must have a log filename'
        state = False

    return state


def isAminoAcidSequence(sequence):
    if sequence:
        count = 0
        list = [ 'a', 't', 'c', 'g', 'A', 'T', 'C', 'G']
        for x in sequence:
            if x in list:
               count+=1
        if count/len(sequence) < 0.80:
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

def main(argv): 
    (opts, args) = parser.parse_args()

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    min_length = opts.min_length
    inputfile = open(opts.input_fasta,'r')
    outfile = open(opts.output_fasta, 'w') 
    logfile = open(opts.log_file, 'w') 
     

    if opts.map_file:
       mapfile = open(opts.map_file, 'w') 
    else:
       mapfile = None

    
    sample_name = opts.input_fasta;
    sample_name = re.sub(r'^.*/','',sample_name, re.I)
    sample_name = re.sub(r'^.*\\','',sample_name, re.I)
    sample_name = re.sub(r'\.fasta$','',sample_name, re.I)
    sample_name = re.sub(r'\.fna$','',sample_name, re.I)
    sample_name = re.sub(r'\.faa$','',sample_name, re.I)
    sample_name = re.sub(r'\.fas$','',sample_name, re.I)

    BEFORE = 'BEFORE'
    AFTER = 'AFTER'
    NUMSEQ = "Number of sequences :"   
    NUMSEQ_SHORTER = "Number of sequences shorter than "
    AV_LENGTH= "Average length of sequences:"
    MIN_LENGTH= "Minimum length of sequences:"
    MAX_LENGTH= "Maximum length of sequences:" 

    stats = { 
              MIN_LENGTH: { 'BEFORE':10000000, 'AFTER':1000000 },  
              MAX_LENGTH: { 'BEFORE': 0, 'AFTER':0 },  
              NUMSEQ : { 'BEFORE' :0, 'AFTER':0},   
              NUMSEQ_SHORTER : { 'BEFORE':0, 'AFTER':0 },
              AV_LENGTH : { 'BEFORE':0, 'AFTER':0 },
            }  

    length_distribution = {}
    length_cumulative_distribution = {}

    for i in range(0,31):
        length_distribution[i]= 0
        length_cumulative_distribution[i]= 0

    seq_count = 0
    allNames= dict()
    outputStr = ""
    outputLines = []
    for record in read_fasta_records(inputfile):
        seqname = record.name
        seq = record.sequence
        length = len(seq) 
        
        index = int(len(seq) / 50);
        if index >= 30:
            index = 30
    #print length($seq) ."\t".$index."\n";
        length_distribution[index] += 1

        if length < stats[MIN_LENGTH][BEFORE] :
            stats[MIN_LENGTH][BEFORE] = length

        if length > stats[MAX_LENGTH][BEFORE] : 
            stats[MAX_LENGTH][BEFORE] = length

        if length < MIN_LENGTH:
            stats[NUMSEQ_SHORTER][BEFORE] += 1

        stats[AV_LENGTH][BEFORE]  =  stats[AV_LENGTH][BEFORE] + length

        seqvalue = filter_sequence(seq)

        stats[NUMSEQ][BEFORE] += 1
        
        seqlen = len(seqvalue)
        if seqlen>= min_length :
           
           stats[NUMSEQ][AFTER] += 1
           stats[AV_LENGTH][AFTER]  =  stats[AV_LENGTH][AFTER] + seqlen
           if mapfile==None:
              fprintf(outfile, "%s\n", seqname)
           else:
               fprintf(outfile, ">%s\n",  sample_name + '_' + str(seq_count) )
               key = re.sub(r'^>','',seqname)
               fprintf(mapfile, "%s\n", sample_name+ '_' + str(seq_count) + '\t' + key)
               seq_count += 1

           fprintf(outfile, "%s\n",seqvalue)

           if  seqlen < stats[MIN_LENGTH][AFTER] :
               stats[MIN_LENGTH][AFTER] = seqlen
             
           if  seqlen > stats[MAX_LENGTH][AFTER] :
               stats[MAX_LENGTH][AFTER] = seqlen

    if stats[NUMSEQ][BEFORE] > 0 :
      stats[AV_LENGTH][BEFORE]  = stats[AV_LENGTH][BEFORE]/stats[NUMSEQ][BEFORE]
    else:
      stats[AV_LENGTH][BEFORE]  = 0
    if stats[NUMSEQ][AFTER] > 0 :
       stats[AV_LENGTH][AFTER]  = stats[AV_LENGTH][AFTER]/stats[NUMSEQ][AFTER]
    else :
       stats[AV_LENGTH][AFTER]  = 0

    outfile.close()
    inputfile.close()
    if mapfile != None:
       mapfile.close()


    fprintf(logfile, "   %s\n", " \tBEFORE\tAFTER");
    fprintf(logfile, "   %s\n", NUMSEQ +'\t' + str(stats[NUMSEQ][BEFORE]) + '\t' + str(stats[NUMSEQ][AFTER]));
    fprintf(logfile, "   %s\n", NUMSEQ_SHORTER + str(MIN_LENGTH) + ':\t'+ str(stats[NUMSEQ_SHORTER][BEFORE]) + '\t' + str(stats[NUMSEQ_SHORTER][AFTER]))
    fprintf(logfile, "   %s\n", AV_LENGTH +'\t' + str(stats[AV_LENGTH][BEFORE]) + '\t'+ str(stats[AV_LENGTH][AFTER] ))
    fprintf(logfile, "   %s\n", MIN_LENGTH + '\t' + str(stats[MIN_LENGTH][BEFORE]) +'\t'+ str(stats[MIN_LENGTH][AFTER]))
    fprintf(logfile, "   %s\n", MAX_LENGTH +'\t'+ str(stats[MAX_LENGTH][BEFORE]) + '\t' +  str(stats[MAX_LENGTH][AFTER]))

    fprintf(logfile, "\n\n");
    fprintf(logfile, "   READ_LENGTH_RANGE\tFREQUENCY\t\tMIN_LENGTH\tCUMULATIVE_FREQUENCY\n");
    fprintf(logfile, "   -----------------\t---------\t\t----------\t--------------------\n");

    i  = 30
    length_cumulative_distribution[i] = length_cumulative_distribution[i];
    i  -= 1
    while i >= 0:
       length_cumulative_distribution[i] = length_cumulative_distribution[i+1] + length_distribution[i];
       i -= 1

    for i in range(0,31):
       fprintf(logfile, "   %s\n", str(i*50) + '-' + str((i+1)*50) + '\t' +\
                str(length_distribution[i]) +'\t\t\t' + str( (i+1)*50) + '\t' + str(length_cumulative_distribution[i]) )

    logfile.close()


# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

