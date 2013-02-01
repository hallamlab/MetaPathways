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
     from os import makedirs, sys, remove
     from sys import path
     
     from optparse import OptionParser
     from python_modules.metapaths_utils  import parse_command_line_parameters, fprintf
     from python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     sys.exit(3)






usage= """./run_pgdb_pipeline.py -i input_fasta_file -o output_file [-B blast_executable -F formatdb_executable] """
parser = OptionParser(usage)
parser.add_option("-i", "--input_file", dest="input_fasta",
                  help='the input fasta file [REQUIRED]')
parser.add_option("-o", "--output_file", dest="output_file",
                  help='the output fasta file [REQUIRED]')
parser.add_option("-B", "--BLAST_EXEUTABLE", dest="blast_executable",
                  help='the BLAST executable  [REQUIRED]')
parser.add_option("-F", "--FORMAT_EXECUTABLE", dest="formatdb_executable",
                  help='the FORMATDB executable file [REQUIRED]')


def check_arguments(opts, args):
    if opts.input_fasta == None or opts.output_file == None:
       return True
    else:
       return False

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

def format_db(formatdb_executable, seq_subset_file):
    cmd='%s -p T -i %s' %(formatdb_executable, seq_subset_file.name)
    result= getstatusoutput(cmd)
    

def blast_against_itself(blast_executable, seq_subset_file, blast_table_out):
    cmd='%s -outfmt 6 -db  %s -query %s -out  %s' %(blast_executable,  seq_subset_file.name, seq_subset_file.name, blast_table_out)
    result= getstatusoutput(cmd)

def add_refscore_to_file(blast_table_out, refscore_file, allNames):
    infile = open( blast_table_out,'r')

    refscores = {}

    lines = infile.readlines()
    for line in lines:
       line=line.rstrip()
       fields = line.split('\t')
       if len(fields) != 12:
          print 'Error in the blastout file'
          sys.exit(1)
       if fields[0].rstrip()==fields[1].rstrip():
      #    fprintf(refscore_file, "%s\t%s\n",fields[0], fields[11])
          if fields[0] in refscores and  refscores[fields[0]] < fields[11]:
             refscores[fields[0]]=fields[11]
          else:
             refscores[fields[0]]=fields[11]

    for key, value in refscores.iteritems():
       allNames[key] = True
       fprintf(refscore_file, "%s\t%s\n",key, value)

    infile.close()
        

# compute the refscores
def compute_refscores(formatdb_executable, blast_executable,seq_subset_file, refscore_file, allNames):
    format_db(formatdb_executable, seq_subset_file)
    blast_table_out = seq_subset_file.name + ".blastout"
    blast_against_itself(blast_executable, seq_subset_file, blast_table_out)
    add_refscore_to_file(blast_table_out,refscore_file, allNames)
    return None

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if check_arguments(opts, args):
       print usage
       sys.exit(0)

    input_fasta = opts.input_fasta
    output_file = opts.output_file
    blast_executable = opts.blast_executable
    formatdb_executable = opts.formatdb_executable
 
    # input file to blast with itself to commpute refscore
    infile = open(input_fasta,'r')
   
    #this file has the refscores of the entire file
    outfile = open(output_file, 'w') 

    count = 0

    allNames= dict()
    for record in read_fasta_records(infile):
        if count %20 == 0:
            if count > 0:
              seq_subset_file.close()
              compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames);

              # now remove the old file
              remove(seq_subset_file.name)
              remove(seq_subset_file.name +".blastout")
              remove(seq_subset_file.name +".phr")
              remove(seq_subset_file.name +".pin")
              remove(seq_subset_file.name +".psq")

            seq_subset_file = open(output_file +'.tmp.'+ str(count) +'.fasta','w')
        allNames[record.name.replace(">","")] = False;    
        fprintf(seq_subset_file, "%s\n", record.name)
        fprintf(seq_subset_file, "%s\n", record.sequence)

        count = count + 1

    #print str(count) + "   "  + "going to blast last sequence "
    if (count) %20 != 0:
       #print str(count) + "   "  + "last sequence "
       seq_subset_file.close()
       compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames);
       remove(seq_subset_file.name)
       remove(seq_subset_file.name +".blastout")
       remove(seq_subset_file.name +".phr")
       remove(seq_subset_file.name +".pin")
       remove(seq_subset_file.name +".psq")

    #print count
    for key in allNames:
        if allNames[key] ==False:
           fprintf(outfile, "%s\t%s\n",key, 1000000)

    outfile.close()

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

