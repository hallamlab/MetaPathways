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
     import sys, os, re, math
     from glob import glob
     from os import makedirs, sys, remove, rename, path
     from optparse import OptionParser

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     sys.exit(3)



DBLIST = ["COG-14-2016-10-20",  "kegg-uniprot-2016-10-20", "metacyc-2016-10-31", "refseq-2016-10-06-rel-78", "eggnog-v4-2016-10-30"]

def printf(fmt, *args):
   sys.stdout.write(fmt % args)
   sys.stdout.flush()

def eprintf(fmt, *args):
   sys.stderr.write(fmt % args)
   sys.stderr.flush()


class FastaRecord():
    def __init__(self, longname, sequence):
      self.longname = longname
      self.sequence = sequence
      fields = [ x.strip() for x in self.longname.split(' ') ]
      if len(fields) > 0:
         self.name = fields[0]
      else:
         self.name = None


class FastaReader():
    """Parses a fasta record from a string or file."""
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""


    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self

    def close(self):
         self.file.close()

    def next(self):
        if self.stop:
          raise StopIteration
        
        try:
           if not self.name: 
               self.name = self.file.readline().strip()
           line = self.file.readline()
        except:
           line = None

        if not line:
           self.stop = True
           raise StopIteration

        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip()) 
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = self.future_name

        if line:
          self.future_name = line.strip()

        self.sequence =''.join(fragments)
        self.seqname = self.name
        
        return FastaRecord(self.name, self.sequence)




usage= sys.argv[0] + """ -i file.fna """

parser = None
def createParser():
    global parser
    epilog = """
This script computes the sequence stats for the fasta files
"""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-f",  dest="folders", action='append',
                      help='add the folder to be examined')

    parser.add_option("-s",  dest="stages",   default=[], action='append', 
                      help=''' INPUT   : 1\n
                               ORFs    : 2\n
                               B/LAST  : 3\n
                               PARSE   : 4\n
                               ANNOT   : 5\n
                               PGDB    : 6\n
                               add the folder to be examined''')

    parser.add_option("-t",  dest="type",   default='1', choices=['1', '2', '3', '4'], 
                      help=''' present    : 1 
                               isNonEmpty : 2
                               maxSize    : 3
                               Size       : 4
                               turns on the cumulative mod''')
    parser.add_option("-c", action="store_false", dest="cumul", default=True,
                       help="print the preceeding stages")


def valid_arguments(opts, args):
    state = True
    if opts.folders == None :
        print 'ERROR: Did not specify any folder'
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

       
def numLines(filename):
    count = 0
    if path.exists(filename):
       for line in open(filename):
          count += 1
        

    return count

def maxSizeFasta(file):
    """ process one fasta sequence at a time """

    fastareader= FastaReader(file)
    max_length=0
    count = 0

    for record in fastareader:
        if count > 10000:
          break
        seqname = record.name
        seq = record.sequence
        length = len(seq)

        count += 1

        if length > max_length:
          max_length =length

    fastareader.close()
    return max_length

def avgSizeFasta(file):
    """ process one fasta sequence at a time """

    fastareader= FastaReader(file)
    tot_length=0
    count = 0
    for record in fastareader:
        if count > 10000:
          break
        seqname = record.name
        seq = record.sequence
        length = len(seq)
        tot_length += length
        count += 1

    fastareader.close()
    avg = tot_length/count
    return avg


def extractSampleName(name):
    sample_name = name
    sample_name = re.sub(r'^.*/','',sample_name, re.I)
    sample_name = re.sub(r'^.*\\','',sample_name, re.I)
    sample_name = re.sub(r'\.fasta$','',sample_name, re.I)
    sample_name = re.sub(r'\.fna$','',sample_name, re.I)
    sample_name = re.sub(r'\.faa$','',sample_name, re.I)
    sample_name = re.sub(r'\.fas$','',sample_name, re.I)
    sample_name = re.sub(r'\.fa$','',sample_name, re.I)

    return sample_name


def add_samples(folder, folders_samples):
     files = glob( folder + '/input/*.fasta')

     for file in files:
       sample_name =  extractSampleName(file)

       if not folder in  folders_samples:
          folders_samples[folder]  = {}

       folders_samples[folder][sample_name]  = {}


def check_file(file):
    return path.exists(file)


def isNotEmpty(file):
    size = 0
    if path.exists(file):
       try:
          size = path.getsize(file)
       except: 
          pass
    return size


def  check_type1(folders_samples, folder,  stages):

    for sample in folders_samples[folder].keys():
      if '1' in stages:
         status = check_file(folder + '/input/' + sample+'.fasta')
         if status:
            folders_samples[folder][sample]['1'] = 'Present'
         else:
            folders_samples[folder][sample]['1'] = 'Absent'

      if '2' in stages:
         filename =folder + '/output/' + sample+ '/orf_prediction/' + sample + '.qced.faa'
         status = check_file(filename)
         if status:
            folders_samples[folder][sample]['2'] = 'Present'
         else:
            folders_samples[folder][sample]['2'] = 'Absent'

      for db in get_db_names(stages, '3'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout"
           status = check_file(filename)
           if status:
              folders_samples[folder][sample]['3:' +db] = 'Present'
           else:
              folders_samples[folder][sample]['3:' + db] = 'Absent'

      for db in get_db_names(stages, '4'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout.parsed.txt"
           status = check_file(filename)
           if status:
              folders_samples[folder][sample]['4:' + db] = 'Present'
           else:
              folders_samples[folder][sample]['4:' + db] = 'Absent'


def  check_type2(folders_samples, folder,  stages):

    for sample in folders_samples[folder].keys():

      if '1' in stages:
         filename = folder + '/input/' + sample+'.fasta'
         size = isNotEmpty(filename)
         folders_samples[folder][sample]['1'] = size

      if '2' in stages:
         filename =folder + '/output/' + sample+ '/orf_prediction/' + sample + '.qced.faa'
         size = isNotEmpty(filename)
         folders_samples[folder][sample]['2'] = size

      for db in get_db_names(stages, '3'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout"
           size = isNotEmpty(filename)
           if size:
              folders_samples[folder][sample]['3:' + db ] = int(size)
           else:
              folders_samples[folder][sample]['3:' + db] =  int(size)

      for db in get_db_names(stages, '4'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout.parsed.txt"
           size = isNotEmpty(filename)
           if size:
              folders_samples[folder][sample]['4:' + db ] = int(size)
           else:
              folders_samples[folder][sample]['4:' + db ] =  int(size)




def check_type3(folders_samples, folder,  stages):

    i = 1
    for sample in sorted(folders_samples[folder].keys()):

      eprintf("  %3d\t%s\n",i, sample)
      i+=1

      if '1' in stages:
         filename = folder + '/input/' + sample+'.fasta'
         size = maxSizeFasta(filename)
         folders_samples[folder][sample]['1'] = int(size)

      if '2' in stages:
         filename =folder + '/output/' + sample+ '/orf_prediction/' + sample + '.qced.faa'
         size = numLines(filename)
         folders_samples[folder][sample]['2'] = size

      for db in get_db_names(stages, '3'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout"
           size = numLines(filename)
           if size:
              folders_samples[folder][sample]['3:' + db ] = int(size)
           else:
              folders_samples[folder][sample]['3:' + db] =  int(size)

      for db in get_db_names(stages, '4'):
           filename =folder + '/output/' + sample+ '/blast_results/' + sample + "." + db +".LASTout.parsed.txt"
           size = numLines(filename)
           if size:
              folders_samples[folder][sample]['4:' + db ] = int(size)
           else:
              folders_samples[folder][sample]['4:' + db ] =  int(size)


def get_db_names(stages, c):

     threePATT = re.compile(r'^' + c + ':')
     stage_suffix = {}
     for stage in stages.keys():
        result = threePATT.search(stage)
        if result:
           fields = [x.strip() for x in stage.split(':') ]
           stage_suffix[fields[1]] = True

     return sorted(stage_suffix.keys())

def check_type(folder_samples, folder,  stages, type):

     if type=='1':
         check_type1(folder_samples, folder,  stages)

     if type=='2':
         check_type2(folder_samples, folder,  stages)

     if type=='3':
         check_type3(folder_samples, folder,  stages)


def print_status(folders_samples, folder, sample,  stage, type) :
      printf('\t' + str(folders_samples[folder][sample][stage]))


# the main function
SIZE = 1000

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print usage
       sys.exit(0)

    stages = {}
    for stage in opts.stages:
        stages[stage] = True
    type = opts.type


    # adding the sampls
    folders_samples = {}
    for folder in opts.folders: 
       add_samples(folder, folders_samples)

 
    # stage processing 
    for folder in opts.folders: 
       eprintf("%s\n",folder)
       check_type(folders_samples, folder,  stages, opts.type) 


    printf("#FOLDER\tSAMPLE") 
    if '1' in stages:
        printf("\tINPUT_FILE") 

    if '2' in stages:
        printf("\tORF_FILE") 

    for db in get_db_names(stages, '3'):
        if '3:'+db in stages:
          printf("\t" + db + ".LASTout") 

    for db in get_db_names(stages, '4'):
        if '4:'+db in stages:
          printf("\t" + db + ".LASTout.parsed.txt") 
    printf("\n") 


    status1 = ''
    status2 = ''
    if type=='1':
       status1 = 'Y/N' 
       status2 = 'Y/N' 

    if type=='2':
       status1 = 'Size' 
       status2 = 'Size' 

    if type=='3':
       status1 = 'Avg Len' 
       status2 = 'Num Lines' 

    printf("#Name\tName") 
    if '1' in stages:
        printf("\t"+ status1) 
   
    if '2' in stages:
        printf("\t"+ status2) 
   
    for db in get_db_names(stages, '3'):
        if '3:'+db in stages:
           printf("\t"+ status2) 
   
    for db in get_db_names(stages, '4'):
        if '4:'+db in stages:
           printf("\t"+ status2) 
    printf("\n") 



    for folder in opts.folders: 
       for sample in sorted(folders_samples[folder].keys()):
          printf("%s\t%s",folder, sample) 

          if '1' in stages:
             print_status(folders_samples, folder, sample,  '1', opts.type) 

          if '2' in stages:
             print_status(folders_samples, folder, sample,  '2', opts.type) 

          for db in get_db_names(stages, '3'):
            if '3:'+db in stages:
               print_status(folders_samples, folder, sample,  '3:'+ db, opts.type) 

          for db in get_db_names(stages, '4'):
            if '4:'+db in stages:
               print_status(folders_samples, folder, sample,  '4:'+ db, opts.type) 

          printf("\n") 

  #  print folders_samples


     


# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

