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
     import sys, os, re, math, scipy
     import traceback
     from os import makedirs, sys, remove, rename
     from sys import path
     from optparse import OptionParser

except:
     print(""" Could not load some user defined  module functions""")
     print(""" Make sure your typed 'source MetaPathwaysrc'""")
     print(""" """)
     print(traceback.print_exc(10))
     sys.exit(3)



usage= sys.argv[0] + """ --list <samples_in_a_file> --folder <folders with samples>"""

def fprintf(file, fmt, *args):
   file.write(fmt % args)
   
def printf(fmt, *args):
   sys.stdout.write(fmt % args)
   sys.stdout.flush()
 
def eprintf(fmt, *args):
   sys.stderr.write(fmt % args)
   sys.stderr.flush()


parser = None
def createParser():
    global parser
    epilog = """
This script computes the sequence stats for the fasta files
"""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-l", "--list", dest="listfile", help='the list of samples in a file; one sample per row [REQUIRED]')

    parser.add_option("-f", "--folder", dest="folder", help='folder where the samples are located [REQUIRED]')


def valid_arguments(opts, args):
    state = True
    if opts.listfile == None :
        print('ERROR: Missing sample list file')
        state = False

    if opts.folder == None :
        print('ERROR: Missing processed samples folder')
        state = False
    return state


def read_list(listfile):
     samples = []
     with open(listfile, 'r') as fh:
         samples = [ x.strip() for x in fh.readlines()]
     return samples
          

def get_orfs(string):
    string = re.sub(r'[\[\]]','', string)  
    orfs = string.split(',')
    return orfs

def read_num_pwys(sample, folder):
    samplecyc = sample.lower() 
    filename = folder + "/" + sample +"/" + 'results/pgdb/' + samplecyc + '.pwy.txt'

    if not os.path.exists(filename):
      return None

    with open(filename, 'r') as fh:
        lines = [ x.strip() for x in fh.readlines()]

    return str(len(lines))

def add_annotation_count(sample, folder):
    filename = folder + "/" + sample +"/" + 'results/annotation_table/' + sample + '.ORF_annotation_table.txt'

    if not os.path.exists(filename):
      return 0

    with open(filename, 'r') as fh:
        lines = [ x.strip() for x in fh.readlines()]
    fh.close()

    return str(len(lines))


def add_annotated_protein_counts(sample, folder):

    filename = folder + "/" + sample +"/" + 'results/annotation_table/' + sample + '.functional_and_taxonomic_table.txt'
    if not os.path.exists(filename):
      return 0

    commentPATT = re.compile(r'^#')
    count = -1 # get rid of th header
    with open(filename, 'r') as fh:
        lines = [ x.strip() for x in fh.readlines()]
        for line in lines:
           if not commentPATT.search(line):
              count += 1
    fh.close()

    filename = folder + "/" + sample + '/run_statistics/' + sample+ '.run.stats.txt'
    with open(filename, 'a') as fh:
       fprintf(fh, "%s\t%s\t%s\n", "8000", "Number of ORFS Annotated in Total", count)


code_to_description = {}
def read_stats_txt(sample, folder):
    global regs
    global code_to_description 
    samplecyc = sample.lower() 
    filename = folder + "/" + sample + '/run_statistics/' + sample+ '.run.stats.txt'

    if not os.path.exists(filename):
      return None

    with open(filename, 'r') as fh:
        lines = [ x.strip() for x in fh.readlines()]

    code_to_count = {}

    stats = {}
    for line in lines:
      fields = line.split('\t')
      code_to_count[int(fields[0])] = fields[2]
      code_to_description[int(fields[0])] = fields[1]
    fh.close()
             
    #if float(stats['PERCENT_ORFS_ANNOTATED']) < 30:
    #  print sample, float(stats['PERCENT_ORFS_ANNOTATED']) 
    #add_annotated_protein_counts(sample, folder)

    return code_to_count

def read_sample_data(sample, folder):
    #printf("%s\t", sample)
    stats = read_stats_txt(sample, folder)

    return stats

def orf_in_sample(pwyorfs):
    count= 0
    for pwy in pwyorfs.keys():
        count +=  pwyorfs[pwy]

    return count
    

def read_sample(samples, folder):
    statsdata = {}

    i= 1
    for sample in samples:
       print("{}\t{}".format(i, sample))
       i += 1
       data =  read_sample_data(sample, folder) 
       if data!=None:
          statsdata[sample] = data
          #count = orf_in_sample(pwydata[sample])
          #print sample, len(pwydata[sample].keys()), count
       else:
          print(sample, 'MISSING')
          sys.exit(0)
    return statsdata

regs=[]

def print_stats(samples, statsdata):
   pass


def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print(usage)
       sys.exit(0)

    samples = read_list(opts.listfile) 
    statsdata = read_sample(samples, opts.folder)

    code_to_des = [ (int(code), desc) for code, desc in  code_to_description.items() ]
    code_to_des.sort(key = lambda x: x[0], reverse=False)
    
   
    for code, desc in code_to_des:
      vals = [str(code),  desc]
      for sample in samples:
        try:
           vals.append(statsdata[sample][code])
        except:
           pass

      print('\t'.join(vals))



# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

