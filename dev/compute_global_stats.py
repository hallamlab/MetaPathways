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


def read_stats_txt(sample, folder):

    global regs
    samplecyc = sample.lower() 
    filename = folder + "/" + sample + '/run_statistics/' + sample+ '.run.stats.txt'

    if not os.path.exists(filename):
      return None

    with open(filename, 'r') as fh:
        lines = [ x.strip() for x in fh.readlines()]

    stats = {}
    for reg in regs:
       for line in lines:
          res = reg[1].search(line)
          if res:
             fields = line.split('\t')
             stats[reg[0]]=fields[2]
       if reg[0] not in stats:
             stats[reg[0]]=0
    fh.close()
             
    try:
       #add_annotated_protein_counts(sample, folder)
       stats['NUM_PWYS'] = read_num_pwys(sample, folder)
       #stats['ORFS_ANNOTATED']= add_annotation_count(sample, folder)
       stats['ORFS_UNANNOTATED']= int( stats['NUM_AMINO (AFTER QC)']) - int( stats['ORFS_ANNOTATED']) 
       stats['PERCENT_ORFS_ANNOTATED']= "%.2f" %( float(stats['ORFS_ANNOTATED'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)
       stats['PERCENT_ANNOT_FROM_RefSeq'] = "%.2f" %( float(stats['ANNOT_FROM_RefSeq'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)
       stats['PERCENT_ANNOT_FROM_MetaCyc'] = "%.2f" %( float(stats['ANNOT_FROM_MetaCyc'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)
       stats['PERCENT_ANNOT_FROM_EggNog'] = "%.2f" %( float(stats['ANNOT_FROM_EggNog'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)
       stats['PERCENT_ANNOT_FROM_COG'] = "%.2f" %( float(stats['ANNOT_FROM_COG'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)
       stats['PERCENT_ANNOT_FROM_Kegg'] = "%.2f" %( float(stats['ANNOT_FROM_Kegg'])/float( stats['NUM_AMINO (AFTER QC)'])*100.00)

       #stats['PERCENT_RATIO1'] = "%.2f" %( float(stats['ANNOT_FROM_COG'])/ float(stats['ANNOT_FROM_RefSeq'])*100 )
       #stats['PERCENT_RATIO2'] = "%.2f" %( float(stats['ANNOT_FROM_MetaCyc'])/ float(stats['ANNOT_FROM_RefSeq'])*100 )
       #stats['PERCENT_RATIO3'] = "%.2f" %( float(stats['ANNOT_FROM_EggNog'])/ float(stats['ANNOT_FROM_RefSeq'])*100 )
       #stats['PERCENT_RATIO4'] = "%.2f" %( float(stats['ANNOT_FROM_Kegg'])/ float(stats['ANNOT_FROM_RefSeq'])*100 )
    except:
       print(traceback.print_exc(10))
       sys.exit(0)

    #if float(stats['PERCENT_ORFS_ANNOTATED']) < 30:
    #  print sample, float(stats['PERCENT_ORFS_ANNOTATED']) 
    add_annotated_protein_counts(sample, folder)

    return stats

def read_sample_data(sample, folder):
    printf("%s\t", sample)
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
       eprintf("%d\t%s\n", i, sample)
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

def print_data(samples, statsdata):
    for sample in samples:
        printf("\t%s", sample)
    printf("\n")
   
    for reg in  regs:
       printf("%25s",reg[0])
       for sample in samples:
            if reg[0] in statsdata[sample]:
               printf("\t%s",statsdata[sample][reg[0]])
            else:
               printf("\t0")
       printf("\n")

regs=[]

def print_stats(samples, statsdata):
    stats = {}    
    #stats['PERCENT_RATIOS1'] =[]
    stats['PERCENT_RATIOS2'] =[]
    stats['PERCENT_RATIOS3'] =[]
    stats['PERCENT_RATIOS4'] =[]

    
    #for sample in samples:
    #   stats['PERCENT_RATIOS1'].append(float(statsdata[sample]['PERCENT_RATIO1']))
    #   stats['PERCENT_RATIOS2'].append(float(statsdata[sample]['PERCENT_RATIO2']))
    #   stats['PERCENT_RATIOS3'].append(float(statsdata[sample]['PERCENT_RATIO3']))
    #   stats['PERCENT_RATIOS4'].append(float(statsdata[sample]['PERCENT_RATIO4']))

    metrics ={}
    #metrics['PERCENT_RATIOS1_AVG'] = scipy.mean(stats['PERCENT_RATIOS1'])
    metrics['PERCENT_RATIOS2_AVG'] = scipy.mean(stats['PERCENT_RATIOS2'])
    metrics['PERCENT_RATIOS3_AVG'] = scipy.mean(stats['PERCENT_RATIOS3'])
    metrics['PERCENT_RATIOS4_AVG'] = scipy.mean(stats['PERCENT_RATIOS4'])

    #metrics['PERCENT_RATIOS1_STD'] = scipy.std(stats['PERCENT_RATIOS1'])
    metrics['PERCENT_RATIOS2_STD'] = scipy.std(stats['PERCENT_RATIOS2'])
    metrics['PERCENT_RATIOS3_STD'] = scipy.std(stats['PERCENT_RATIOS3'])
    metrics['PERCENT_RATIOS4_STD'] = scipy.std(stats['PERCENT_RATIOS4'])


    a = 0.50
    i = 1
    headers=[ ('Sample', '30'), ('#COG', '8'), ('#MetaCyc', '8'), ('#EggNog', '8'), ('#Kegg', '8'), ('#RefSeq', '8'), ('#RefSeq(%)', '8')]
    
    eprintf("%5s", "#")
    for header in headers:
      eprintf("\t%" + header[1] + "s", header[0])
    eprintf("\n")

    for sample in samples:
       vals =[]
       flag = False
       count = 0
       if False and  abs( float(statsdata[sample]['PERCENT_RATIO1']) - metrics['PERCENT_RATIOS1_AVG'] ) >  a*metrics['PERCENT_RATIOS1_STD']:
          vals.append('COG')
          count += 1
          flag = True
       else:
          vals.append(' ')

       if abs( float(statsdata[sample]['PERCENT_RATIO2']) - metrics['PERCENT_RATIOS2_AVG'] ) >  a*metrics['PERCENT_RATIOS2_STD']:
          vals.append('MetaCyc')
          count += 1
          flag = True
       else:
          vals.append(' ')

       if abs( float(statsdata[sample]['PERCENT_RATIO3']) - metrics['PERCENT_RATIOS3_AVG'] ) >  a*metrics['PERCENT_RATIOS3_STD']:
          vals.append('EggNog')
          count += 1
          flag = True
       else:
          vals.append(' ')

       if abs( float(statsdata[sample]['PERCENT_RATIO4']) - metrics['PERCENT_RATIOS4_AVG'] ) >  a*metrics['PERCENT_RATIOS4_STD']:
          vals.append('Kegg')
          count += 1
          flag = True
       else:
          vals.append(' ')

       if count==4:
          vals.append('RefSeq')
       else:
          vals.append(' ')

       if flag == True:
          eprintf("%5d\t%30s", i,sample)
          i += 1
          for val in vals:
             eprintf("\t%8s", val)

          eprintf("\t%8s", statsdata[sample]['PERCENT_ANNOT_FROM_RefSeq'])
          eprintf("\n")


def setup_fields():
    global regs

    regs = [ 
              ['NUM_CONTIGS (BFORE QC)', re.compile(r'1000\tNumber of')],
              ['NUM_CONTIGS (AFTER QC)', re.compile(r'1005\tNumber of')],
              ['NUM_AMINO (BEFORE QC)', re.compile(r'2000\tNumber of translated')],
              ['NUM_AMINO (AFTER QC)', re.compile(r'2005\tNumber')],
              ['HITS_IN_COG', re.compile(r'hits in COG')],
              ['HITS_IN_MetaCyc', re.compile(r'hits in metacyc')],
              ['HITS_IN_Kegg', re.compile(r'hits in kegg')], 
              ['HITS_IN_EggNog', re.compile(r'hits in eggnog')], 
              ['HITS_IN_RefSeq', re.compile(r'hits in refseq')],
              ['ANNOT_FROM_COG', re.compile(r'Protein Annotations from COG')],
              ['ANNOT_FROM_MetaCyc', re.compile(r'Protein Annotations from metacyc')],
              ['ANNOT_FROM_Kegg', re.compile(r'Protein Annotations from kegg')],
              ['ANNOT_FROM_EggNog', re.compile(r'Protein Annotations from eggnog')],
              ['ANNOT_FROM_RefSeq', re.compile(r'Protein Annotations from refseq')],
         #     ['ORFS_ANNOTATED', re.compile(r'Note sure what Protein Annotations from refseq-2016-10-06-rel-78')],
              ['PERCENT_ANNOT_FROM_COG', re.compile(r'Is calculated COG')],
              ['PERCENT_ANNOT_FROM_MetaCyc', re.compile(r'Is calculated MetaCyc')],
              ['PERCENT_ANNOT_FROM_Kegg', re.compile(r'Is calculated Kegg')],
              ['PERCENT_ANNOT_FROM_EggNog', re.compile(r'Is calculated EggNog')],
              ['PERCENT_ANNOT_FROM_RefSeq', re.compile(r'Is calculated refseq')],
              ['NUM_PWYS', re.compile(r'Not supposed to match anything')],
              ['ORFS_ANNOTATED', re.compile(r'Number of ORFS Annotated in Total') ],
              ['ORFS_UNANNOTATED', re.compile(r'No supposed to match') ],
              ['PERCENT_ORFS_ANNOTATED', re.compile(r'No supposed to match') ],
              ['PERCENT_RATIO1', re.compile(r'No supposed to match') ],
              ['PERCENT_RATIO2', re.compile(r'No supposed to match') ],
              ['PERCENT_RATIO3', re.compile(r'No supposed to match') ],
              ['PERCENT_RATIO3', re.compile(r'No supposed to match') ],
           ]  



def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    if not valid_arguments(opts, args):
       print(usage)
       sys.exit(0)

    setup_fields()
    samples = read_list(opts.listfile) 

    statsdata = read_sample(samples, opts.folder)

    print_data(samples, statsdata)
    
    print_stats(samples, statsdata)


# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

