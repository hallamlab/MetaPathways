#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re, shutil
     import sys, os
     from optparse import OptionParser, OptionGroup
     from glob import glob
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)



PATHDELIM = "/"
script_name = sys.argv[0]
usage= script_name + """ --in-folder  <input file> --out-folder <output files>  -c <combined_sample> -s <sample> [ -s <sample> ... -s <samplen> ]"""
parser = OptionParser(usage)

parser.add_option( "--in-folder", dest="input_folder", default =None, 
                  help='the input folder where the MetaPathways output reside')

parser.add_option( "--out-folder", dest="output_folder", default =None, 
                  help='the output folder where the output is going')

parser.add_option( "-s", "--single-sample", dest="samples", default =[], action='append',
                  help='individual sample snames [defaulty is empty]')

parser.add_option( "-c", "--combined-sample", dest="combined_sample", default = sys.stdout, 
                  help='the name of the combined samples')


def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)



PATHDELIM="/"

def read_in_file(filel, readdir):

     _fastqfiles = glob(readdir + PATHDELIM +  'E*.fastq')

     fastqfiles =[]
     for _f in _fastqfiles:
         f = re.sub(r'^.*[//]','', _f)
         fastqfiles.append(f)

     fastq_groups = get_fastq_pairs(fastqfiles)
     #print fastq_groups

     with open(filel,'r') as infile:
        lines = infile.readlines()

     for line in lines:
         fields = line.strip().split('\t')
         if len(fields)!=2:
            continue
         samplefastqs = fields[1].split('|')
         print 
         print fields[0], samplefastqs
         paired_files = get_paired_files(fastq_groups, samplefastqs)
         b = 1
         for f1 in paired_files:
             for  f2 in paired_files[f1]: 
                 src = f1 + "_"+ f2 + ".fastq"
                 dest = fields[0] + "_"+ f2 + ".b" +str(b) + ".fastq"
                 #printf("%s/ %s %s",readdir,  src, dest)
                 os.rename(readdir +"/" + src, readdir +"/"+ dest)
             b +=1
             print 
         
         

     return fields

def get_fastq_pairs(fastqfiles):
    fastq_structure ={}
    for fastq in fastqfiles: 
        f = re.sub('_[1-2].fastq','',fastq)
        if not f in fastq_structure:
           fastq_structure[f]=[]

        if fastq == f +"_1.fastq":
          fastq_structure[f].append('1')
        else:
          fastq_structure[f].append('2')

    return fastq_structure


def get_paired_files(fastq_groups, samplefastqs):
    pairedfiles={}
    for sfastq in samplefastqs:
       if sfastq in fastq_groups  :
          pairedfiles[sfastq] = fastq_groups[sfastq]
      


    return pairedfiles

def write_new_file(lines, output_file):
    
    print "Fixing file " + output_file 
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print "ERROR :Cannot open output file "  + output_file
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def check_arguments(opts):
    if opts.input_folder==None:
       return False
    if opts.output_folder==None:
       return False

    if len(opts.samples)==0:
       return False
    if opts.combined_sample==None:
       return False
    return True




def read_0pf_file(functions, sample_data, input_folder, sample):
     fpfile = input_folder + PATHDELIM + sample + PATHDELIM + "ptools" + PATHDELIM + "0.pf" 
     with open(fpfile,'r') as infile:
        lines = infile.readlines()


     sample_data[sample]['0.pf']=[]
     orf=[]
     function =''
     orf_name = ''
     ec = ''
     for _line in lines:
         line  = _line.strip()
         fields = _line.strip().split('\t')
      
         if fields[0]=="ID" or fields[0]=="NAME":
            orf_name = fields[1].replace('O_',sample+"_") 
            line = fields[0] + "\t" + orf_name

         if fields[0]=="FUNCTION":
            function=fields[1]

         if fields[0]=="EC":
            ec=fields[1]

         orf.append(line)
         if fields[0]=="//":
            #decide on orf orf.append(orf)
            if function in functions:
                sample_data['reduced.txt'][orf_name] = functions[function]
                pass
            else:
               functions[function]=orf_name
               sample_data[sample]['0.pf'].append('\n'.join(orf))

            function = ''
            ec = ''
            orf=[]
 
         #sample_data[sample]['0.pf'].append(line)

def read_reduced_file(sample_data, input_folder, sample):
     fpfile = input_folder + PATHDELIM + sample + PATHDELIM + "ptools" + PATHDELIM + "reduced.txt" 

     with open(fpfile,'r') as infile:
        lines = infile.readlines()

     for _line in lines:
         line  = _line.strip()
         fields = _line.strip().split('\t')
         field0 = fields[0].replace('O_',sample+"_") 
         field1 = fields[1].replace('O_',sample+"_") 
         newline = field0 +'\t' +field1 
         sample_data['reduced.txt'][field0]=[field1]


def load_sample_data(functions, sample_data, input_folder, sample):
     sample_data[sample] = {}

     read_reduced_file(sample_data, input_folder, sample)

     read_0pf_file(functions, sample_data, input_folder, sample)


def  write_0pf_file(sample_data, fout):
     for line in sample_data['0.pf']:
        fprintf(fout, "%s\n",line)

def  write_reduced_file(sample_data, fout):
     for key, value in sample_data['reduced.txt'].iteritems():
        fprintf(fout, "%s\t%s\n",key, value)


def  write_genetic_elements_file(filename):
    line = "ID\t0\nNAME\t0\nTYPE\t:CONTIG\nANNOT-FILE\t0.pf\n//"

    with open(filename, 'w') as f:
        fprintf(f, "%s\n",line)


def  write_organism_params(filename, sample_name):
    contents=[ 
               ["ID", sample_name], ["STORAGE",	"FILE"],\
               ["NAME", sample_name], ["ABBREV-NAME", sample_name],\
               ["STRAIN", "1"], ["RANK", "|species|"],\
               ["NCBI-TAXON-ID", "12908"] 
             ]

    with open(filename, 'w') as f:
       for line in contents:
           fprintf(f, "%s\t%s\n",line[0], line[1])
       

def write_combined_sample_data(sample_data, output_folder, combined_sample):

    samplefolder= output_folder + PATHDELIM + combined_sample 
    if os.path.exists(samplefolder):
       shutil.rmtree(samplefolder)
    os.mkdir(samplefolder)
    os.mkdir(samplefolder + PATHDELIM + "ptools")

    fpfile = samplefolder + PATHDELIM + "ptools"+ PATHDELIM + "0.pf" 
    try:
       fout = open(fpfile,'w')
    except:
       print "ERROR: cannot write output 0.pf file sample"

    for sample in sample_data:
       if sample!='reduced.txt':
         write_0pf_file(sample_data[sample], fout)
    fout.close()


    fpfile = samplefolder + PATHDELIM + "ptools"+ PATHDELIM + "reduced.txt" 
    try:
       fout = open(fpfile,'w')
    except:
       print "ERROR: cannot write output 0.pf file sample"

    #for sample in sample_data:
    write_reduced_file(sample_data, fout)

    write_genetic_elements_file(samplefolder + PATHDELIM + "ptools"+ PATHDELIM + "genetic-elements.dat")

    write_organism_params(samplefolder + PATHDELIM + "ptools"+ PATHDELIM + "organism-params.dat", combined_sample)

    fout.close()
     

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts):
        print "Usage: ",  usage
        sys.exit(0)
     
    sample_data = {}
    functions = {} 
    sample_data['reduced.txt']={}

    for sample in opts.samples:
       load_sample_data(functions, sample_data, opts.input_folder, sample)
    
    write_combined_sample_data(sample_data, opts.output_folder, opts.combined_sample)
    
    #filenames=read_in_file(opts.input, opts.folder)


# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

