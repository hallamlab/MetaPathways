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
     import re
     import sys
     from optparse import OptionParser, OptionGroup
     from glob import glob
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

script_name = sys.argv[0]
usage= script_name + """ --pwy-file <pwyfile> --reduced  <reduced file> -o/--output  <output>"""
parser = OptionParser(usage)

parser.add_option( "--pwy-file", dest="pwy_file",  help='the pwy file')

parser.add_option( "--reduced", dest="reduced",  help='the reduced file')

parser.add_option( "-o", "--output", dest="output",  help='the output file')


def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def check_arguments(opts, args):
    if opts.pwy_file == None:
         return False
    if opts.reduced== None:
         return False
    return True


def read_equivalent_file(filename):
     reduced = {}
     with open(filename,'r') as infile:
      for line in infile:
         line = line.strip()
         fields = [ x.strip() for x in line.split('\t') ]
         if len(fields) ==2:
             if not  fields[1] in reduced:
                reduced[fields[1]] = []
             reduced[fields[1]].append(fields[0])
     return reduced
          

def fix_pwy_file(pwy_txt_file,  output_file, equiv):
    orfPATT = re.compile(r'\[(.*)\]')
    with open(pwy_txt_file,'r') as infile, open(output_file, 'w') as outfile:
       for _line in infile:
          line = _line.strip()
          res = orfPATT.search(line)
          if res:
             orf_str = ''
             orfs = [ x.strip() for x in res.group(1).split(',') ]
             for orf in orfs:
               if orf_str=='':
                  orf_str = orf
               else:
                  orf_str += ','+ orf
               if orf in equiv:
                  new_str = ','.join(equiv[orf])
                  orf_str += ','+ new_str
             orf_str = "[" + orf_str + "]"  
             new_line = re.sub(orfPATT, orf_str, line)
             fprintf(outfile, "%s\n", new_line)
          else:
             fprintf(outfile, "%s\n", line)


# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)


# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    reduced_equiv = read_equivalent_file(opts.reduced)

    fix_pwy_file(opts.pwy_file, opts.output, reduced_equiv)

    
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

