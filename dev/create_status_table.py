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

except:
     print(""" Could not load some user defined  module functions""")
     print(""" Make sure your typed 'source MetaPathwaysrc'""")
     print(""" """)
     sys.exit(3)




usage= sys.argv[0] + """ -s <sample_name_list_file>  -m <metapathways_output_folder> """

parser = None
def createParser():
    global parser
    epilog = """This script generates a table of the stats for a list of samples in MP output folder."""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-s", "--sample-list-file", dest="sample_list_file",
                      help='a file with sample list [REQUIRED]')

    parser.add_option("-m", "--mp-output_folder", dest="mp_output_folder",
                      help='output file with mp output for several samples[REQUIRED]')

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def main(argv, errorlogger = None, runstatslogger = None): 
    global parser
    (opts, args) = parser.parse_args(argv)

    
    
    with open(opts.sample_list_file, 'r') as fp:
       samples = [x.strip() for x in fp.readlines()]
       
    stats = {}
    for sample in samples:
       stats[sample] = read_runs_stats_file( opts.mp_output_folder + "/" + sample+"/run_statistics/" + sample + ".run.stats.txt" )

    keys={}
    for sample in samples:
      for key in stats[sample]:
        keys[key]=(int(key), stats[sample][key][0])

    key_list = sorted([ (x, y[0], y[1]) for x, y in keys.items()], key=lambda x: x[1])
    #print(key_list)

 
    print('\t\t' +'\t'.join(samples))

    for key, _,  des  in key_list:
       row = [key, des]
       for sample in samples:
          if not key in stats[sample]:
             row.append( '0' )
          else:
             row.append( stats[sample][key][1] )

       print('\t'.join(row))





def read_runs_stats_file(file_name):
    with open(file_name, 'r') as fp:
       lines = [x.strip() for x in fp.readlines()]

    stats = {}
    for line in lines:
       fields = [x.strip() for x in line.split('\t') ]
       stats[fields[0]] = [ fields[1], fields[2]]
    return stats

# the main function of metapaths
if __name__ == "__main__":
    createParser()
    main(sys.argv[1:])

