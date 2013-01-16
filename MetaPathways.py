#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2012, The metapaths Project"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

import sys
from os import makedirs, sys, listdir, environ, path
import re 
import inspect
from commands import getstatusoutput
cmd_folder = path.abspath(path.split(inspect.getfile( inspect.currentframe() ))[0])



#print cmd_folder
#if not sys.platform.startswith('win'):
#    res =getstatusoutput('source  '+ cmd_folder +'/'+'.metapathsrc')
#    if( int(res[0])==0 ): 
#       print 'Ran ' + cmd_folder +'/'+'.metapathsrc ' + ' file successfully!'
#    else:
#       print 'Error : ' + res[1] 
#       print 'while running  ' + cmd_folder +'/'+'.metapathsrc ' + ' file!'

#sys.path.insert(0,cmd_folder + "/libs/python_modules/")
#sys.path.insert(1, cmd_folder + "/libs/")
#print sys.path


from optparse import OptionParser
import shutil 

from libs.python_modules import metapaths_utils
from libs.python_modules.metapaths_utils  import parse_command_line_parameters
from libs.python_modules.parse  import parse_metapaths_parameters
from libs.python_modules.metapaths_pipeline import print_commands, call_commands_serially, print_to_stdout, no_status_updates

from libs.python_modules.metapaths import run_metapathways
from libs.python_modules.annotate import *

#config = load_config()
metapaths_config = """template_config.txt""";

script_info={}
script_info['brief_description'] = """A workflow script for making PGDBs from metagenomic sequences"""
script_info['script_description'] = """ takes a sequence file and performs all processing steps through building the OTU table.
             REQUIRED: You must have a fas and an yaml file  and  a custom parameters file:"""
script_info['script_usage'] = []

usage= """./run_pgdb_pipeline.py  -i input_file -o outdir  -p parameters.txt """
parser = OptionParser(usage)
parser.add_option("-i", "--input_file", dest="input_fp",
                  help='the input fasta file/input dir [REQUIRED]')
parser.add_option("-o", "--output_dir", dest="output_dir",
                  help='the input fasta file/input dir [REQUIRED]')
parser.add_option('-p','--parameter_fp', dest="parameter_fp",
                   help='path to the parameter file [REQUIRED]')
parser.add_option("-c", "--config_filer", dest="config_file",
                  help='pipeline_configuratin file [OPTIONAL,  default : \"MetaPathways/template_config.txt\"]')
parser.add_option('-r','--run-type', dest="run_type", default='safe',
                   choices=['safe', 'overlay', 'overwrite','dry-run'], 
                   help= '\n(a) \'overwrite\' -- wipes out the previous runs with the same name\n'+
                         '\n(b)\'overlay\' -- recomputes the steps that are not present \n' +
                         '\n(c)\'dry-run\' -- shows the steps that are going to be computed or not\n' +
                         '\n(d)\'safe\' -- safe mode does not run on an existing run folder\n')
#ith out of order completion \ time-stamps in the \'workflow_log.txt\' 
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="print lots of information on the stdout [default]")
parser.add_option("-P", "--print-only",
                  action="store_true", dest="print_only", default=False,
                  help="print only  the commands [default False]")


# checks if the supplied arguments are adequate
def check_arguments(opts, args):
    if (opts.input_fp == None and opts.output_dir ==None )  or\
     opts.output_dir == None or opts.parameter_fp == None :
       return True
    else:
       return False


def create_input_fasta_yaml_pairs(input_dir, output_dir):
    fileslist =  listdir(input_dir)
    input_files = {}
    for file in fileslist:
       shortname = re.sub('.(fasta|fas|fna|faa)','',file) 
       if re.search('.(fasta|fas|fna|faa)',file):
         input_files[file] = shortname

    paired_input = {} 
    for key, value in input_files.iteritems():
            paired_input[input_dir+'/'+key]=output_dir + '/'+ value

    return paired_input

def openGrades():
    pass

def openRank():
    pass

# main function




def main(argv):
    (opts, args) = parser.parse_args()
    if check_arguments(opts, args):
       print usage
       sys.exit(0)


    # initialize the input directory or file

    input_fp = opts.input_fp 

    output_dir = opts.output_dir
    verbose = opts.verbose
    print_only = opts.print_only

    run_type = opts.run_type.strip()

    if run_type == 'overwrite':
       force_remove_dir= True
    else:
       force_remove_dir=None

    if opts.config_file:
       config_file= opts.config_file
    else:
       config_file = cmd_folder + "/" + metapaths_config
    

    # try to load the parameter file    
    try:
        parameter_f = open(opts.parameter_fp)
    except IOError:
        raise IOError,\
         "Can't open parameters file (%s). Does it exist? Do you have read access?"\
         % opts.parameter_fp


    if force_remove_dir:
        try: 
           if path.exists(output_dir):
              shutil.rmtree(output_dir)
        except OSError:
           print "ERROR: Cannot remove directory: " + output_dir
           sys.exit(1)

    try:
       if (run_type in ['overlay', 'safe'] or force_remove_dir) and not path.exists(output_dir):
             makedirs(output_dir)
       elif run_type in ['safe'] and path.exists(output_dir):
        print ""
        print "ERROR: Cannot create output directory \"" + output_dir + "\"\n"+\
              "       Perhaps output directory already exists.\n" +\
              "       Please choose a different directory, or \n" +\
              "       run with the option \"-r  overwrite\" to force overwrite it."
        sys.exit(1)
    except OSError:
        print ""
        print "ERROR: Cannot create output directory \"" + output_dir + "\"\n"+\
              "       Perhaps directory \"" + output_dir  + "\" already exists.\n" +\
              "       Please choose a different directory, or \n" +\
              "       run with the option \"-r  overwrite\" to force overwrite it."
        sys.exit(1)
        
    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially
    
    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates
    

    command_line_params={}
    command_line_params['verbose']= opts.verbose


    # load the sample inputs  it expects either a fasta file or  a directory 
    # containing fasta and yaml file pairs

    input_output_list = {}
    if path.isfile(input_fp):   # check if it is a file
       input_output_list[input_fp] = output_dir
    else:
       if path.exists(input_fp):   # check if dir exists
          input_output_list = create_input_fasta_yaml_pairs(input_fp, output_dir)
       else:   # must be an error
          print "No valid input sample file or directory containing samples exists .!"
          print "As provided as arguments in the -in option.!"
          sys.exit(1)
   

    config_params=parse_metapaths_parameters(parameter_f)

    # add check the config parameters 

    for input_fp, output_dir in input_output_list.iteritems():
        
        if run_type=='overwrite' and  path.exists(output_dir):
           shutil.rmtree(output_dir)
           makedirs(output_dir)
        if not  path.exists(output_dir):
           makedirs(output_dir)
        
        run_metapathways(
           input_fp, 
           output_dir,
           command_handler=command_handler,
           command_line_params=command_line_params,
           config_params=config_params,
           metapaths_config=metapaths_config,
           status_update_callback=status_update_callback,
           config_file=config_file,
           run_type = run_type
        )

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])    

