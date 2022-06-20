#!/usr/bin/python

"""This script run the pathologic """

try:
   import optparse, sys, re, csv, traceback
   from optparse import OptionGroup
   import pickle
   import math
   from libs.python_modules.taxonomy.LCAComputation import *
   import operator

   from os import path, _exit, remove, rename
   import logging.handlers
   from glob import glob
   from libs.python_modules.utils.sysutil import pathDelim
   from libs.python_modules.utils.metapathways_utils  import fprintf, printf, eprintf,  exit_process
   from libs.python_modules.utils.sysutil import getstatusoutput

   from libs.python_modules.utils.pathwaytoolsutils import *

except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed 'source MetaPathwaysrc'"""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)

PATHDELIM=pathDelim()

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)

def files_exist( files , errorlogger = None):
    status = True    
    for file in files:
       if not path.exists(file):
          if errorlogger:
             errorlogger.write( 'ERROR\tCould not find ptools input  file : ' +  file )
          status = False
    return not status



usage = sys.argv[0] + """ -s sample -p pgdb_dir --ptoolsExec pathwaytools_executable """
parser = None
def createParser():
    global parser

    epilog = """The flat file extraction script"""

    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)
    standard_options_group = OptionGroup(parser, "Standard Ptools group" )
    # Input options
    standard_options_group.add_option('-s', '--sample', dest='sample_name', default=None,
                           help='sample name')

    standard_options_group.add_option('-p', '--pgdb', dest='pgdbdir', default=None,
                           help='folder of the PGDB')

    standard_options_group.add_option('--ptoolsExec', dest='ptoolsExec', default=None,
                           help='PathoLogic Executable')

    standard_options_group.add_option("-o", "--output-pwy-table", dest="table_out",
        help='the output table for the pathways [REQUIRED]')

import os, signal
TIME = 10

def __StopPathwayTools():
    processPATT = re.compile(r'pathway-tools-runtime')
    for line in os.popen("ps xa"):
        fields = line.split()
        pid = fields[0]
        process = fields[4]
        result = processPATT.search(process)
        if result :
            os.kill(int(pid), signal.SIGHUP)


def StopPathwayTools():
  try:
     __StopPathwayTools()
     time.sleep(TIME)
     __StopPathwayTools()
     time.sleep(TIME)

     if path.exists("/tmp/ptools-socket"): 
        remove("/tmp/ptools-socket")
  except:
    pass


def main(argv, errorlogger = None, runcommand = None, runstatslogger = None):
    global parser

    options, args = parser.parse_args(argv)

    # is there a pathwaytools executable installed
    if False and not path.exists(options.ptoolsExec):
       eprintf("ERROR\tPathwayTools executable %s not found!\n", options.ptoolsExec)
       if errorlogger:
          errorlogger.printf("ERROR\tPathwayTools executable %s not found!\n",  options.ptoolsExec)
       exit_process("ERROR\tPathwayTools executable %s not found!\n" %(options.ptoolsExec))


    # command to build the ePGDB
    command = "%s "  %(options.ptoolsExec)
    command += " -api"

    pythonCyc = startPathwayTools(options.sample_name.lower(), options.ptoolsExec, True)
    #resultLines = pythonCyc.getReactionListLines()
    resultLines = pythonCyc.getFlatFiles()
    StopPathwayTools()
    try:
      if False:
         pythonCyc = startPathwayTools(options.sample_name.lower(), options.ptoolsExec, True)
         pythonCyc.setDebug() # disable pathway debug statements
         printf("INFO\tExtracting the reaction list from ePGDB " + options.sample_name + "\n")
         resultLines = pythonCyc.getReactionListLines()
         #pythonCyc.stopPathwayTools()
         reaction_list_file = open(options.reactions_list + ".tmp", 'w')
         for line in resultLines:
          fprintf(reaction_list_file,"%s\n",line.strip())
         reaction_list_file.close()
         StopPathwayTools()

    except:
           print traceback.print_exc(10)
           eprintf("ERROR\tFailed to run extract pathways for %s : \n" %(options.sample_name))
           eprintf("INFO\tKill any other PathwayTools instance running on the machine and try again")
           if errorlogger:
               errorlogger.write("ERROR\tFailed to run extract pathways for %s : " %(options.sample_name))
               errorlogger.write("INFO\tKill any other PathwayTools instance running on the machine and try again\n")
           StopPathwayTools()


def startPathwayTools(organism, ptoolsExec, debug):
    StopPathwayTools()
    pythonCyc = PythonCyc()
    pythonCyc.setDebug(debug = debug)
    pythonCyc.setOrganism(organism)
    pythonCyc.setPToolsExec(ptoolsExec)
    pythonCyc.startPathwayTools()
    return pythonCyc


def runPathologicCommand(runcommand = None):
    if runcommand == None:
      return False
    result = getstatusoutput(runcommand)
    return result[0]


# this is the portion of the code that fixes the name

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes



def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
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


def cleanup(string):
    """
    Cleans up pathway long-names for presentation.
    :param string:
    :return:
    """
    string = re.sub("|", "", string) # vertical bar
    string = re.sub("&", "", string) # ampersand
    string = re.sub(";", "", string) # semicolon
    string = re.sub("<[^<]+?>", '', string) # HTML tags
    string = re.sub("\'", "", string) # remove quotes

    return string

def get_preferred_taxa_name(taxa_id, megan_map, id_to_name):
    """
    Helper function to format NCBI IDs into preferred names. First checks for MEGAN name,
    if not found moves to current taxonomy in loaded NCBI taxonomy tree, failing that
    gives the taxonomy of 'Unknown', but still provides the id, e.g., 'Unknown (12345)'.
    :param taxa_id: numeric taxa id to translate
    :param megan_map: preferred megan mapping hash
    :param id_to_name: local ncbi tree hash
    :return: "perferred name (id)"
    """
    taxa_id = str(taxa_id)
    if taxa_id in megan_map:
        taxa = megan_map[ taxa_id ] + " (" + taxa_id + ")"
    elif taxa_id in id_to_name:
        taxa = id_to_name[ taxa_id ] + " (" + taxa_id + ")"
    else:
        taxa = "Unknown" + " (" + taxa_id + ")"

    return taxa

def MetaPathways_run_pathologic(argv, extra_command = None, errorlogger = None, runstatslogger =None): 
    if errorlogger != None:
       errorlogger.write("#STEP\tBUILD_PGDB\n")
    createParser()
    main(argv, errorlogger = errorlogger, runcommand= extra_command, runstatslogger = runstatslogger)
    return (0,'')

if __name__ == '__main__':
    createParser()
    main(sys.argv[1:])

