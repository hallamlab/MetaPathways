#!/usr/bin/env python
# File created on 30 Dec 2009.



from __future__ import division
from subprocess import Popen, PIPE, STDOUT
from os import makedirs, listdir

from glob import glob
from os.path import split, splitext, join, dirname, abspath
from datetime import datetime

from metapaths_utils import printf, eprintf
from sysutil import getstatusoutput
import sys


from BlastGrid import *

from optparse import OptionParser

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2010, The metapaths Project"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""
This file contains the metapaths workflow functions which string together 
independent scripts. 
"""

## Start utilities used by the pipeline functions
def generate_log_fp(output_dir,
                    basefile_name='meta_paths_run_log',
                    suffix='txt',
                    timestamp_pattern=''):
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s.%s' % (basefile_name,suffix)
    return join(output_dir,filename)

def generate_steps_log_fp(output_dir,
                    basefile_name='meta_paths_steps_log',
                    suffix='txt'):
    filename = '%s.%s' % (basefile_name,suffix)

    return join(output_dir,filename)

class WorkflowError(Exception):
    pass

class WorkflowLogger(object):
    
    def __init__(self,log_fp=None,params=None,metapaths_config=None,open_mode='w'):
        if log_fp:
            self._f = open(log_fp,open_mode)
        else:
            self._f = None
        self._filename = log_fp
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('#Logging started at %s\n\n' % start_time)
        self.writemetapathsConfig(metapaths_config)
        self.writeParams(params)
    def get_log_filename(self): 
        return self._filename

    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass
    
    def writemetapathsConfig(self,metapaths_config):
        if metapaths_config == None:
            self.write('#No metapaths config provided.\n')
        else:
            self.write('#metapaths_config values:\n')
            for k,v in metapaths_config.items():
                if v:
                    self.write('%s\t%s\n' % (k,v))
            self.write('\n')
            
    def writeParams(self,params):
        if params == None:
            self.write('#No params provided.\n')
        else:
            self.write('#parameter file values:\n')
            for k,v in params.items():
                for inner_k,inner_v in v.items():
                    val = inner_v or 'True'
                    self.write('%s:%s\t%s\n' % (k,inner_k,val))
            self.write('\n')
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

def print_commands(commands,
                   status_update_callback,
                   logger):
    """Print list of commands to run """
    #logger.write("Printing commands only.\n\n")
    #for c in commands:
    #    for e in c:
    #        status_update_callback('#%s' % e[0])
    #        print '%s' % e
    #        logger.write('# %s command\n%s\n\n' % e)
            
def call_commands_serially(commands, status_update_callback, logger, stepslogger, params ):
    """Run list of commands, one after another """
    #logger.write("Executing commands.\n\n")
    from sys import exit
    make = False
    for c in commands:
        if c[3]=='stop':
           print "Stopping!"
           return 
           break

        if params['verbose']:
             eprintf("\n\n\nIssuing Command : %s\n", c[1])

        eprintf("%s" %(c[0]))
        if c[3] in ['yes', 'redo' ] and c[4]:
           result = getstatusoutput(c[1])
           if result[0] == 0 :
             if c[3] in ['redo']:
                eprintf('..... Redo Success!\n')
             else:
                eprintf('..... Success!\n')
           else:
             eprintf('..... Error!\n')
             #print c
        elif c[3]=='grid' and c[4]:
           blastgrid(c[1])
           continue  
        elif c[3]=='skip':
           eprintf('..... Skipping!\n')
           continue
        else:
           eprintf('..... Already Computed!\n')
           continue

        #status_update_callback('%s\n%s' %(result[0], result[0]))

        logger.write('COMMAND : \n%s \n' %( c[1]))
        if result[0] == 0 :
            logger.write('Success!:\n%s\n' %( result[1]))
        else:
            print result[1]

            printf('Error! : %s\n' %( result[1]))
            logger.write('Error!:\n%s\n' %( result[1]))
            break 

        timestamp_pattern='%Y-%m-%d %H:%M:%S'
        timestamp = datetime.now().strftime(timestamp_pattern)
        stepslogger.write('%s\t%s\n' %( c[2], timestamp))


            #proc = Popen(e[1],shell=True,universal_newlines=True, stdout=PIPE,stderr=PIPE)
            # communicate pulls all stdout/stderr from the PIPEs to 
            # avoid blocking -- don't remove this line!
            #stdout, stderr = proc.communicate()
            #print stderr

            #if proc.returncode==0:
            #     print """ ....done Successfully!\n""" 
            #else:
            #     print """ ....done Uncuccessfully!\n""" 
            #return_value = proc.returncode
            
            #return_value = 0 #KMK
            #if return_value != 0:
                #msg = "\n\n*** ERROR RAISED DURING STEP: %s\n" % e[0] +\ "Command run was:\n %s\n" % e[1] +\ "Command returned exit status: %d\n" % return_value +\ "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
    #            logger.write(msg)
    #            logger.close()
#            raise WorkflowError, msg

    #logger.close()


#def execute_command():



def print_to_stdout(s):
    print s
    
def no_status_updates(s):
    pass

def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append('--%s' % (param_id))
        if param_value != None:
            result.append(param_value)
    return ' '.join(result)


## End  workflow and related functions
