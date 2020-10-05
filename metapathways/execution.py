"""This file contains the metapaths orkflow logic that
   executes and reports status
"""

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import sys
    import os
    import re

    from subprocess import Popen, PIPE, STDOUT
    from os import makedirs, listdir, _exit
    from glob import glob
    from optparse import OptionParser
    from os.path import split, splitext, join, dirname, abspath
    from datetime import datetime

    import sysutil as sysutils
    import general_utils as gutils
    import scripts as python_scripts

except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)

def execute_pipeline_stage(pipeline_command, \
     extra_command=None, errorlogger=None, runstatslogger=None):

    argv = [x.strip() for x in pipeline_command.split()]

    funcname = re.sub(r".py$", "", argv[0])
    funcname = re.sub(r"^.*/", "", funcname)
    args = argv[1:]
    if hasattr(python_scripts, funcname):
        methodtocall = getattr(getattr(python_scripts, funcname), funcname)
        if extra_command == None:
            result = methodtocall(
                args, errorlogger=errorlogger, runstatslogger=runstatslogger)
        else:
            result = methodtocall(
                args,
                errorlogger=errorlogger,
                extra_command=extra_command,
                runstatslogger=runstatslogger,
            )
    else:
        result = sysutils.getstatusoutput(pipeline_command)
    return result


def printMissingList(missingList):
    gutils.eprintf("MISSING INPUT LIST:\n")
    for missingItem in missingList:
        gutils.eprintf("     %s\n", missingItem)

def execute_tasks(s, verbose=False, block=0):
    """Run list of commands, one after another """
    # logger.write("Executing commands.\n\n")
    contextBlocks = s.getContextBlocks()

    contextBlock = contextBlocks[block]

    for c in contextBlock:
        if c.status == "stop":
            print("Stopping!")
            s.stepslogger.write("%s\t%s\n" % (c.name, "STOPPED"))
            return (0, "")

        if verbose:
            gutils.eprintf("\n\n\nEXECUTED COMMAND : %s\n", ", ".join(c.commands))

        gutils.eprintf("%s" % (c.message))

        if c.status in ["redo"]:
            c.removeOutput(s)
            status, status_messages = c.isInputAvailable(errorlogger=s.errorlogger)
            if status:
                s.stepslogger.write("%s\t%s\n" % (c.name, "RUNNING"))
                try:
                    result = execute(s, c)
                    if result != None and  result[0] != 0:
                        print('ERROR: ', result)

                except:
                    s.errorlogger.printf("ERROR\t%s\n", result[1])
                    gutils.eprintf(traceback.print_exc(10))
                    result[0] = 1

                if result[0] == 0:
                    gutils.eprintf("..... Redo Success!\n")
                    s.stepslogger.write("%s\t%s\n" % (c.name, "SUCCESS"))
                else:
                    gutils.eprintf("..... Failed!\n")
                    # eprintf('%s result \n',  result )
                    s.stepslogger.write("%s\t%s\n" % (c.name, "FAILED"))
            else:
                gutils.eprintf("..... Skipping [NO INPUT]!\n")
                if verbose:
                    missingList = c.getMissingList(errorlogger=s.errorlogger)
                    printMissingList(missingList)
                s.stepslogger.write("%s\t%s\n" % (c.name, "MISSING_INPUT"))

            if verbose:
                for status_message in status_messages:
                   gutils.eprintf("\t{}\t{}\t{}\n".format(status_message[0], status_message[1], status_message[2]))
        elif c.status in ["yes"]:
            if not c.isOutputAvailable():
                status, status_messages = c.isInputAvailable(errorlogger=s.errorlogger)
                if status:
                    s.stepslogger.write("%s\t%s\n" % (c.name, "RUNNING"))
                    try:
                        result = execute(s, c)
                        if result != None and  result[0] != 0:
                           print('ERROR: ', result)
                    except:
                        s.errorlogger.printf("ERROR\t%s\n", result[1])
                        gutils.eprintf(traceback.print_exc(10))
                        result[0] = 1

                    if result[0] == 0:
                        gutils.eprintf("..... Success!\n")
                        s.stepslogger.write("%s\t%s\n" % (c.name, "SUCCESS"))
                    else:
                        gutils.eprintf("..... Failed!\n")
                        s.stepslogger.write("%s\t%s\n" % (c.name, "FAILED"))
                else:
                    gutils.eprintf("..... Skipping [NO INPUT]!\n")
                    if verbose:
                        missingList = c.getMissingList(errorlogger=s.errorlogger)
                        printMissingList(missingList)

                    s.stepslogger.write("%s\t%s\n" % (c.name, "SKIPPED"))

                if verbose:
                   for status_message in status_messages:
                       gutils.eprintf("\t{}\t{}\t{}\n".format(status_message[0], status_message[1], status_message[2]))

            else:
                gutils.eprintf("..... Already Computed!\n")
                s.stepslogger.write("%s\t%s\n" % (c.name, "ALREADY_COMPUTED"))


        elif c.status in ["skip"]:
            gutils.eprintf("..... Skipping!\n")
            s.stepslogger.write("%s\t%s\n" % (c.name, "SKIPPED"))

def execute(s, c):
    result = [1, "Error while executing " + c.name]
    try:
        if len(c.commands) == 2:
            result = execute_pipeline_stage(
                c.commands[0],
                extra_command=c.commands[1],
                errorlogger=s.errorlogger,
                runstatslogger=s.runstatslogger,
            )
        else:
            result = execute_pipeline_stage(
                c.commands[0],
                errorlogger=s.errorlogger,
                runstatslogger=s.runstatslogger,
            )
    except:
        print('ERROR :', traceback.print_exc(10))
        result = [1, "Error while executing " + c.name]

        pass
    return result

def get_params_str(params):
    result = []
    for param_id, param_value in params.items():
        result.append("--%s" % (param_id))
        if param_value != None:
            result.append(param_value)
    return " ".join(result)


## End  workflow and related functions
