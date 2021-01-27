"""Defines a class that hold the context of a step of the pipeline.
   Specifically, objects create from the Context class carries information
   about  step:
       a) The input files expected
       b) The output files expected
       c) the command itself
"""

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import os
    import re
    import sys
    from shutil import rmtree
    from optparse import make_option
    from os import path, _exit, remove

    from metapathways import general_utils as gutils
    from metapathways import sysutil as sysutils
except:
    print("Cannot load some modules")
    print(traceback.print_exc(10))
    sys.exit(0)

PATHDELIM = sysutils.pathDelim()

class Context:
    """ This class holds the context of a stage """

    def __init__(self):
        self.outputs = {}
        self.outputs1 = {}
        self.inputs = {}
        self.inputs1 = {}
        self.name = None
        self.status = None
        self.commands = []
        self.message = "Message not set"
        pass

    def isOutputAvailable(self):
        return gutils.doFilesExist(self.outputs.values(), gz=True)

    def isInputAvailable(self, errorlogger=None):
        status = True
        status_messages = []
        for key, _file in self.inputs.items():
            if not gutils.doesFileExist(_file) and not gutils.doesFileExist(_file + ".gz"):
                if errorlogger != None:
                    errorlogger.printf("#STEP\t%s\n", self.name)
                    errorlogger.printf("ERROR\tMissing input %s\n", _file)
                status = False
                status_messages.append([key, _file, '[MISSING]'])
            else:
                status_messages.append([key, _file, '[AVAILABLE]'])

        return status, status_messages

    def getMissingList(self, errorlogger=None):
        missingList = []
        status = True
        for file in self.inputs.values():
            if not gutils.doesFileExist(file):
                missingList.append(file)
                if errorlogger != None:
                    errorlogger.printf("ERROR\tMissing input %s\n", file)
                status = False
        return missingList

    def removeOutput(self, errorlogger=None):
        annotationPATT = re.compile(r"annotation_table")
        for item in self.outputs.values():
            if not path.exists(item):
                continue

            if path.isdir(item):
                if annotationPATT.search(item):
                    pass
                else:
                    rmtree(item)
            else:
                remove(item)
