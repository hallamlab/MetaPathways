#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    import re, sys, os, traceback
    from shutil import rmtree
    from os import getenv, makedirs, path, remove
    from operator import itemgetter
    from os.path import abspath, exists, dirname, join, isdir
    from collections import defaultdict
    from optparse import make_option
except:
    print("Cannot load some modules")
    sys.exit(0)
   
class Command():

    inputs = []
    outputs = []
    commands = []
    flags = []

    def __init__(self):
        pass


class Commands():
    """Contains the list of commands for different steps of the pipeline """

    commands = {}


    def __init__(self):
        pass


    def addCommand(self, stepName, command):
        self.commands[stepName] = command
        return True 

if __name__=="__main__":
    pass
