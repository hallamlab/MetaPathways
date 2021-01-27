"""Contains general error message for the MetaPathways """

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import re
    import sys
    import os
    from os import getenv, makedirs, path, remove
except:
    print("Cannot load some modules")
    print(traceback.print_exc(10))
    sys.exit(0)


class Configuration:

    knownErrors = {}

    def __init__(self):
        try:
            self.file = self.initializeErrorMessages()
        except IOError:
            print('ERROR : Cannot  initialize "Configuration" in file ' + sys.argv[0])

    def isKnownError(self, errorid):
        """ Checks if the errorid is known or not"""
        _errorid = str(errorid)
        if _errorid in self.knownErrors:
            return True
        return False

    def getErrorMessage(self, errorid):
        _errorid = str(errorid)
        if self.isKnownError(errorid):
            return self.knownErrors[_errorid]

        return self.knownErrors["0"]

    def initializeErrorMessages(self):
        self.knownErrors = {"0": """ Error occured but the reason is not known"""}


if __name__ == "__main__":
    v = Configuration()
