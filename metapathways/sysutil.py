"""Utility functions that are reliant mostly  on sys and os"""

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
    from datetime import date
except:
    print("Cannot load some modules")
    print(traceback.print_exc(10))
    sys.exit(0)

DATEMAP = [
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
]

def genbankDate():

    D = date.today()

    if int(D.day) < 10:
        return "0" + str(D.day) + "-" + DATEMAP[D.month - 1] + "-" + str(D.year)
    else:
        return str(D.day) + "-" + DATEMAP[D.month - 1] + "-" + str(D.year)


def os_type():

    x = sys.platform
    if x:

        hits = re.search(r"darwin", x, re.I)
        if hits:
            return "mac"

        hits = re.search(r"win", x, re.I)
        if hits:
            return "win"

        hits = re.search(r"linux", x, re.I)
        if hits:
            return "linux"


def pathDelim():
    ostype = os_type()
    if ostype == "win":
        return "\\"

    if ostype in ["linux", "mac"]:
        return "/"


def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + " 2>&1", "r")
    text = pipe.read()
    sts = pipe.close()
    if sts is None:
        sts = 0
    if text[-1:] == "\n":
        text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    if mswindows:
        cmd = "RMDIR " + path + " /s /q"
    else:
        cmd = "rm -rf " + path
    result = getstatusoutput(cmd)
    if output[0] != 0:
        raise RuntimeError(output[1])
