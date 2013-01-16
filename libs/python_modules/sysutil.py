import sys
import os

mswindows = (sys.platform == "win32")

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
 #   if not mswindows:
 #       return commands.getstatusoutput(cmd)
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    if mswindows: 
        cmd = "RMDIR "+ path +" /s /q"
    else:
        cmd = "rm -rf "+path
    result = getstatusoutput(cmd)
    if(output[0]!=0):
        raise RuntimeError(output[1])
