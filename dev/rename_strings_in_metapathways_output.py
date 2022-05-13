#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re
     import sys, os
     from optparse import OptionParser, OptionGroup
     import optparse
     from glob import glob
except:
     print """ Could not load some user defined  module functions"""
     sys.exit(3)



class CustomOptionParser(optparse.OptionParser):
   def format_epilog(self, formatter):
        return self.epilog


script_name = sys.argv[0]
what_i_do = """This script is used for changing the format of a MetaPathways 2.5 to 3.0. 
As an example consider the following situation: We have a sample named X and we would like to change its files,
 as if it was named Y; and also suppose we have a string S that appears in the files (note that it should somewhat unique otherwise
 the changes may not be meaningful) and we want to change it to some other string that is T. and the folder myfiles/output/X contains the 
 output.  Essentially, it does the following  two tasks to achieve this goal

a) change the sample names, by changing the names of the files (located under the folder myfiles/output/X )  appropriately from old sample name, say  X,  to new sample name, say Y.  Example. the file for the QAd orf are changed from  X.qced.faa to Y.qced.faa, etc 
Note that it does not change the names of folders, only files. 

 b) changes all occurances of the string S to T in every file located under the folder myfiles/output/X. A useful situtation is when when 
 we want to make the orfs names in an output or PGDB from X_23_0 to simply 23_0, where 23 is some contig and 0 is the first ORF there.
 Note that this script does not know about the format of MetaPathways and works in general for any such folder and files with the above rules, 
 so use it with that knolwedge and reckless usage may harm the output 

 Example 1.  Change the sample name from N1 to M100 located in folder N1
 command : python rename_strings_in_metapathways_output.py -f N1 --old_sample_name N1 --new_sample_name M100 

 Example 2.  Change the Contig and  ORF style as from X_d to d  and X_d_e, respectively, where d refers to a contig number and e to an ORF number in contig d
 command : python rename_strings_in_metapathways_output.py -f N1 --old_sample_name N1 --new_sample_name M100 
 python rename_strings_in_metapathways_output.py -f N1 --ostring N1 --nstring "" or 
 python rename_strings_in_metapathways_output.py -f N1 --ostring N1


 Example 3: Achieve the above two tasks together.
 python rename_strings_in_metapathways_output.py -f N1 --old_sample_name N1 --new_sample_name M100 --ostring N1 
 or  rename_strings_in_metapathways_output.py -f N1 --old_sample_name N1 --new_sample_name M100 




"""

usage= script_name + """ -f <folder> --ostring <old_sample_name> --new_sample_name <new_sample_name> 
                           [ --ostring <old_string> ] [ --nstring <new_string> ]"""
parser = OptionParser(usage, description= what_i_do)

parser.add_option( "-f", dest="folder",  
                  help='the folder the sample folder resides , such as, myfolder/output/X,  but not myfolder/output')

parser.add_option( "--ostring", dest="ostring", default = "", 
                  help='old string you want to replace [default is empty]')

parser.add_option( "--nstring", dest="nstring", default = "",  
                  help='new string  [default is empyt string]')

parser.add_option( "--old_sample_name", dest="osample", default = None, 
                  help='old sample name')

parser.add_option( "--new_sample_name", dest="nsample", default = None, 
                  help='old sample name')


def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)



def check_arguments(opts, args):
    if opts.folder == None:
         return False
    return True

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

     return attributes


def  fix_pgdb_input_files(pgdb_folder, pgdbs = []):
     pgdb_list = glob(pgdb_folder + '*cyc/*/input/organism.dat')     

     for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[0]+'\t' + id
     

def getID(line):
     fields = line.split('\t')
     if len(fields)==2:
        return fields[1]
     
def process_organism_file(filel):
     patternsToFix = [ re.compile(r'NAME\tunclassified sequences'), re.compile(r'ABBREV-NAME\tu. sequences') ]
     patternID =  re.compile(r'^ID\t.*')
     try:
         orgfile = open(filel,'r')
     except IOError:
         print "ERROR : Cannot open organism file" + str(filel)
         return 

     lines = orgfile.readlines()
     newlines = []

     needsFixing = False

     id = None
     for line in lines:
         line = line.strip()
         if len(line)==0:
            continue
         flag = False

         result = patternID.search(line)
         if result:   
             id = getID(line)
          
         for patternToFix in patternsToFix:
             result = patternToFix.search(line)
             if result and id:
                 newline = fixLine(line, id)
                 newlines.append(newline)
                 flag= True
                 needsFixing = True

         if flag==False:
            newlines.append(line)

     orgfile.close()
     if needsFixing:
       write_new_file(newlines, filel)


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


def  modify_and_rename_files(filerenames, source, target):


    for file, newname in filerenames.items():
        if len(source)>0 and  source!=target:
            with open(file) as f:
               lines = f.readlines()
               outputfile = open(file + '.tmp','w')
               for line in lines:
                   newline = line.replace(source, target)
                   fprintf(outputfile, "%s", newline)

               outputfile.close();
               os.remove(file)
               os.rename(file+'.tmp', newname)
        else:
           if file!=newname:
              os.rename(file, newname)


def get_files_to_modify(folder):

   all_files = []
   for root, directories, filenames in os.walk( folder +'/'):
     for filename in filenames: 
        all_files.append(os.path.join(root,filename))

   return all_files
   

def  get_new_file_names(files, source, target):

     if source!=None:
       sourceRe = re.compile(source)
     filerename ={}

     for file in files:
         if source!=None and target != None:
           if sourceRe.search(file):
             newname = file.replace(source+'.', target + '.', 2)
             filerename[file] = newname 
         else:
            filerename[file] = file 
     return filerename

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    files = get_files_to_modify(opts.folder)

    filerenames= get_new_file_names(files, opts.osample, opts.nsample)

    filerenames=modify_and_rename_files(filerenames, opts.ostring, opts.nstring)
    

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

