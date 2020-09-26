#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

import re
import sys, os
from optparse import OptionParser, OptionGroup
from glob import glob

script_name = sys.argv[0]
usage= script_name + """ --pin <input pf file> --pout <output pf file>  \n --rin <reduced input file> 
                         --rout <reduced output file>
                          This script takes the 0.pf and reduced.txt files to create a new set of 
                          0.pf and  reduced.txt
                     """
parser = OptionParser(usage)

parser.add_option( "--pin", dest="pfin", default =None, 
                  help='the input 0.pf file')

parser.add_option( "--pout", dest="pfout", default = None, 
                  help='the output .pf file')

parser.add_option( "--rin", dest="redin", default = None, 
                  help='reduced input file')

parser.add_option( "--rout", dest="redout", default = None, 
                  help='reduced out file')


parser.add_option( "--level", dest="level", default = '0',  choices=['0', '1'],
                  help='level of equivalennce    \n0 : if ECs and FUNCTIONS are the similar, \n1: if ECs are the similar')

def fprintf(file, fmt, *args):
    file.write(fmt % args)

def printf(fmt, *args):
    sys.stdout.write(fmt % args)



def check_arguments(opts, args):
    if opts.pfin == None:
         return False
    if opts.redin == None:
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
     
def read_pf_file(filel):
     patternID =  re.compile(r'^ID\t.*')

     annotations = []
     ecs= {}
     try:
         orgfile = open(filel,'r')
     except IOError:
         print("ERROR : Cannot open organism file" + str(filel))
         return 
     lines = orgfile.readlines()
     orfids = []
     annotation={}
     products ={}

     for line in lines:
         fields = line.strip().split('\t')
         if len(fields)!=2:
            continue

         if fields[0]=="ID":
             id = fields[1]
             orfids.append(id)
             annotation[id] = []

         annotation[id].append(line.strip())

         if fields[0]=="EC":
            ecs[id] = fields[1]

         if fields[0]=="FUNCTION":
            products[id] = fields[1]

     orgfile.close()

     return orfids, products, annotation, ecs

     
     #for id in orfids:
     #   for line in annotation[id]:
     #      print line

def read_reduced_file(filel):

     try:
         orgfile = open(filel,'r')
     except IOError:
         print("ERROR : Cannot open organism file" + str(filel))
         return 
     lines = orgfile.readlines()
     equivalent={}

     for line in lines:
         fields = line.strip().split('\t')
         if len(fields)!=2:
            continue
         equivalent[fields[0]] =fields[1]

     orgfile.close()

     return equivalent


def write_new_file(lines, output_file):
    
    print("Fixing file " + output_file )
    try:
       outputfile = open(output_file,'w')
       pass
    except IOError:
         print("ERROR :Cannot open output file "  + output_file)
   
    for line in lines:
       fprintf(outputfile, "%s\n", line)

    outputfile.close()


def  modify_files(filerenames, source, target):

    for file, newname in filerenames.iteritems():
        print("<<<<" , file)
        with open(file) as f:
            lines = f.readlines()
            outputfile = open(file + '.tmp','w')
            for line in lines:
                newline = line.replace(source, target)
                fprintf(outputfile, "%s", newline)
            outputfile.close();
            os.remove(file)
            os.rename(file+'.tmp', newname)


def get_files_to_modify(folder):

   all_files = []
   for root, directories, filenames in os.walk( folder +'/'):
     for filename in filenames: 
        all_files.append(os.path.join(root,filename))

   return all_files
   

def  get_new_file_names(files, source, target):
     sourceRe = re.compile(source)
     filerename ={}

     for file in files:
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
       print(usage)
       sys.exit(0)

    orfids, products, annotations, ecs = read_pf_file(opts.pfin)
    equivalent  = read_reduced_file(opts.redin)

    seenec ={}
    seenproducts ={}
    neworfs =[]

    for orfid in orfids:
       if not orfid in products:
          continue

       if not orfid in ecs:
          if products[orfid] in seenproducts:
             equivalent[orfid] = seenproducts[products[orfid]]
          else:
             seenproducts[products[orfid]] = orfid
             neworfs.append(orfid)
          continue


       if ecs[orfid] in seenec: 
          if not orfid in equivalent:
             equivalent[orfid] = seenec[ecs[orfid]]
       else:
          neworfs.append(orfid)
          seenec[ecs[orfid]] = orfid
          
    write_new_pf_file(neworfs, annotations,opts.pfout)

#   print  len(orfids), len(neworfs), len(equivalent)
    print("Original # ORFs   : {}".format(len(orfids)))
    print("Compressed # ORFs : {}".format(len(neworfs)))
    print("Equiv # ORFS      : {}".format(len(equivalent)))
    print("Annotations #     : {}".format(len(annotations)))

    write_new_reduced_file(equivalent, opts.redout)


def write_new_reduced_file(equivalent, filename):
    outputfile = open(filename,'w')
    for orfid in equivalent: 
       if orfid !=  equivalent[orfid]:
          fprintf(outputfile, "%s\t%s\n",orfid, equivalent[orfid])
    outputfile.close()

def write_new_pf_file(neworfs, annotations, filename):
    outputfile = open(filename,'w')
    for orfid in neworfs: 
       for line in annotations[orfid]:
           fprintf(outputfile, "%s\n",line)
       fprintf(outputfile, "//\n")
    outputfile.close()

    #print len(equivalent.keys()), len(orfids), len(ecs)
    #for ec in ecs:
    #   print ecs[ec]


# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

