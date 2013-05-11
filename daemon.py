#!/usr/bin/python
# File created on 27 Jan 2012.
from __future__ import division

__author__ = "Kishori M Konwar, Niels W Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = [""]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar, Niels W Hanson"
__status__ = "Release"

import subprocess
import sys
from os import remove, makedirs, sys, listdir, environ, path
import re 
import inspect
from commands import getstatusoutput
from optparse import OptionParser
import shutil 
import traceback
from glob import glob



#config = load_config()
metapaths_config = """template_config.txt""";

script_info={}
script_info['script_usage'] = []
usage= """./daemon  --home-dir cwd   """
parser = OptionParser(usage)
parser.add_option("--home-dir", dest="home_dir", default ='\'\'',
                  help='home dir [REQUIRED]')
parser.add_option("--does-sample-dir-exist", dest="does_sample_dir_exist", default='',
                  help='does sample dir exist')
parser.add_option("--does-file-exist", dest="does_file_exist", default='',
                  help='does file dir exist')

parser.add_option("--does-file-patt-exist", dest="does_file_patt_exist", default='',
                  help='does file with the pattern exist')

parser.add_option("--create-sample-dir", dest="create_sample_dir",  default='', 
                  help='create sample dir')

parser.add_option("--remove-sample-dir", dest="remove_sample_dir",  default='', 
                  help='remove sample dir')

parser.add_option("--number-of-sequences-in-file", dest="number_of_sequences_in_file",  default='', 
                  help='count the number of sequences in file')

parser.add_option("--split-into-batches", dest="split_into_batches",  default='', 
                  help='count the number of sequences in file')

parser.add_option("--format-database", dest="format_database",  default='', 
                  help='formats the database')

parser.add_option("--batch-size", dest="batch_size",  default=500, 
                  help='batch size')

parser.add_option("--os-type", dest="os_type", action= 'store_true',  default=False,   
                  help='return OS type')

parser.add_option("--cpu-type", dest="cpu_type",  default='',   
                  help='return CPU type')

parser.add_option("--algorithm", dest="algorithm", choices = ['BLAST', 'LAST'], default = "BLAST",
                  help='the algorithm used for computing homology [DEFAULT: BLAST]')


parser.add_option("--submit-job", dest="submit_job",  default='',   
                  help='submit job')

parser.add_option("--memory", dest="memory",  default='10gb',   
                  help='memory size request')

parser.add_option("--walltime", dest="walltime",  default='10:00:00',   
                  help='wall time request')

parser.add_option("--database-file", dest="database_file",  default='',   
                  help='database file name')

parser.add_option("--database-files", dest="database_files",  default=[], action='append',   
                  help='database file names')

parser.add_option("--dbname", dest="dbname",  default='',   
                  help='dbname')

parser.add_option("--dbnames", dest="dbnames",  default=[], action='append',   
                  help='dbnames, an array')

parser.add_option("--is-complete", dest="is_complete",  default='',   
                  help='is consolidate')

parser.add_option("--consolidate", dest="consolidate",  default='',   
                  help='consolidate with dbname and sample name')

parser.add_option("--get-number-of-running-jobs", dest="get_number_of_running_jobs", default='', 
                  help='get the number of running jobs')

parser.add_option("--get-number-of-samples", dest="get_number_of_samples", default='', 
                  help='get the number of samples')

parser.add_option("--get-number-of-completed", dest="get_number_of_completed", default='', 
                  help='get the number of completed  samples')


def fprintf(file, fmt, *args):
    file.write(fmt % args)
  
def printf(fmt, *args):
     sys.stdout.write(fmt % args)

def get_sequence_name(line):
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     #print name
     return name

def get_number_of_completed(sample_name):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_list(namePrefix + 'completed.txt', samples_dictionary)
     print str(len(samples_dictionary))

def get_number_of_samples(sample_name):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_dictionary, col=1)
     print str(len(samples_dictionary))
     return len(samples_dictionary)

def _get_number_of_lines_in_file(filename):
    try:
       listfile = open(filename, 'r')
       lines = listfile.readlines()
       listfile.close()
       return len(lines)
    except:
       return 0


def  get_number_of_running_jobs(login):
     args = [ 'qstat', '-u', login]
     p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     result = p.communicate() 
     lines = result[0].strip().split('\n')
     num_running_jobs = 0 
     loginPattern = re.compile( login )

     for line in lines:
        if loginPattern.match(line): 
           num_runnning_jobs += 1 
         
     if num_running_jobs==0:   
        num_running_jobs = str(len(lines)-4) 

     print str(num_running_jobs)


def  remove_sample_dir(sample_dir):
     files = glob(sample_dir + '*')
     for  f in files:
        remove(f)
     if path.exists(sample_dir):
        shutil.rmtree(sample_dir)

def create_the_sequence_files(sequence_file_name, samplefilename, sample_name, database_file, dbname,  size):
     try:
        sequencefile = open(sequence_file_name, 'r')
     except IOError:
        print "Cannot read file " + sequence_file_name + " !"

     sequence_lines = sequencefile.readlines()

     sequencefile.close()
     fragments= []
     name=""

     seq_dictionary={}
     seq_beg_pattern = re.compile(">")
     for line in sequence_lines:
        line = line.strip()
        if seq_beg_pattern.search(line):
          if len(name) > 0:
             sequence=''.join(fragments)
             seq_dictionary[name]=sequence
             fragments = []
          name=get_sequence_name(line)
        else:
          fragments.append(line)

     if len(name) > 0:
        sequence=''.join(fragments)
        seq_dictionary[name]=sequence


     samplefile = open(samplefilename , 'w')

     count =0
     filecount=0
     smallfilename = sample_name +'_' + str(filecount) +'.faa'
     fprintf(samplefile,'%s_%s_%s\t%s\t%s\t%s\n',dbname, database_file, smallfilename, dbname,\
             database_file, smallfilename)

     smallfile = open(smallfilename , 'w')
     for name in seq_dictionary:
       if count %size==0  and count>0:
         smallfile.close()
         filecount+=1
         smallfilename = sample_name +'_' + str(filecount) +'.faa'
         smallfile = open(smallfilename , 'w')
         fprintf(samplefile,'%s_%s_%s\t%s\t%s\t%s\n',dbname, database_file, smallfilename,\
                 dbname, database_file, smallfilename)

       fprintf(smallfile,'>%s\n',name)
       fprintf(smallfile,'%s\n',seq_dictionary[name])
       count+=1

     smallfile.close()
     samplefile.close()

def already_split_for_dbname(submittedfilename, dbname): 
     dbnames_dictionary={}
     read_one_column(submittedfilename, dbnames_dictionary, col=1)
     if dbname in dbnames_dictionary:
        return True
     else:
        return False


def split_into_batches(sequence_file_name, database_files,  dbnames, size):
     sample_name = re.sub(r'[.]qced[.]faa','',sequence_file_name)
     submittedfilename = sample_name +'_submitted.txt'
     completedfilename = sample_name +'_completed.txt'
     samplefilename = sample_name +'_samples.txt'

     if not path.exists(samplefilename):
        if len(dbnames) > 0 and len(database_files)>0:
           create_the_sequence_files(sequence_file_name, samplefilename, sample_name, database_files[0], dbnames[0],  size)
           submittedfile = open(submittedfilename , 'w')
           submittedfile.close()
           completedfile = open(completedfilename , 'w')
           completedfile.close()


     if not path.exists(submittedfilename): 
        submittedfile = open(submittedfilename , 'w')
        submittedfile.close()

     if not path.exists(completedfilename):
        completedfile = open(completedfilename , 'w')
        completedfile.close()


     for database_file, dbname in zip(database_files, dbnames):
         if already_split_for_dbname(samplefilename, dbname): 
            #print "already split for " + dbname
            continue

         samples_filename_dictionary={}
         
         read_one_column(samplefilename, samples_filename_dictionary, col=3)

         samplefile = open(samplefilename , 'a')
         for smallfilename in samples_filename_dictionary:
#            print smallfilename
            fprintf(samplefile,'%s_%s_%s\t%s\t%s\t%s\n',dbname, database_file, smallfilename, dbname, database_file, smallfilename)
            #printf('%s_%s_%s\t%s\t%s\t%s\n',dbname, database_file, smallfilename, dbname, database_file, smallfilename)
         samplefile.close()


     return _get_number_of_lines_in_file(samplefilename)

def  read_one_column(listfilename, dictionary, col=0) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
          dictionary[fields[col]] = True
    listfile.close()
  except:
    traceback.print_exception()

# col begin with 0
def  read_list(listfilename, dictionary, col=1) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col:
        dictionary[fields[0]] = fields[col]
    listfile.close()
  except:
    traceback.print_exception()

def  read_list_reverse(listfilename, dictionary) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields)==1:
          dictionary[fields[0]] = 1
       if len(fields)==2:
          dictionary[fields[1]] = fields[0]
    listfile.close()
  except:
    traceback.print_exception()
         
def add_to_listfile(listfilename, key, id ):
    listfile = open(listfilename, 'a')
    fprintf(listfile, "%s\t%s\n", key, id);
    listfile.close()

def read_completed_task(statsDir, jobid_dictionary):
    files = glob(statsDir +'*')
    for file in files:
       shortfile = re.sub(r'^.*/','',file)
       jobid_dictionary[shortfile] = 1
       #print 'a' + file 
       remove(file)
    #print jobid_dictionary
    return

def format_database(database, algorithm):
     if algorithm=='LAST':
       dbformatter = 'MetaPathways/executables/lastdb' 
       args = [ dbformatter, '-p', '-c',  database, database ]

     if algorithm=='BLAST':
       dbformatter = 'MetaPathways/executables/formatdb' 
       args = [ dbformatter, '-p', 'T', '-i', database ]

     p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
     result = p.communicate()
     if result[1].strip()=='':
        print 'yes'
        return True
     else:
        print 'no'
        return False


def  _create_dictionary_of_arrays(listfilename, dictionary, col1=1, col2=2) :
  try:
    listfile = open(listfilename, 'r')
    lines = listfile.readlines()
    for line in lines:
       fields = [ x.strip() for x in line.strip().split('\t') ]
       if len(fields) > col1 and len(fields) > col2:
           if not fields[col1] in  dictionary:
              dictionary[fields[col1]] = []
           dictionary[fields[col1]].append(fields[col2])
    listfile.close()
  except:
    traceback.print_exception()

def  append_file_content_from_to(sourcefilename, targetfilename):
    sourcefile = open(sourcefilename, 'r')
    targetfile = open(targetfilename, 'a')
    sourcelines = sourcefile.readlines()
    for line in sourcelines:
        fprintf(targetfile, "%s\n", line.strip())
    sourcefile.close()
    targetfile.close()


def consolidate(sample_name, dbname):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     by_dbnames={}
     _create_dictionary_of_arrays(namePrefix + 'samples.txt', by_dbnames, col1=1, col2=3)
     consolidatedfile = 'MetaPathways/' + sample_name + '/' + sample_name +'.' + dbname + '.blastout'
     if path.exists(consolidatedfile):
       remove(consolidatedfile)
     for smallfile in by_dbnames[dbname]:
        #print smallfile + ' ' + consolidatedfile
        append_file_content_from_to(smallfile + '.' + dbname + '.blastout', consolidatedfile)

     print 'yes'
     return True

def is_complete(sample_name, dbname):
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_dictionary, col=1)
     completed_dictionary={}
     read_list(namePrefix + 'completed.txt', completed_dictionary, col=1)
     if len(samples_dictionary) == len(completed_dictionary):
         print 'yes'
         return True
     else:
         print 'no'
         return False


def submit_job(sample_name, mem, walltime, algorithm):

   validWallTime = re.compile(r'[0-9]+:[0-9]+:[0-9]+')
   validMem = re.compile(r'[0-9]+gb')

   if not validWallTime.match(walltime):  
      walltime = '11:00:00'
   if not validMem.match(mem):  
      mem = '11gb'

   try:
     namePrefix = 'MetaPathways/' + sample_name + '/' + sample_name +'_' 
     samples_filename_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_filename_dictionary, col=3)
     samples_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_dictionary, col=2)
     samples_dbname_dictionary={}
     read_list(namePrefix + 'samples.txt', samples_dbname_dictionary, col=1)

     submitted_dictionary={}
     read_list(namePrefix + 'submitted.txt', submitted_dictionary, col=1)
     completed_dictionary={}
     read_list(namePrefix + 'completed.txt', completed_dictionary, col=1)

     for key in samples_dictionary:
       if not key  in submitted_dictionary:

         commandfile = open('command_job_qsub.txt', 'w')

         databasefile = samples_dictionary[key] 
         dbname = samples_dbname_dictionary[key] 
         seq_file_name = samples_filename_dictionary[key] 
          
         if algorithm=='BLAST': 
              command = 'MetaPathways/executables/blastp -num_threads 1  -max_target_seqs 5  -outfmt 6' +\
                        ' -db ' +  ('MetaPathways/databases/' + databasefile)  +\
                        ' -query ' + seq_file_name  + ' -evalue 0.000001 '+\
                        ' -out '  +  (seq_file_name + "." + dbname +".blastout")

         if algorithm=='LAST':
              command = 'MetaPathways/executables/lastal' +\
                        ' -o '  +  (seq_file_name + "." + dbname +".blastout") +\
                        ' -f 0 '  +  ('MetaPathways/databases/' + databasefile)  +\
                        ' ' + seq_file_name



         #fprintf(commandfile, "echo %s\n",key)
         fprintf(commandfile, "%s\n",command)

         commandfile.close()
         args = [ 'qsub', '-l', 'walltime='+walltime, '-m',  'ea',  '-j',  'eo',  '-e',
             'MetaPathways/' + sample_name +'/.qstatdir/$PBS_JOBID',  '-l',  'procs=1',\
             '-l', 'mem='+ mem, 'command_job_qsub.txt']
         p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
         result = p.communicate()
         if result[1].strip()=='':
            add_to_listfile(namePrefix + 'submitted.txt', key, result[0].strip())
         break
   except:
       traceback.print_exception()

   reverse_submitted_dictionary={}
   read_list_reverse(namePrefix + 'submitted.txt', reverse_submitted_dictionary)
   #print reverse_submitted_dictionary
  

   statsDir = 'MetaPathways/' + sample_name + '/.qstatdir/'
   jobid_dictionary={}
   read_completed_task(statsDir, jobid_dictionary)
   #print jobid_dictionary
   for key in submitted_dictionary:
  #    print submitted_dictionary[key] 
      if submitted_dictionary[key]  in jobid_dictionary:
   #       print key + '\t' + submitted_dictionary[key] 
         if not key  in completed_dictionary:
              add_to_listfile(namePrefix + 'completed.txt', key, submitted_dictionary[key])


def os_type():
    try:
       from platform import machine, platform
    except:
       print "unknown"
       return
    
    if platform():
       print platform()

def number_of_sequences_in_file(sequence_file_name):
     try:
        sequencefile = open(sequence_file_name, 'r')
     except IOError:
        print "Cannot read file " + sequence_file_name + " !"
     sequence_lines = sequencefile.readlines()
     sequencefile.close()
     name=""
     count =0
     seq_beg_pattern = re.compile(">")
     for line in sequence_lines:
        line = line.strip()
        if seq_beg_pattern.search(line):
          count+=1
     return count


# checks if the supplied arguments are adequate
def isValid(opts, args):
    if (opts.home_dir==None):
       return True
    else:
       return False


# in the folder and then delete the folder too
def removeSampleDir( origFolderName):
    folderName = origFolderName + '/*'
    files = glob(folderName)
    for  f in files:
       remove(f)
    if path.exists(origFolderName):
      shutil.rmtree(origFolderName)


def doesFileExist(file):
    if path.exists(file):
       return True
    else:
       return False

def doesFilePattExist(filepatt):
    files = glob(filepatt)
    if len(files)>0:
       return True
    else:
       return False

def doesSampleDirExists(sampledir):
    dir = sampledir
    if path.exists(dir):
       return True
    else:
       return False


def extendPath( folders):
    newdir = ''
    for folder in folders:
       if len(folder)!=0:
        newdir =  newdir + folder  + '/'
    return newdir

def main(argv):
    (opts, args) = parser.parse_args()
    if isValid(opts, args):
       print usage
       sys.exit(0)


    # initialize the input directory or file
    if len(opts.does_sample_dir_exist):
        sampledir = extendPath([opts.home_dir, opts.does_sample_dir_exist])
        if not doesSampleDirExists(sampledir):
           print 'no'
        else:
           print 'yes'

    # create the sample directory
    if len(opts.create_sample_dir)>0:
        try:
            sampledir = extendPath([opts.home_dir, opts.create_sample_dir])
            makedirs(sampledir)
            print 'yes'
        except:
            print 'no'

    # remove the sample directory
    if len(opts.remove_sample_dir)>0:
        try:
            sampledir = extendPath([opts.home_dir, opts.remove_sample_dir])
            remove_sample_dir(sampledir)
            print 'yes'
        except:
            print 'no'

    # create the sample directory
    if len(opts.does_file_exist)>0:
        try:
            #file = extendPath([ opts.home_dir, opts.create_sample_dir,opts.dies_file_exist] )
            if doesFileExist(opts.does_file_exist): 
               print 'yes'
            else:
               print 'no'
        except:
            print 'no'
     
    # does file with the pattern exist
    if len(opts.does_file_patt_exist)>0:
        try:
            if doesFilePattExist(opts.does_file_patt_exist): 
               print 'yes'
            else:
               print 'no'
        except:
            print 'no'

    if len(opts.number_of_sequences_in_file) > 0:
        count = number_of_sequences_in_file(opts.number_of_sequences_in_file)
        print str(count)
    #print opts


    if len(opts.split_into_batches) > 0:
       try:
          batch_count = split_into_batches(opts.split_into_batches, opts.database_files, opts.dbnames,  int(opts.batch_size))
          print batch_count
       except:
          return 0

    if len(opts.submit_job) > 0:
       try:
          submit_job(opts.submit_job, opts.memory, opts.walltime, opts.algorithm)
       except:
          return 0

    if opts.os_type==True:
       try:
          os_type()
       except:
          return 0

    # format database
    if len(opts.format_database) > 0:
       try:
          format_database(opts.format_database, opts.algorithm)
       except:
          return 0

    #get the number of running jobs
    if len(opts.get_number_of_running_jobs)>0:
       try:
          get_number_of_running_jobs(opts.get_number_of_running_jobs)
       except:
          return 0

    #get number of samples 
    if len(opts.get_number_of_samples)>0:
       try:
          get_number_of_samples(opts.get_number_of_samples)
       except:
          return 0

    #get number of completed samples 
    if len(opts.get_number_of_completed)>0:
       try:
          get_number_of_completed(opts.get_number_of_completed)
       except:
          return 0

    # check is completed 
    if len(opts.is_complete) :
      try:
         is_complete(opts.is_complete, opts.dbname)
      except:
         return 0

    # check is completed 
    if len(opts.consolidate)>0 and len(opts.dbname)>0: 
      try:
         consolidate(opts.consolidate, opts.dbname)
      except:
         return 0

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])    

