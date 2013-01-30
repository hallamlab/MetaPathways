#!/usr/bin/python
# File created on 27 Jan 2012.

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2010, The metapaths Project"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

from optparse import make_option
from os import makedirs, path, listdir, remove
import sys
from sysutil import getstatusoutput
#from metapaths_utils  import parse_command_line_parameters
from parse  import parse_metapaths_parameters
from metapaths_pipeline import print_commands, call_commands_serially, WorkflowLogger, generate_log_fp, generate_steps_log_fp
import shutil
import re
from metapaths_utils import fprintf, printf, remove_existing_pgdb
import  errno
from glob import glob
from datetime import date

def copyFile(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

def dry_run_status( commands ):
    for command in commands:
        printf("%s", command[0])
        if command[4] == True:
           printf("%s", " Required")
        else:
           printf("%s", " Not Required")
    printf("\n")


def get_refdb_name( dbstring ):
    dbstring = dbstring.rstrip()
    dbstring = dbstring.lstrip()
    dbstring = dbstring.lower() 
    if re.search("kegg",dbstring):
        dbname = "kegg"
        return dbname
    if re.search("cog",dbstring):
        dbname = "cog"
        return dbname
    if re.search("metacyc",dbstring):
        dbname = "metacyc"
        return dbname
    if re.search("refseq",dbstring):
        dbname = "refseq"
        return dbname

    if re.search("InnateDB_human",dbstring):
        dbname = "innatehuman"
        return dbname
     
    if re.search("InnateDB_mus",dbstring):
        dbname = "innatemus"
        return dbname

    if re.search("silva",dbstring):
        dbname = "silva"
        return dbname
    if re.search("greengenes",dbstring):
        dbname="greengenes"
        return dbname

    dbname= re.sub("[^a-zA-Z].*","", dbstring)
    return dbname


def format_db(formatdb_executable, seqType, refdb_sequence_file):
     cmd='%s -p %s -i %s' %(formatdb_executable, seqType, refdb_sequence_file)
     result= getstatusoutput(cmd)
     if result[0]==0:
        return True
     else:
        return False



def create_quality_check_cmd(min_length, log, input, output, mapping, config_settings):
    cmd = "%s  --min_length %d --log_file %s  -i %s -o  %s -m %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PREPROCESS_FASTA']), min_length, log, input,  output, mapping) 
    return cmd


def create_orf_prediction_cmd(options, log, input, output, config_settings):
    cmd = "%s %s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['ORF_PREDICTION']), options, input, output) 
    return cmd

# convert an input gbk file to fna faa and gff file
def  convert_gbk_to_fna_faa_gff(input_gbk, output_fna, output_faa, output_gff, config_settings):
    cmd = "%s  -g %s --output-fna %s --output-faa %s --output-gff %s" %((config_settings['METAPATHWAYS_PATH'] \
                 + config_settings['GBK_TO_FNA_FAA_GFF']), input_gbk, output_fna, output_faa, output_gff) 
    return cmd

# convert an input gff file to fna faa and gff file
def  convert_gff_to_fna_faa_gff(inputs, outputs,  config_settings):
    cmd = "%s " %(config_settings['METAPATHWAYS_PATH']+ config_settings['GFF_TO_FNA_FAA_GFF'])
    for  source, target in zip(inputs, outputs):
       cmd += ' --source ' + source + ' --target ' + target
    return cmd


#converts the orfs into 
#def create_orf_converter_cmd(options, input, output, config_settings):
#    cmd = "%s  %s  %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['ORF_CONVERTER']), options,  input) 
#    return cmd

def create_aa_orf_sequences_cmd(input_gff, input_fasta, output_fna, output_faa, output_gff,config_settings):
    cmd = "%s -g  %s  -n %s --output_nuc %s --output_amino %s --output_gff %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GFF_TO_FASTA']), input_gff, input_fasta, output_fna, output_faa, output_gff) 
    return cmd

# create genbank files
def create_genebank_file_cmd(output, input_orf, input_fas, sample_name, config_settings) :
    cmd = "%s  --out %s --orphelia %s --fasta %s  --sample-name %s"\
           %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_GENBANK']), output, input_orf, input_fas, sample_name) 
    return cmd

# extract from the genbank files a fasta file
def create_genebank_2_fasta_cmd(options, input, output, config_settings):
    cmd = "%s %s --o %s %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_TO_FASTA']), options, output, input) 
    return cmd

def create_create_filtered_amino_acid_sequences_cmd(input, output, logfile, min_length, config_settings):
    cmd = "%s  --min_length %d --log_file %s  -i %s -o  %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PREPROCESS_FASTA']), min_length, logfile, input,  output) 
    return cmd 

#compute refscores
def create_refscores_compute_cmd(input, output, config_settings):
    cmd = "%s  -F %s -B %s -o %s -i %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['COMPUTE_REFSCORE']), \
           (config_settings['METAPATHWAYS_PATH'] + config_settings['FORMATDB_EXECUTABLE']), \
           (config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTP_EXECUTABLE']), \
           output, input) 
    return cmd

def formatted_db_exists(dbname):
    fileList = glob(dbname) 
    if len(fileList) > 0:
       return True
    else: 
       return False

def check_if_db_exists(dbname):
    if path.exists(dbname):
       return True
    else: 
       return False


# Makes sure that the ref database is formatted for blasting 
def check_an_format_refdb(dbname, seqType,  config_settings): 
    suffix=[]
    if seqType=='T':
      suffix.append('*.phr')
      suffix.append('*.psq')
      suffix.append('*.pin')
    else:
      suffix.append('*.nhr')
      suffix.append('*.nsq')
      suffix.append('*.nin')
    
    fullRefDbName = config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['REFDBS'] +"/" + dbname


    if not (formatted_db_exists(fullRefDbName+suffix[0]) and formatted_db_exists(fullRefDbName + suffix[1])\
             and formatted_db_exists(fullRefDbName + suffix[2]) ):
        print "WARNING : You do not seem to have Database " + dbname + " formatted!"
        if  check_if_db_exists(fullRefDbName):
            print "          Found raw sequences for  Database " + dbname + "!"
            print "          Trying to format on the fly ....!"
            result =format_db(config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['FORMATDB_EXECUTABLE'], seqType, fullRefDbName)
            if result ==True:
                print "          Formatting successful!"
                return 
            else:
                print "          Formatting failed! Please consider formatting manually or do not try to annotated with this database!"
                sys.exit(0)
        print "ERROR : You do not even have the raw sequence for Database " + dbname + " to format!"
        print "        Please put the appropriate files in folder \"blastDB\""
        sys.exit(0)
  

#    fullRefDbMapName = config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['REFDBS'] +"/" + dbname + '-names.txt'
#    if not doFilesExist( [fullRefDbMapName ] ):
#        print 'ssss' + fullRefDbMapName
#        print 'no map file' 
#        sys.exit(0)

   

def  make_sure_map_file_exists(dbmapfile):
    if not doFilesExist( [dbmapfile ] ):
         print 'WARNING: ' + 'Creating the database map file'
         fullRefDbName = re.sub(r'-names.txt','',dbmapfile)
         mapfile = open(dbmapfile,'w')
         fullRefDbFile = open(fullRefDbName,'r')
         for line in fullRefDbFile:
             if re.match(r'>', line):
                 fprintf(mapfile, "%s\n",line.strip())
         mapfile.close()
         fullRefDbFile.close()


#        print 'ssss' + fullRefDbMapName
#        print 'no map file' 
#        sys.exit(0)

# creates the command to blast the sample orf sequences against the reference databases for the 
#purpose of annotation
def create_blastp_against_refdb_cmd(input, output, output_dir, sample_name, dbfilename, config_settings, config_params,  run_command): 
    max_evalue = float(get_parameter(config_params, 'annotation', 'max_evalue', default=0.000001))
    system =    get_parameter(config_params,  'metapaths_steps', 'BLAST_REFDB', default='yes')
    max_hits =    get_parameter(config_params,  'annotation', 'max_hits', default=5)


    dbname = get_refdb_name(dbfilename);
    if run_command:
       check_an_format_refdb(dbfilename, 'T',  config_settings)

    if system=='grid':
       batch_size = get_parameter(config_params,  'grid_engine', 'batch_size', default=500)
       max_concurrent_batches = get_parameter(config_params,  'grid_engine', 'max_concurrent_batches', default=500)
       user = get_parameter(config_params,  'grid_engine', 'user', default='')
       server = get_parameter(config_params,  'grid_engine', 'server', default='')

       cmd = {'sample_name':output_dir,  'dbnames':[dbname], 'database_files':[dbfilename],\
               'run_type':'overlay', 'batch_size':batch_size,'max_parallel_jobs':max_concurrent_batches,\
               'user':user, 'server':server }  
    else:
       cmd="%s -num_threads 16  -max_target_seqs %s  -outfmt 6  -db %s -query  %s -evalue  %f  -out %s"\
         %((config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTP_EXECUTABLE']), max_hits,\
          config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['REFDBS'] +"/" + dbfilename, input, max_evalue, output) 
    
    return cmd


# Command for parsing the blast flie snd create the parse blast files
#     input -- blastoutput
#     output -- parsed files
#     refscorefile   -- refscore file
#     min_bsr   -- minimum bsr ratio for accepting into annotation
#     max_evalue  -- max evalue
#     min_score   -- min score
#     min_length  -- min_length in amino acids, typcically 100 amino acids should be minimum
#     dbmane      -- name of the refrence databases
def create_parse_blast_cmd(input, refscorefile, dbname, dbmapfile,  config_settings, config_params): 
    min_bsr = float(get_parameter(config_params, 'annotation', 'min_bsr', default=0.4))
    min_score = float(get_parameter(config_params, 'annotation', 'min_score', default=0.0))
    min_length = float(get_parameter(config_params, 'annotation', 'min_length', default=100))
    max_evalue = float(get_parameter(config_params, 'annotation', 'max_evalue', default=1000))
    
     
    cmd="%s -d %s  -b %s -m %s  -r  %s  --min_bsr %f  --min_score %f --min_length %f --max_evalue %f"\
         %((config_settings['METAPATHWAYS_PATH'] + config_settings['PARSE_BLAST']), dbname, input, dbmapfile,  refscorefile, min_bsr, min_score, min_length, max_evalue);

    #if dbname=='refseq':
    #   cmd = cmd + " --remove_tax"

    #if not path.exists( dbmapfile ):
    #    print "ERROR: You do not have the required  annotation map file "
        #print "      " + dbmapfile
    #    sys.exit(0)

    return cmd

# this function creates the command that is required to make the annotated genbank file
def create_annotate_genebank_cmd(sample_name, input_unannotated_gff, output_annotated_gff,\
           blast_results_dir,  options, refdbs, output_comparative_annotation, config_settings):
    cmd="%s --input_gff  %s -o %s  %s --output-comparative-annotation %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['ANNOTATE']),input_unannotated_gff, output_annotated_gff,  options, output_comparative_annotation) 

    for refdb in refdbs:
        cmd = cmd + " -b " + blast_results_dir + "/" + sample_name + "." + refdb+ ".blastout.parsed.txt -d " + refdb + " -w 1 "

    return cmd

def create_genbank_ptinput_sequin_cmd(input_annotated_gff, nucleotide_fasta, amino_fasta, outputs, config_settings):
    
    cmd="%s -g %s -n %s -p %s " %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_FILE']), input_annotated_gff, nucleotide_fasta, amino_fasta) ; 

    if 'gbk' in outputs:
       cmd += ' --out-gbk ' + outputs['gbk']

    if 'ptinput' in outputs:
       cmd += ' --out-ptinput ' + outputs['ptinput']


    return cmd

def  create_report_files_cmd(dbs, input_dir, input_annotated_gff,  sample_name, output_dir, config_settings):
    db_argument_string = ''
    for db in dbs:
        dbname = get_refdb_name(db)
        db_argument_string += ' -d ' + dbname
        db_argument_string += ' -b ' + input_dir + sample_name +'.' + dbname + '.blastout.parsed.txt'
     
    cmd="%s %s --input-annotated-gff %s  --input-kegg-maps  %s  --input-cog-maps %s --output-dir %s"\
           %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_REPORT_FILES']), db_argument_string,
            input_annotated_gff, config_settings['REFDBS'] +'/KO_classification.txt',\
            config_settings['REFDBS'] + '/COG_categories.txt',  output_dir)
    return cmd


def create_genbank_2_sequin_cmd(input_gbk_file, output_sequin_file, config_settings):
    cmd="%s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['GENBANK_2_SEQUIN']),input_gbk_file, output_sequin_file); 
    return cmd


def create_reduce_genebank_cmd(input_gbk_file, output_reduced_gbk_file, config_settings):
    cmd="%s -i %s -o %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['REDUCE_ANNOTATED_FILE']),input_gbk_file, output_reduced_gbk_file ); 
    return cmd

#creates the fasta and the pf files to be processed by pathologic and  a subsequent step will 
# create the PGDB
def create_fasta_and_pf_files_cmd(sample_name, input_reduced_annotated_gbk_file_pat, output_fasta_pf_dir, config_settings):
    cmd="%s -s %s -p %s -a %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['PATHOLOGIC_INPUT']), sample_name, output_fasta_pf_dir, input_reduced_annotated_gbk_file_pat)
    return cmd

def create_create_genetic_elements_cmd(output_fasta_pf_dir, config_settings):
    cmd="%s -p %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_GENETIC_ELEMENTS']),output_fasta_pf_dir)
    return cmd

# this function creates an organism file for pathway tools
def create_create_organism_params_cmd(yaml_file, output_fasta_pf_dir, pgdb_ID, config_settings):
    cmd="%s --yaml-file %s --ptools-dir %s %s" %((config_settings['METAPATHWAYS_PATH'] + config_settings['CREATE_ORGANISM_PARAM']), yaml_file, output_fasta_pf_dir, pgdb_ID)
    return cmd

# this is where the command for PGDB creation is called, by using 
def create_pgdb_using_pathway_tools_cmd(output_fasta_pf_dir, taxonomic_pruning_flag, config_settings):
    cmd="%s -patho %s/" %(config_settings['PATHOLOGIC_EXECUTABLE'],output_fasta_pf_dir)
    if taxonomic_pruning_flag=='no':
        cmd= cmd + " -no-taxonomic-pruning "
        #cmd= cmd + " -no-taxonomic-pruning -no-cel-overview -no-web-cel-ov"
    cmd= cmd + " -no-web-cel-overview"
    return cmd
 

def create_scan_rRNA_seqs_cmd(input_fasta, output_blast, refdb, config_settings, command_Status):

    if command_Status:
       check_an_format_refdb(refdb, 'F',  config_settings)
    else:
       sys.exit(0)

#    if not check_if_db_exists(config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['REFDBS'] +"/" + refdb) and run_command=="yes":
#       print "Error : Reference Database " + refdb + " does not exist!"
#       sys.exit(0)

    cmd="%s -outfmt 6 -num_threads 16  -query %s -out %s -db %s -max_target_seqs 5"%((config_settings['METAPATHWAYS_PATH'] + config_settings['BLASTN_EXECUTABLE']),input_fasta,output_blast, config_settings['METAPATHWAYS_PATH'] + '/' + config_settings['REFDBS'] + '/' + refdb)
    return cmd


# create the command to do rRNA scan statististics
def create_rRNA_scan_statistics(blastoutput, refdbname, bscore_cutoff, eval_cutoff, identity_cutoff, config_settings, rRNA_stat_results):
    cmd= "%s -o %s -b %s -e %s -s %s"  %((config_settings['METAPATHWAYS_PATH'] + config_settings['STATS_rRNA']), rRNA_stat_results, bscore_cutoff, eval_cutoff, identity_cutoff)


    cmd =  cmd +  " -i "  + blastoutput + " -d " + config_settings['METAPATHWAYS_PATH'] +\
              config_settings['REFDBS'] +  refdbname 

    return cmd

# create the command to do tRNA scan statististics
def create_tRNA_scan_statistics(in_file,stat_file, fasta_file,  config_settings):
    cmd= "%s -o %s -F %s  -i %s -T %s  -D %s"  %((config_settings['METAPATHWAYS_PATH'] + config_settings['SCAN_tRNA']) ,\
           stat_file, fasta_file, in_file,\
          (config_settings['METAPATHWAYS_PATH'] + config_settings['EXECUTABLES_DIR'])+'/TPCsignal',\
          (config_settings['METAPATHWAYS_PATH'] + config_settings['EXECUTABLES_DIR'])+'/Dsignal')
    return cmd

# create the command to make the MLTreeMap Images
def create_MLTreeMap_Imagemaker(mltreemap_image_output, mltreemap_final_outputs, config_settings):
    executable_path = config_settings['MLTREEMAP_IMAGEMAKER']
    if not path.isfile( executable_path):
        executable_path = config_settings['METAPATHWAYS_PATH'] + executable_path
    cmd= "%s -i %s -o %s -m a"  %(executable_path, mltreemap_final_outputs, mltreemap_image_output) 
    return cmd

# create the command to do mltreemap calculations
def create_MLTreeMap_Calculations(mltreemap_input_file, mltreemap_output_dir, config_settings):
    executable_path = config_settings['MLTREEMAP_CALCULATION']
    if not path.isfile( executable_path):
        executable_path = config_settings['METAPATHWAYS_PATH'] + executable_path
    cmd= "%s -i %s -o %s"  %(executable_path, mltreemap_input_file, mltreemap_output_dir) 
    return cmd

#gets the parameter value from a category as specified in the 
# parameter file
def get_parameter(config_params,category, field, default = None):
    if category in config_params:
        if field in config_params[category]:
            return config_params[category][field]
        else:
            return default
    return default


# parameter file
def get_make_parameter(config_params,category, field, default = False):
    if category in config_params:
        if field in config_params[category]:
            return config_params[category][field]
        else:
            return default
    return default

def get_pipeline_steps(steps_log_file):
    #print steps_log_file
    try:
       logfile = open(steps_log_file, 'r')
    except IOError:
       print "Did not find " + logfile + "!" 
       print "Try running in \'complete\' run-type"
    else:
       lines = logfile.readlines()
#       for line in lines:
#           print line
    pipeline_steps = None
    return pipeline_steps


def write_run_parameters_file(fileName, parameters):
    try:
       paramFile = open(fileName, 'w')
    except IOError:
       print "Cannot write run parameters to file " + fileName + " !"

#       16s_rRNA      {'min_identity': '40', 'max_evalue': '0.000001', 'min_bitscore': '06', 'refdbs': 'silva_104_rep_set,greengenes_db_DW'}
    paramFile.write("\nRun Date : " + str(date.today()) + " \n")

    paramFile.write("\n\nNucleotide Quality Control parameters\n")
    paramFile.write( "  min length" + "\t" + str(parameters['quality_control']['min_length']) + "\n")

    paramFile.write("\n\nORF prediction parameters\n")
    paramFile.write( "  min length" + "\t" + str(parameters['orf_prediction']['min_length']) + "\n")
    paramFile.write( "  algorithm" + "\t" + str(parameters['orf_prediction']['algorithm']) + "\n")


    paramFile.write("\n\nAmino acid quality control and annotation parameters\n")
    paramFile.write( "  min bit score" + "\t" + str(parameters['annotation']['min_score']) + "\n")
    paramFile.write( "  min seq length" + "\t" + str(parameters['annotation']['min_length']) + "\n")
    paramFile.write( "  annotation reference dbs" + "\t" + str(parameters['annotation']['dbs']) + "\n")
    paramFile.write( "  min BSR" + "\t" + str(parameters['annotation']['min_bsr']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['annotation']['max_evalue']) + "\n")

    paramFile.write("\n\nPathway Tools parameters\n")
    paramFile.write( "  taxonomic pruning " + "\t" + str(parameters['ptools_settings']['taxonomic_pruning']) + "\n")

    paramFile.write("\n\n16s_rRNA search/match parameters\n")
    paramFile.write( "  min identity" + "\t" + str(parameters['16s_rRNA']['min_identity']) + "\n")
    paramFile.write( "  max evalue" + "\t" + str(parameters['16s_rRNA']['max_evalue']) + "\n")
    paramFile.write( "  16s reference dbs" + "\t" + str(parameters['16s_rRNA']['refdbs']) + "\n")

    paramFile.close()


# checks if the necessary files, directories  and executables really exists or not
def check_config_settings(config_settings, file):
   essentialItems= ['METAPATHWAYS_PATH', 'EXECUTABLES_DIR']
   missingItems = []
   for key, value in  config_settings.items():
      # make sure  MetaPathways directory is present
      if key in ['METAPATHWAYS_PATH' ]:
         if not path.isdir(config_settings[key]) :
            print "ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"" %(key, file)  
            print "ERROR: Currently it is set to \"%s\"\n" %( config_settings[key] )  
            missingItems.append(key) 
         continue


      # make sure  REFDB, EXECUTABLES_DIR directories are present
      if key in [ 'REFDBS', 'EXECUTABLES_DIR']:
         if not path.isdir( config_settings['METAPATHWAYS_PATH'] + "/" + config_settings[key]) :
            print "ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"" %(key, file)  
            print "ERROR: Currently it is set to \"%s\"\n" %( config_settings[key] )  
            missingItems.append(key) 
         continue

      # make sure  MetaPaths directory is present
      if key in ['PERL_EXECUTABLE',  'PYTHON_EXECUTABLE' , 'PATHOLOGIC_EXECUTABLE' ]:
         if not path.isfile( config_settings[key]) :
            print "ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"" %(key, file)  
            print "ERROR: Currently it is set to \"%s\"\n" %( config_settings[key] )  
            missingItems.append(key) 
         continue

      
      # check if the desired file exists, if not, then print a message
      if not path.isfile( config_settings['METAPATHWAYS_PATH'] +  value ) :
          if not path.isfile( value ) :
              print "ERROR: Path for \"%s\" is NOT set properly in configuration file \"%s\"" %(key, file)  
              print "ERROR: Currently it is set to \"%s\"\n" %( config_settings['METAPATHWAYS_PATH'] + value ) 
              missingItems.append(key) 
          continue
     
   stop_execution = False
   for item in missingItems:
      if item in essentialItems:
         print """ ERROR: Essential field in setting """ +  item + " is missing in configuration file!"
         stop_execution = True

   if stop_execution ==True:
      print """  """
      print  """ Terminating execution due to missing essentiail  fields in configuration file!"""
      sys.exit(0)

   

# This function reads the pipeline configuration file and sets the 
# paths to differenc scripts and executables the pipeline call
def read_pipeline_configuration( file ):
    try:
       configfile = open(file, 'r')
    except IOError:
       print "Did not find pipeline config " + file + "!" 
    else:
       lines = configfile.readlines()

    config_settings = {}
    for line in lines:
        if not re.match("#",line) and len(line.strip()) > 0 :
           line = line.strip()
           line = re.sub('\t+', ' ', line)
           line = re.sub('\s+', ' ', line)
           line = re.sub('\'', '', line)
           line = re.sub('\"', '', line)
           fields = re.split('\s', line)
           if not len(fields) == 2:
              print """     The following line in your config settings files is set set up yet"""
              print """     Please rerun the pipeline after setting up this line"""
              print """ Error ine line :""" + line
              sys.exit(-1);
              
#           print fields[0] + "\t" + fields[1]
           config_settings[fields[0]] = fields[1]   
    config_settings['METAPATHWAYS_PATH'] = config_settings['METAPATHWAYS_PATH'] + "/"
    config_settings['REFDBS'] = config_settings['REFDBS'] + "/"
    
    check_config_settings(config_settings, file);

    return config_settings

# decide if a command should be run if it is overlay,
# when the expected outputs are present
def shouldRunStep1(run_type, dir , expected_outputs):
    if  run_type =='overlay'  and  doFilesExist(expected_outputs, dir =  dir):
        return False
    else:
        return True


# decide if a command should be run if it is overlay,
# when results are alread computed decide not to run
def shouldRunStep(run_type, expected_output):
    if  run_type =='overlay'  and  path.exists(expected_output):
        return False
    else:
        return True
    #print enable_flag

# has the results to use
def hasResults(expected_output):
    if  path.exists(expected_output):
        return True
    else:
        return False


# has the results to use
def hasResults1(dir , expected_outputs):
    if  doFilesExist(expected_outputs, dir =  dir):
        return True
    else:
        return False


# if the directory is empty then there is not precomputed results
# and so you should decide to run the command
def shouldRunStepOnDirectory(run_type, dirName):
    dirName = dirName + '/*'
    files = glob(dirName)
    if len(files)==0:
      return True
    else:
      return False

# if the command is "redo" then delete all the files
# in the folder and then delete the folder too
def removeDirOnRedo(command_Status, origFolderName):
    if command_Status=='redo' and path.exists(origFolderName) :
       folderName = origFolderName + '/*'
       files = glob(folderName)
       for  f in files:
         remove(f)
       if path.exists(origFolderName): 
         shutil.rmtree(origFolderName)

# if the command is "redo" then delete the file
def removeFileOnRedo(command_Status, fileName):
    if command_Status=='redo' and path.exists(fileName) :
        remove(fileName)
        return True
    else:
        return False


# remove all the files in the directory
def cleanDirectory( folderName):
    folderName = folderName + '/*'
    files = glob(folderName)
    for  f in files:
       remove(f)

# if folder does not exist then create one
def checkOrCreateFolder( folderName ):
    if not path.exists(folderName) :
        makedirs(folderName)
        return False
    else:
        return True

#does the file Exist?
def doFilesExist( fileNames, dir="" ):
    for fileName in fileNames:
       file = fileName
       if dir!='':
         file = dir + "/" + fileName
       if not path.exists(file):
          return False
    return True


# check if all of the metapaths_steps have 
# settings from the valid list [ yes, skip stop, redo]
def checkMetapaths_Steps(config_params):
     choices = {}
     choices['INPUT_TYPE']   = ['fasta', 'gbk_unannotated', 'gbk_annotated', 'gff_unannotated', 'gff_annotated']
     choices['PREPROCESS_FASTA']   = ['yes', 'skip', 'stop', 'redo']
     choices['ORF_PREDICTION']  = ['yes', 'skip', 'stop', 'redo']
     choices['GFF_TO_AMINO']    = ['yes', 'skip', 'stop', 'redo']
     choices['FILTERED_FASTA']  = ['yes', 'skip', 'stop', 'redo']
     choices['COMPUTE_REFSCORE']    = ['yes', 'skip', 'stop', 'redo']
     choices['BLAST_REFDB'] = ['yes', 'skip', 'stop', 'redo', 'grid']
     choices['PARSE_BLAST'] = ['yes', 'skip', 'stop', 'redo']
     choices['SCAN_rRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['STATS_rRNA']  = ['yes', 'skip', 'stop', 'redo']
     choices['ANNOTATE']    = ['yes', 'skip', 'stop', 'redo']
     choices['PATHOLOGIC_INPUT']    = ['yes', 'skip', 'stop', 'redo']
     choices['GENBANK_FILE']    = ['yes', 'skip', 'stop', 'redo']
     choices['CREATE_SEQUIN_FILE']  = ['yes', 'skip', 'stop', 'redo']
     choices['CREATE_REPORT_FILES']  = ['yes', 'skip', 'stop', 'redo']
     choices['SCAN_tRNA']   = ['yes', 'skip', 'stop', 'redo']
     choices['MLTREEMAP_CALCULATION']   = ['yes', 'skip', 'stop', 'redo']
     choices['MLTREEMAP_IMAGEMAKER']    = ['yes', 'skip', 'stop', 'redo']
     choices['PATHOLOGIC']  = ['yes', 'skip', 'stop', 'redo']


     if config_params['metapaths_steps']:
        for step, status  in config_params['metapaths_steps'].iteritems():
            if step in choices:
              if not status in choices[step]:
                print "ERROR: Incorrect setting in your parameter file"
                print "       for step " + step + " as  " + status
                sys.exit(0)
            else:
               print "ERROR: Incorrect step " + step
               sys.exit(0)

def  copy_fna_faa_gff_orf_prediction( source_files, target_files, config_settings) :

     for source, target in zip(source_files, target_files):  
         #print source + ' ' + target
         sourcefile = open(source, 'r')
         targetfile = open(target, 'w')
         sourcelines = sourcefile.readlines()
         for line in sourcelines:
            fprintf(targetfile, "%s\n", line.strip())

         sourcefile.close()
         targetfile.close()


def run_metapathways(input_fp, output_dir, command_handler, command_line_params, config_params, metapaths_config, status_update_callback, config_file, run_type):
#    print input_fp
#    print output_dir
#    print command_handler
#    print params
#    print metapaths_config
#    print status_update_callback

    config_settings = read_pipeline_configuration( config_file )


    #yaml_file = re.sub('.(fasta)|(fna)|(fas)','',input_fp,1, re.M|re.I)

    # SCHEDULE THE PIPELINE STEPS TO RUN
    logger = WorkflowLogger(generate_log_fp(output_dir))
    stepslogger = WorkflowLogger(generate_steps_log_fp(output_dir))
    #pipeline_steps = get_pipeline_steps(stepslogger.get_log_filename());
    
    ####################  IMPORTANT VARIABLES ########################
    commands = []
    
    checkMetapaths_Steps(config_params)

    preprocessed_dir = output_dir + "/preprocessed/"
    checkOrCreateFolder( preprocessed_dir) 

    orf_prediction_dir =  output_dir + "/orf_prediction/"
    checkOrCreateFolder( orf_prediction_dir) 

    genbank_dir =  output_dir + "/genbank/"
    checkOrCreateFolder( genbank_dir)

    output_run_statistics_dir = output_dir + "/run_statistics/" 
    checkOrCreateFolder( output_run_statistics_dir) 

    blast_results_dir =  output_dir + "/blast_results/"
    checkOrCreateFolder( blast_results_dir) 

    output_mltreemap_calculations_dir = output_dir + "/mltreemap_calculations/" 
    checkOrCreateFolder( output_mltreemap_calculations_dir) 

    output_results = output_dir + "/results/" 
    checkOrCreateFolder( output_results) 

    output_results_annotation_table_dir  = output_results + "/annotation_table/" 
    checkOrCreateFolder( output_results_annotation_table_dir) 

    output_results_megan_dir  = output_results + "/megan/" 
    checkOrCreateFolder( output_results_megan_dir) 
   
    output_results_mltreemap_dir  = output_results + "/mltreemap/" 
    checkOrCreateFolder( output_results_mltreemap_dir) 

    mltreemap_image_output = output_results_mltreemap_dir  +"/tables_and_figures/" 
    checkOrCreateFolder(mltreemap_image_output) 

    output_fasta_pf_dir=  output_dir + "/ptools/"
    checkOrCreateFolder(output_fasta_pf_dir) 

    output_results_pgdb_dir  = output_results + "/pgdb/" 
    checkOrCreateFolder( output_results_pgdb_dir) 

    output_results_rRNA_dir  = output_results + "/rRNA/" 
    checkOrCreateFolder( output_results_rRNA_dir) 

    output_results_tRNA_dir  = output_results + "/tRNA/" 
    checkOrCreateFolder( output_results_tRNA_dir) 

    #sys.exit(0)
    # Here the sample name variable has the sample name

    
    sample_name = re.sub(r'[.][a-zA-Z]*$','',input_fp)
    sample_name = path.basename(sample_name)
    sample_name = re.sub('[.]','_',sample_name)

    # IF the format is gff then copy the gff, fna and faa to the orf_prediction folder
    
    if config_params['INPUT']['format'] in ['gff-annotated', 'gff-unannotated' ]:
       input_gff = re.sub('[.][a-zA-Z]*$','',input_fp) + ".gff"
       input_fasta = re.sub('[.][a-zA-Z]*$','',input_fp) + ".fasta"
       input_faa = re.sub('[.][a-zA-Z]*$','',input_fp) + ".faa"

       output_fasta = preprocessed_dir + sample_name + ".fasta"
       output_faa = orf_prediction_dir + sample_name + ".faa"
       if config_params['INPUT']['format'] in ['gff-unannotated' ]:
             output_gff = orf_prediction_dir + sample_name + ".unannot.gff"
       if config_params['INPUT']['format'] in ['gff-annotated' ]:
             output_gff = genbank_dir + sample_name + ".annot.gff"

       command_Status =  get_parameter( config_params,'metapaths_steps','PREPROCESS_FASTA')
       removeFileOnRedo(command_Status, output_fasta)
       removeFileOnRedo(command_Status, output_faa)
       removeFileOnRedo(command_Status, output_gff)

       enable_flag = shouldRunStep(run_type, output_fna) or shouldRunStep(run_type, output_faa)  or shouldRunStep(run_type, output_gff)  

       input_gff_to_fasta_faa_gff_cmd=convert_gff_to_fna_faa_gff([input_gff, input_fasta, input_faa], [ output_gff, output_fasta, output_faa], config_settings) 
       commands.append(["\n1. Converting gff file to fna, faa and gff files ....", input_gff_to_fasta_faa_gff_cmd,\
                       'PREPROCESS_FASTA', command_Status, enable_flag])



    #################################
    
    # IF the format is gbk then convert them  to gff faa and fna files
    if config_params['INPUT']['format'] in ['gbk-annotated', 'gbk-unannotated' ]:
       input_gbk = re.sub('[.][a-zA-Z]*$','',input_fp) + ".gbk"
       output_fna = preprocessed_dir + sample_name + ".fasta"
       output_faa = orf_prediction_dir + sample_name + ".faa"

       if config_params['INPUT']['format'] in ['gbk-unannotated' ]:
             output_gff = orf_prediction_dir + sample_name + ".unannot.gff"
       if config_params['INPUT']['format'] in ['gbk-annotated' ]:
             output_gff = genbank_dir + sample_name + ".annot.gff"

       genbank_to_fna_faa_gff_cmd = convert_gbk_to_fna_faa_gff(input_gbk, output_fna, output_faa, output_gff, config_settings) 
    
       command_Status =  get_parameter( config_params,'metapaths_steps','PREPROCESS_FASTA')
       removeFileOnRedo(command_Status, output_fna)
       removeFileOnRedo(command_Status, output_faa)
       removeFileOnRedo(command_Status, output_gff)
       enable_flag = shouldRunStep(run_type, output_fna) or shouldRunStep(run_type, output_faa)  or shouldRunStep(run_type, output_gff)  

       commands.append(["\n1. Converting genbank file to fna, faa and gff files ....", genbank_to_fna_faa_gff_cmd,\
                       'PREPROCESS_FASTA', command_Status, enable_flag])

    #################################
    if config_params['INPUT']['format'] in ['fasta']:
         output_fas = preprocessed_dir + "/" + sample_name + ".fasta" 
         write_run_parameters_file( output_run_statistics_dir + "/run_parameters.txt", config_params)
         nuc_stats_file = output_run_statistics_dir + "/" + sample_name + ".nuc.stats" 
     
         min_length = float(get_parameter(config_params, 'quality_control', 'min_length', default=300))
         
         input_file = input_fp 
         mapping_txt =  preprocessed_dir + "/" + sample_name + ".mapping.txt" 
         quality_check_cmd = create_quality_check_cmd(min_length, nuc_stats_file, input_file, output_fas, mapping_txt, config_settings) 
     
     
         command_Status =  get_parameter( config_params,'metapaths_steps','PREPROCESS_FASTA')
         removeFileOnRedo(command_Status, output_fas)
         removeFileOnRedo(command_Status, nuc_stats_file)
         # Niels - changed from shouldRunStep1() to shouldRunStep()
         enable_flag = shouldRunStep(run_type, output_fas) or shouldRunStep(run_type, nuc_stats_file )
         
         commands.append( ["\n1. Running Quality Check ....", quality_check_cmd, 'PREPROCESS_FASTA', command_Status, enable_flag])
     
         #################################
     
         # ORF PREDICTION STEP
         #orf_prediction_dir =  output_dir + "/orf_prediction/"
         input_file = output_fas
         output_gff = orf_prediction_dir + sample_name + ".gff"
         orf_prediction_cmd = create_orf_prediction_cmd(" -m -p meta -f gff", "log_file",  input_file, output_gff,config_settings) 
     
          # make folder for orf prediction procesing
         command_Status=  get_parameter( config_params,'metapaths_steps','ORF_PREDICTION')
         removeFileOnRedo(command_Status, output_gff)
         enable_flag = shouldRunStep(run_type,  output_gff) 
     
         #add command
         commands.append(["\n2. Running Orf Prediction ....", orf_prediction_cmd, 'ORF_PREDICTION',command_Status, enable_flag])
         #################################


         # GFF to Fasta
         input_gff = output_gff
         input_fasta= output_fas
         output_faa = orf_prediction_dir + '/' +  sample_name + ".faa"   
         output_fna = orf_prediction_dir + '/' +  sample_name + ".fna"   
         output_gff = orf_prediction_dir + '/' +  sample_name + ".unannot.gff"   
         command_Status= get_parameter( config_params,'metapaths_steps','GFF_TO_AMINO')
         removeFileOnRedo(command_Status, output_faa)
         enable_flag=shouldRunStep1(run_type, orf_prediction_dir,  [ sample_name + '.unannot.gff',  sample_name + '.faa' , sample_name + '.fna'] )  
         create_aa_sequences_cmd = create_aa_orf_sequences_cmd(input_gff, input_fasta, output_fna, output_faa, output_gff, config_settings) 
         commands.append(["\n3. Creating the Amino Acid sequences ....", create_aa_sequences_cmd, 'GFF_TO_AMINO',command_Status, enable_flag])
    #################################

    
    #------ COMPUTE THE STATISTICS OF THE AMINO ACID SEQUENCES
    output_faa = orf_prediction_dir + '/' +  sample_name + ".faa"   

    input_gbk_faa = output_faa
    output_filtered_faa = orf_prediction_dir + '/' +  sample_name + ".qced.faa"   
    amino_stats_file = output_run_statistics_dir + '/' + sample_name + ".amino.stats"

    min_length = float(get_parameter(config_params, 'orf_prediction', 'min_length', default=100))

    create_filtered_amino_acid_sequences_cmd = create_create_filtered_amino_acid_sequences_cmd(input_gbk_faa, 
                                               output_filtered_faa, amino_stats_file, min_length, config_settings)
    command_Status=  get_parameter( config_params,'metapaths_steps','FILTERED_FASTA')
    removeFileOnRedo(command_Status, output_filtered_faa)
    removeFileOnRedo(command_Status, amino_stats_file)
    enable_flag=shouldRunStep(run_type, output_filtered_faa)  or shouldRunStep(run_type, amino_stats_file)  

    #add command
    commands.append( ["\n4. Creating the filtered fasta sequences ....", create_filtered_amino_acid_sequences_cmd, 'FILTERED_FASTA', command_Status, enable_flag])
    #################################


    # COMPUTE THE REFSCORES FROM THE FASTA FILE OF AMINO ACID SEQUENCES
    #blast_results_dir =  output_dir + "/blast_results"
    input_faa =  output_filtered_faa 
    output_refscores =  blast_results_dir + "/" + sample_name + ".refscores" 
    if config_params['INPUT']['format'] in ['fasta', 'gbk-unannotated', 'gff-unannotated' ]:
       refscores_compute_cmd = create_refscores_compute_cmd(input_faa, output_refscores, config_settings) 
       command_Status=  get_parameter( config_params,'metapaths_steps','COMPUTE_REFSCORE')
       removeFileOnRedo(command_Status, output_refscores)
       enable_flag=shouldRunStep(run_type, output_refscores)  

       # add command
       commands.append( ["\n5. Computing refscores for the ORFs ....", refscores_compute_cmd,'COMPUTE_REFSCORE', command_Status, enable_flag])
    #################################

    
    # BLAST THE ORFs AGAINST THE REFERENCE DATABASES  FOR FUNCTIONAL ANNOTATION
    dbstring = get_parameter(config_params, 'annotation', 'dbs', default=None)
    dbs= dbstring.split(",")
    input_faa =  output_filtered_faa 

    if config_params['INPUT']['format'] in ['fasta', 'gbk-unannotated', 'gff-unannotated' ]:
       command_Status=  get_parameter( config_params,'metapaths_steps','BLAST_REFDB')
       message = "\n6. Blasting  ORFs against reference database - "
       for db in dbs:
             dbname = get_refdb_name(db);
   
             blastoutput = blast_results_dir +"/" + sample_name + "." + dbname +  ".blastout"
             removeFileOnRedo(command_Status, blastoutput)
             enable_flag=shouldRunStep(run_type, blastoutput)  
              
             blast_against_refdb_cmd=create_blastp_against_refdb_cmd(input_faa, blastoutput,\
                                       output_dir,  sample_name, db,\
                                      config_settings, config_params, enable_flag)
   
             commands.append([message + db + " ....", blast_against_refdb_cmd, 'BLAST_REFDB_' + dbname, command_Status, enable_flag])
             message = "\n                                               "
    #################################

    # PARSE THE BLAST RESULTS 
    #input_gbk_faa =  output_gbk_faa 
    message = "\n7. Parsing blast outputs for reference database - "
    if config_params['INPUT']['format'] in ['fasta', 'gbk-unannotated', 'gff-unannotated' ]:
       for db in dbs:
          dbname = get_refdb_name(db);
   
          input_db_blastout = blast_results_dir + "/" + sample_name + "." + dbname + ".blastout"
          refscorefile = blast_results_dir + "/" + sample_name + ".refscores"
          output_db_blast_parsed = blast_results_dir +"/" + sample_name + "." + dbname + ".blastout.parsed.txt"

          command_Status=  get_parameter( config_params,'metapaths_steps','PARSE_BLAST')
          removeFileOnRedo(command_Status, output_db_blast_parsed)
          enable_flag=shouldRunStep(run_type, output_db_blast_parsed)  
   
          dbmapfile = config_settings['METAPATHWAYS_PATH'] + config_settings['REFDBS'] + "/" +  db + "-names.txt"

          make_sure_map_file_exists(dbmapfile)
       
   
          parse_blast_cmd = create_parse_blast_cmd(input_db_blastout, refscorefile,  dbname, dbmapfile,  config_settings, config_params)
   
          #add command
       
          commands.append( [message  + dbname + " ....",\
                 parse_blast_cmd, 'PARSE_BLAST_'+dbname, get_parameter( config_params,'metapaths_steps','PARSE_BLAST'), enable_flag])
          message = "\n                                                  "
    #################################


    
    # BLAST AGAINST rRNA REFERENCE DATABASES
    input_fasta = preprocessed_dir + "/" + sample_name + ".fasta" 

    refdbstring = get_parameter(config_params, '16s_rRNA', 'refdbs', default=None)
    refdbnames=refdbstring.split(',')


    message = "\n8. Scan for rRNA sequences in reference database - "
    for refdbname in refdbnames:
       rRNA_blastout= blast_results_dir + "/" + sample_name + ".rRNA." + get_refdb_name(refdbname) + ".blastout"

       command_Status=  get_parameter( config_params,'metapaths_steps','SCAN_rRNA')
       scan_rRNA_seqs_cmd = create_scan_rRNA_seqs_cmd(input_fasta, rRNA_blastout, refdbname, config_settings, command_Status)

       
       removeFileOnRedo(command_Status, rRNA_blastout)
       enable_flag=shouldRunStep(run_type, rRNA_blastout)  

       commands.append( [message + refdbname + "...." , scan_rRNA_seqs_cmd, 'SCAN_rRNA_' + refdbname, command_Status, enable_flag])
       message = "\n                                                   "
    #################################

    
    # COMPUTE rRNA STATISTICS 
    bscore_cutoff = get_parameter(config_params,'16s_rRNA', 'min_bitscore', default=27)
    eval_cutoff = get_parameter(config_params, '16s_rRNA', 'max_evalue', default=6)
    identity_cutoff = get_parameter(config_params, '16s_rRNA', 'min_identity', default=40)

    message = "\n9. Gathering rRNA stats .... "
    for refdbname in refdbnames:
        rRNA_stat_results= output_results_rRNA_dir + sample_name + "." + get_refdb_name(refdbname) + ".rRNA.stats.txt"
        rRNA_blastout= blast_results_dir + "/" + sample_name + ".rRNA." + get_refdb_name(refdbname) + ".blastout"

        rRNA_stats_cmd = create_rRNA_scan_statistics(rRNA_blastout, refdbname, bscore_cutoff, eval_cutoff, identity_cutoff, config_settings, rRNA_stat_results)

        command_Status=  get_parameter( config_params,'metapaths_steps','STATS_rRNA')
        removeFileOnRedo(command_Status, rRNA_stat_results)
        enable_flag=shouldRunStep(run_type, rRNA_stat_results)  

        #add command
        commands.append( [message + get_refdb_name(refdbname), rRNA_stats_cmd, 'STATS_' + refdbname, command_Status, enable_flag])
        message = "\n                             "
    #################################
    
    #SCAN FOR tRNA genes
    tRNA_stats_output = output_results_tRNA_dir + "/" + sample_name +  ".tRNA.stats.txt"
    tRNA_fasta_output = output_results_tRNA_dir + "/" + sample_name +  ".tRNA.fasta"
    input_fasta = preprocessed_dir + "/" + sample_name + ".fasta" 

    scan_tRNA_seqs_cmd = create_tRNA_scan_statistics(input_fasta, tRNA_stats_output, tRNA_fasta_output,  config_settings)

    command_Status=  get_parameter( config_params,'metapaths_steps','SCAN_tRNA')
    removeFileOnRedo(command_Status, tRNA_stats_output)
    removeFileOnRedo(command_Status, tRNA_fasta_output)
    enable_flag=shouldRunStep1(run_type, output_results_tRNA_dir, [sample_name + ".tRNA.stats.txt", sample_name + ".tRNA.fasta"] )

    #add command
    commands.append( ["\n10. Scanning for tRNA using tRNAscan 1.4  ", scan_tRNA_seqs_cmd, 'SCAN_tRNA', command_Status, enable_flag])
    #################################

    # ANNOTATION STEP 
    input_unannotated_gff = orf_prediction_dir +"/" + sample_name+".unannot.gff"
    output_annotated_gff  = genbank_dir +"/" + sample_name+".annot.gff"
    output_comparative_annotation  =  output_results_annotation_table_dir + '/' + sample_name

    if config_params['INPUT']['format'] in ['fasta', 'gbk-unannotated', 'gff-unannotated' ]:
        parse_blasts = [] 
        options =''
        for db in dbs:
            parse_blasts.append( get_refdb_name(db) );
         
    
        rRNArefdbstring = get_parameter(config_params, '16s_rRNA', 'refdbs', default=None)
        rRNArefdbs=rRNArefdbstring.split(',')
        for rRNArefdb in rRNArefdbs:
            rRNA_stat_results= output_results_rRNA_dir + sample_name + '.' + get_refdb_name(rRNArefdb) + '.rRNA.stats.txt' 
            isrRNA_stat_available=  hasResults(rRNA_stat_results)  
            if  isrRNA_stat_available:
                options += " --rRNA_16S " +  rRNA_stat_results 
    
        tRNA_stat_results= output_results_tRNA_dir + '/' + sample_name + '.tRNA.stats.txt' 
        istRNA_stat_available=  hasResults(tRNA_stat_results)  
        if  istRNA_stat_available:
           options += " --tRNA " +  tRNA_stat_results 
    
        # create the command
        annotate_gbk_cmd = create_annotate_genebank_cmd(sample_name, input_unannotated_gff,\
                             output_annotated_gff,  blast_results_dir, options,\
                              parse_blasts, output_comparative_annotation, config_settings)
    
        command_Status=  get_parameter( config_params,'metapaths_steps','ANNOTATE')
        removeFileOnRedo(command_Status, output_annotated_gff)
        enable_flag=shouldRunStep(run_type, output_annotated_gff)  
    
        commands.append( ["\n11. Annotate gff files    ....", annotate_gbk_cmd,'ANNOTATE',command_Status, enable_flag])
    #################################
 
    
    # CREATE A GENBANK FILE, PTOOLS INPUT and SEQUIN FILE
    input_annot_gff =  output_annotated_gff
    input_nucleotide_fasta = preprocessed_dir + "/" + sample_name + ".fasta" 
    input_amino_acid_fasta = output_filtered_faa

    output_annot_gbk= genbank_dir + '/' + sample_name +  '.gbk'
    gbk_command_Status=  get_parameter(config_params,'metapaths_steps','GENBANK_FILE')
    removeFileOnRedo(gbk_command_Status, output_annot_gbk)
    enable_gbk_flag=shouldRunStep(run_type, output_annot_gbk)  

    outputs = {}
    if enable_gbk_flag:
       outputs['gbk'] = output_annot_gbk

    files = [ '0.pf', '0.fasta', 'genetic-elements.dat', 'organism-params.dat']
    ptinput_command_Status=  get_parameter( config_params,'metapaths_steps','PATHOLOGIC_INPUT')
    removeDirOnRedo(ptinput_command_Status, output_fasta_pf_dir)
    enable_ptinput_flag=shouldRunStep1(run_type, output_fasta_pf_dir, files)

    # Niels - Issue #12 removed some strange output appearing to standard out
    
    if enable_ptinput_flag:
        outputs['ptinput'] = output_fasta_pf_dir

    if enable_ptinput_flag or enable_gbk_flag:
       enable_flag = True
    else:
       enable_flag = False

    if ptinput_command_Status=='yes' or gbk_command_Status=='yes':
       command_Status = 'yes'
    else: 
       command_Status = 'skip'
    

    genbank_annotation_table_cmd = create_genbank_ptinput_sequin_cmd(input_annot_gff, input_nucleotide_fasta,  input_amino_acid_fasta, outputs, config_settings)
    #add command
    commands.append( ["\n12. Create genbank/sequin/ptools input command    ....", genbank_annotation_table_cmd,'GENBANK_FILE', command_Status, enable_flag])
    #################################
 
    
    # CREATE A REPORT TABLES WITH FUNCTION AND TAXONOMY
    input_annotated_gff = output_annotated_gff
    output_annot_table = output_results_annotation_table_dir + '/functional_and_taxonomic_table.txt'
    if config_params['INPUT']['format'] in ['fasta', 'gbk-unannotated', 'gff-unannotated' ]:
       command_Status=  get_parameter(config_params,'metapaths_steps','CREATE_REPORT_FILES')
       removeFileOnRedo(command_Status, output_annot_table)
       enable_flag=shouldRunStep(run_type, output_annot_table)  
       message = "\n13. Creating KEGG and COG hits table and LCA based taxonomy table  "
       input_dir = blast_results_dir +"/" #D+ sample_name + "." + dbname + ".blastout.parsed"
       create_report_cmd=create_report_files_cmd(dbs, input_dir, input_annotated_gff, sample_name, output_results_annotation_table_dir, config_settings)
       commands.append( [message  + " ....", create_report_cmd, 'CREATE_REPORT_FILES', command_Status, enable_flag])

    #################################

    # IS REFSEQ BLAST OUTPUT PRESENT?
    # Niels - Issue #12: change this to detect if RefSeq is scheduled or not rather than just 
    # having the results here.
    # refseqblastoutput = blast_results_dir + "/" + sample_name + "." + 'refseq' + ".blastout"
    blast_dbs = (get_parameter( config_params,'annotation','dbs'))
    
    # if not  (path.exists(refseqblastoutput) :
    if not (re.search(".*refseq.*",blast_dbs, re.IGNORECASE)):
        print "WARNING: Refseq annotation is not scheduled!"
        print "         Taxonomic information will not be found in the annotation table."
    #################################

    # CREATE THE SEQUIN FILE
    input_annot_gbk= genbank_dir + '/' + sample_name +  '.gbk'
    output_sequin_file = genbank_dir + '/' + sample_name + '.sequin'
    command_Status=  get_parameter( config_params,'metapaths_steps','CREATE_SEQUIN_FILE');
    genbank_2_sequin_cmd = create_genbank_2_sequin_cmd(input_annot_gbk, output_sequin_file, config_settings)

    removeFileOnRedo(command_Status, output_sequin_file)
    enable_flag=shouldRunStep(run_type, output_sequin_file)  
    commands.append( ["\n11b. Create SEQUIN file    ....", genbank_2_sequin_cmd,'CREATE_SEQUIN_FILE', command_Status, enable_flag])
    #################################


    #MLTreeMap Calculations 
    mltreemap_input_file = input_fp

    command_Status=  get_parameter( config_params,'metapaths_steps','MLTREEMAP_CALCULATION')
    removeDirOnRedo(command_Status, output_mltreemap_calculations_dir+'/various_outputs')
    removeDirOnRedo(command_Status, output_mltreemap_calculations_dir+'/final_outputs')
    removeDirOnRedo(command_Status, output_mltreemap_calculations_dir+'/final_RAxML_outputs')
    enable_flag=shouldRunStepOnDirectory(run_type, output_mltreemap_calculations_dir)

    mltreemap_calculation_cmd = create_MLTreeMap_Calculations(mltreemap_input_file, output_mltreemap_calculations_dir, config_settings)
    commands.append( ["\n16. MLTreeMap Calculations  ", mltreemap_calculation_cmd, 'MLTREEMAP_CALCULATION', command_Status, enable_flag])


    #MLTreeMap Image Creating 
    mltreemap_final_outputs = output_mltreemap_calculations_dir + "/final_outputs/" 

    mltreemap_imagemaker_cmd = create_MLTreeMap_Imagemaker(mltreemap_image_output, mltreemap_final_outputs, config_settings)
    command_Status=  get_parameter( config_params,'metapaths_steps','MLTREEMAP_IMAGEMAKER')
    enable_flag=shouldRunStep(run_type, mltreemap_image_output)

    #add command
    commands.append( ["\n17. Making MLTreeMap Images  ", mltreemap_imagemaker_cmd, 'MLTREEMAP_IMAGEMAKER', command_Status, enable_flag])
    #################################


    # GENERATE PGDB using pathologic
    # note "redo" is a forceful way of recomputing the pgdb 
    command_Status=  get_parameter( config_params,'metapaths_steps','PATHOLOGIC')

    enable_flag = False
    if (run_type in ['overwrite', 'overlay']  and  command_Status== 'yes'):
       enable_flag = True
       #remove the pathologic dir for the pgdb with the same sample name since a new pgdb is requested
       remove_existing_pgdb( sample_name, config_settings['PATHOLOGIC_EXECUTABLE'])


    taxonomic_pruning_flag = get_parameter(config_params,'ptools_settings', 'taxonomic_pruning')
    create_pgdb_cmd = create_pgdb_using_pathway_tools_cmd(output_fasta_pf_dir, taxonomic_pruning_flag, config_settings)
 #   print create_pgdb_cmd 
 #   print taxonomic_pruning_flag
    commands.append( ["\n18. Create PGDB  ....", create_pgdb_cmd, 'PATHOLOGIC', command_Status, enable_flag])
    #################################



    print """  """
    print """  ********************************************************** """
    print """  ********************   Running  metapaths ******************* """
    print """  ********************************************************** """
    print """              """  +  input_fp + """                       """
    print """  ********************************************************** """
    if run_type == 'dry-run':
         dry_run_status(commands)
    else:
         command_handler(commands, status_update_callback, logger, stepslogger, command_line_params)


