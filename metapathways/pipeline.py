"""The main script that calls the pipeline """

try:
    import sys
    import traceback
    import re
    import inspect
    import shutil
    import glob
    import getpass

    from os import makedirs, sys, listdir, environ, path, _exit, system
    import argparse

    from metapathways import errorcodes as errormod
    from metapathways import metapathways_utils  as mputils
    from metapathways import parse as parsemod
    from metapathways import sysutil as sysutils
    from metapathways import metapathways_steps as mpsteps
    from metapathways import parameters as paramsmod
    from metapathways import diagnoze as diagnoze
    from metapathways import sampledata as sampledata
    from metapathways import general_utils as gutils
    from metapathways import _version
except:
   print("""Could not load some user defined  module functions""")
   print(traceback.print_exc(10))
   sys.exit(3)

__author__ = _version.__author__
__version__ = _version.__version__
__copyright__ = _version.__copyright__
__maintainer__ = _version.__maintainer__
__status__ = _version.__status__

cmd_folder = path.abspath(path.split(inspect.getfile( inspect.currentframe() ))[0])


PATHDELIM =  sysutils.pathDelim()

#config = load_config()
metapaths_param = """config/template_param.txt""";

script_info={}
script_info['brief_description'] = """A workflow script for making PGDBs from metagenomic sequences"""
script_info['script_description'] = \
    """ This script starts a MetaPathways pipeline run. It requires an input directory of fasta or genbank files
    containing sequences to process, an output directory for results to be placed. It also requires the
    configuration files, template_config.txt and template_param.txt in the config/ directory, to be updated with the
    location of resources on your system.
    """
script_info['script_usage'] = []



def runParser():
    parser = argparse.ArgumentParser(description='MetaPathways command-line tool for annotation of contigs.')
    subparsers = parser.add_subparsers(dest="command")
    run_parser = subparsers.add_parser(
        'run',
        description='Minimum REQUIRED Command:\n'
                    'MetaPathways run -i INPUT_FILE -o OUTPUT_DIR -d REFDB_DIR\n\n',
        usage='Metapathways run [options]', formatter_class=argparse.RawTextHelpFormatter)
    
    # Minimum required args
    req_args = run_parser.add_argument_group('Minimum Required Arguments')
    req_args.add_argument("-i", "--input_file", required=True,
                        help='path to the input fasta file/input dir [REQUIRED]')
    req_args.add_argument("-o", "--output_dir", required=True,
                        help='path to the output directory [REQUIRED]')
    req_args.add_argument("-d", "--refdb_dir", required=True,
                        help="path to the reference DB [REQUIRED]")
    
    
    ### OPTIONAL ARGUMENTS ###

    # Quality Control
    qc_args = run_parser.add_argument_group('Quality Controls Arguments')
    qc_args.add_argument('--input_format', type=str, default='fasta', choices=['fasta', 'fasta-amino'],
                         help='Input format, FASTA support only [fasta]')
    qc_args.add_argument('--qc_min_length', type=int, default=180,
                         help='Minimum length for quality control [180]')
    qc_args.add_argument('--qc_delete_replicates', type=str, default='yes', choices=['yes', 'no'],
                         help='Delete replicates in quality control [yes]')

    # ORF prediction
    orf_args = run_parser.add_argument_group('ORF Prediction Arguments')
    orf_args.add_argument('--orf_strand', type=str, default='both', choices=['pos', 'neg', 'both'],
                          help='Strand for ORF prediction [both]')
    orf_args.add_argument('--orf_algorithm', type=str, default='prodigal',
                          help='Algorithm for ORF prediction, Prodigal support only [prodigal]')
    orf_args.add_argument('--orf_min_length', type=int, default=60,
                          help='Minimum ORF length [60]')
    orf_args.add_argument('--orf_translation_table', type=int, default=11,
                          help='Translation table for ORF prediction, see Prodigal for translation tables [11] ')
    orf_args.add_argument('--orf_mode', type=str, default='meta', choices=['single', 'meta'],
                          help='Mode for ORF prediction [meta]')

    # Functional annotation
    func_args = run_parser.add_argument_group('Functional Annotation Arguments')
    func_args.add_argument('--annotation_algorithm', type=str, default='FAST', choices=['FAST', 'BLAST'],
                           help='Algorithm for ORF annotation [FAST]')
    func_args.add_argument('--annotation_dbs', nargs='+', default=['swissprot'],
                           help='Database(s) for annotation, space-separated list [swissprot]')
    func_args.add_argument('--annotation_min_bsr', type=float, default=0.4,
                           help='Minimum BSR for annotation [0.4]')
    func_args.add_argument('--annotation_max_evalue', type=float, default=0.000001,
                           help='Maximum e-value for annotation [0.000001]')
    func_args.add_argument('--annotation_min_score', type=int, default=20,
                           help='Minimum score for annotation [20]')
    func_args.add_argument('--annotation_min_length', type=int, default=45,
                           help='Minimum length for annotation [45]')
    func_args.add_argument('--annotation_max_hits', type=int, default=5,
                           help='Maximum hits for annotation [5]')
    func_args.add_argument('--annotation_run_mode', type=str, default='pervol', choices=['default', 'pervol'],
                           help='Run mode for annotation, FAST only [pervol]')

    # rRNA annotation
    rrna_args = run_parser.add_argument_group('rRNA Annotation Arguments')
    rrna_args.add_argument('--rRNA_refdbs', type=str, nargs='+',
                           default=['SILVA_138.1_LSURef_NR99_tax_silva_trunc', 'SILVA_138.1_SSURef_NR99_tax_silva_trunc'],
                           help='Reference databases for rRNA annotation, space-separated list\n' \
                           '[SILVA_138.1_LSURef_NR99_tax_silva_trunc SILVA_138.1_SSURef_NR99_tax_silva_trunc]')
    rrna_args.add_argument('--rRNA_max_evalue', type=float, default=0.000001,
                           help='Maximum e-value for rRNA annotation [0.000001]')
    rrna_args.add_argument('--rRNA_min_identity', type=int, default=20,
                           help='Minimum identity for rRNA annotation [20]')
    rrna_args.add_argument('--rRNA_min_bitscore', type=int, default=50,
                           help='Minimum bitscore for rRNA annotation [50]')

    # Pathway tools
    ptools_args = run_parser.add_argument_group('Pathway Tools Preprocessing Arguments')
    ptools_args.add_argument('--ptools_taxonomic_pruning', type=str, default='no', choices=['yes', 'no'],
                             help='Taxonomic pruning in pathway tools [no]')

    # Read mapping
    reads_args = run_parser.add_argument_group('Read Mapping Arguments')
    reads_args.add_argument("-r", "--readsdir",
                        help="location of the raw fastq data for RPKM and TPM [optional]")

    # Pipeline execution flags
    pipe_args = run_parser.add_argument_group('Pipeline Step Arguments')
    pipe_args.add_argument('--PREPROCESS_INPUT', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: PREPROCESS_INPUT [yes]')
    pipe_args.add_argument('--ORF_PREDICTION', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: ORF_PREDICTION [yes]')
    pipe_args.add_argument('--FILTER_AMINOS', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: FILTER_AMINOS [yes]')
    pipe_args.add_argument('--SCAN_rRNA', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: SCAN_rRNA [yes]')
    pipe_args.add_argument('--SCAN_tRNA', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: SCAN_tRNA [yes]')
    pipe_args.add_argument('--FUNC_SEARCH', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: FUNC_SEARCH [yes]')
    pipe_args.add_argument('--PARSE_FUNC_SEARCH', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: PARSE_FUNC_SEARCH [yes]')
    pipe_args.add_argument('--ANNOTATE_ORFS', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: ANNOTATE_ORFS [yes]')
    pipe_args.add_argument('--GENBANK_FILE', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: GENBANK_FILE [yes]')
    pipe_args.add_argument('--CREATE_ANNOT_REPORTS', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: CREATE_ANNOT_REPORTS [yes]')
    pipe_args.add_argument('--PATHOLOGIC_INPUT', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: PATHOLOGIC_INPUT [yes]')
    pipe_args.add_argument('--COMPUTE_TPM', type=str, default='yes', choices=['yes', 'skip', 'redo'],
                           help='Step: COMPUTE_TPM [yes]')
    pipe_args.add_argument('--force_redo', action="store_true", default=False,
                           help="Redo all steps [False]")

    # Other arguments
    misc_args = run_parser.add_argument_group('Miscellaneous Arguments')
    misc_args.add_argument("-s", "--samples", nargs='+', action="append", default=[],
                        help="process only specific samples, space-separated list")
    misc_args.add_argument("-t", "--threads", default=4, type=int,
                        help="max number of cores to use in multithreaded steps [4]")
    misc_args.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="print more information on the stdout")

    return parser


def dbParser():
    parser = argparse.ArgumentParser(description='MetaPathways command-line tool for building databases.')
    subparsers = parser.add_subparsers(dest="command")
    build_parser = subparsers.add_parser(
        'build_db',
        description='Minimum REQUIRED Command:\n'
                    'Metapathways build_db -o output_dir -b db_type\n\n',
        usage='Metapathways build_db [options]', formatter_class=argparse.RawTextHelpFormatter)

    build_parser.add_argument("-o", "--output_dir", dest="output_dir", required=True,
                        help='path to save reference database [REQUIRED]')
    build_parser.add_argument('-b', '--db_type', dest="db_type", required=True,
                        help="Build version of database [REQUIRED].\nacceptable values are:" \
                            "\n\tstage_fast_full,\n\tstage_fast_lite,\n\tstage_fast_noMeta," \
                            "\n\tstage_blast_full,\n\tstage_blast_lite")
    build_parser.add_argument("-t", "--threads", dest="num_cpus", default='all',
                        help="max number of cores to use in multithreaded steps [DEFAULT all]")

    return parser

def msParser():
    parser = argparse.ArgumentParser(description='MetaPathways command-line tool for MAG splitting.')
    subparsers = parser.add_subparsers(dest="command")
    mag_parser = subparsers.add_parser(
        'mag_split',
        description='Minimum REQUIRED Command:\n'
                    'Metapathways mag_split -o output_dir -m contig_mag_map\n\n',
        usage='Metapathways mag_split [options]', formatter_class=argparse.RawTextHelpFormatter)

    mag_parser.add_argument("-o", "--output_dir", dest="output_dir", required=True,
                        help='path where MP output was saved [REQUIRED]')
    mag_parser.add_argument("-m", "--contig_mag_map", dest="mag_map", required=True,
                        help="TSV file that contains contig-to-MAG mapping [REQUIRED]")

    return parser


def ptParser():
    parser = argparse.ArgumentParser(description='MetaPathways command-line tool for ptools.')
    subparsers = parser.add_subparsers(dest="command")
    ptools_parser = subparsers.add_parser(
        'ptools',
        description='Minimum REQUIRED Command:\n'
                    'Metapathways ptools -o output_dir\n\n',
        usage='Metapathways ptools [options]', formatter_class=argparse.RawTextHelpFormatter)

    ptools_parser.add_argument("-o", "--output_dir", dest="output_dir", required=True,
                        help='path where MP output was saved [REQUIRED]')
    ptools_parser.add_argument("--tag", dest="tag",
                        help="Custom name for ePGDB [optional]")
    ptools_parser.add_argument("--container", action="store_true", dest="container", default=False,
                        help="Flag only used in containerized env [special flag]")

    return parser


def derive_sample_name(filename):
    basename = path.basename(filename)

    shortname = re.sub('[.]gbk$','',basename, re.IGNORECASE)
    shortname = re.sub('[.](fasta|fas|fna|faa|fa|fna.gz)$','',shortname, re.IGNORECASE)
    return shortname

def remove_unspecified_samples(input_output_list, sample_subset,  globalerrorlogger = None):
   """ keep only the samples that are specified  before processing  """
   shortened_names = {}
   input_sample_list = list(input_output_list.keys())

   for sample_name in input_sample_list:
      short_sample_name = derive_sample_name(sample_name)
      if len(short_sample_name) > 35:
         gutils.eprintf("ERROR\tSample name %s must not be longer than 35 characters!\n",short_sample_name)
         if globalerrorlogger:
             globalerrorlogger.printf("ERROR\tSample name %s must not be longer than 35 characters!\n",short_sample_name)
      if sample_subset and  not derive_sample_name(sample_name) in sample_subset:
         del input_output_list[sample_name]



def check_for_error_in_input_file_name(shortname, globalerrorlogger=None):

    """  creates a list of  input output pairs if input is  an input dir """
    clean = True
    if re.search(r'[.]',shortname):
         gutils.eprintf("ERROR\tSample name %s contains a '.' in its name!\n",shortname)
         if globalerrorlogger:
            globalerrorlogger.printf("ERROR\tSample name %s contains a '.' in its name!\n",shortname)
         clean = False

    if clean:
         return clean

    errmessage = """Sample names (input assembly names before extension) must consist only of alphanumeric characters and and underscores"""
    gutils.eprintf("ERROR\t%s\n",errmessage)
    if globalerrorlogger:
        globalerrorlogger.printf("ERROR\t%s\n",errmessage)
    #    mputils.exit_process(errmessage + "Exiting!" + "\n", logger=globalerrorlogger)
    return False


def create_an_input_output_pair(input_file, output_dir,  globalerrorlogger=None):
    """ creates an input output pair if input is just an input file """

    input_output = {}

    if not re.search(r'.(fasta|fas|fna|faa|gbk|gff|fa|fna.gz)$',input_file, re.IGNORECASE):
       return input_output

    shortname = None
    shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE)
    shortname = re.sub('[.](fasta|fas|fna|faa|fa|fna.gz)$','',input_file, re.IGNORECASE)
    shortname = re.sub(r'.*' + PATHDELIM ,'',shortname)

    if  check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger):
       input_output[input_file] = path.abspath(output_dir) + PATHDELIM + shortname

    return input_output


def create_input_output_pairs(input_dir, output_dir,  globalerrorlogger=None):
    """  creates a list of  input output pairs if input is  an input dir """
    fileslist =  listdir(input_dir)

    gbkPatt = re.compile('[.]gbk$',re.IGNORECASE)
    fastaPatt = re.compile('[.](fasta|fas|fna|faa|fa|fna.gz)$',re.IGNORECASE)
    gffPatt = re.compile('[.]gff$',re.IGNORECASE)

    input_files = {}
    for input_file in fileslist:

       shortname = None
       result = None

       result =  gbkPatt.search(input_file)
       if result:
         shortname = re.sub('[.]gbk$','',input_file, re.IGNORECASE)

       if result==None:
          result =  fastaPatt.search(input_file)
          if result:
             shortname = re.sub('[.](fasta|fas|fna|faa|fa|fna.gz)$','',input_file, re.IGNORECASE)

       if shortname == None:
          continue

       if re.search('.(fasta|fas|fna|faa|gff|gbk|fa|fna.gz)$',input_file, re.IGNORECASE):
          if check_for_error_in_input_file_name(shortname, globalerrorlogger=globalerrorlogger):
             input_files[input_file] = shortname

    paired_input = {}
    for key, value in input_files.items():
       paired_input[input_dir + PATHDELIM + key] = path.abspath(output_dir) + PATHDELIM + value

    return paired_input

def removeSuffix(sample_subset_in):
    sample_subset_out = []
    for sample_name in sample_subset_in:
       mod_name = re.sub('.(fasta|fas|fna|faa|gff|gbk|fa|fna.gz)$','',sample_name)
       sample_subset_out.append(mod_name)

    return sample_subset_out

def halt_on_invalid_input(input_output_list, filetypes, sample_subset):
    for samplePath in input_output_list.keys():
       sampleName =  path.basename(input_output_list[samplePath])

       if filetypes[samplePath][0]=='UNKNOWN':
          gutils.eprintf("ERROR\tIncorrect input sample %s. Check for bad characters or format\n!", samplePath)
          return False

       ''' in the selected list'''
       if not sampleName in sample_subset:
          continue

    return True

def report_missing_filenames(input_output_list, sample_subset, logger=None):
    foundFiles = {}
    for samplePath in input_output_list.keys():
       sampleName =  path.basename(input_output_list[samplePath])
       foundFiles[sampleName] =True

    for sample_in_subset in sample_subset:
       if not sample_in_subset in foundFiles:
          gutils.eprintf("ERROR\tCannot find input file for sample %s\n!", sample_in_subset)
          if logger:
             logger.printf("ERROR\tCannot file input for sample %s!\n", sample_in_subset)


def create_arg_dict(parser):
    # Nested dictionary to hold the structure
    arg_structure = {}

    # Go through each action in the parser
    for action in parser._actions:
        # Check if it is an argument group by type
        if isinstance(action, argparse._ArgumentGroup):
            group_name = action.title
            arg_structure[group_name] = {}

            # For each argument in this group
            for arg in action._group_actions:
                # Extract the primary argument name (removing '--' prefix for clarity)
                arg_name = arg.option_strings[0].lstrip('-')
                arg_structure[group_name][arg_name] = None  # Initial placeholder, you can replace it with actual values later

    return arg_structure


def run():
    argv = sys.argv
    gutils.eprintf("%-10s:%s\n" % ('COMMAND', ' '.join(argv)))
    parser = runParser()
    args = parser.parse_args()
    if args.force_redo:
        steps_list = ['PREPROCESS_INPUT', 'ORF_PREDICTION', 'FILTER_AMINOS', 'SCAN_rRNA',
                      'SCAN_tRNA', 'FUNC_SEARCH', 'PARSE_FUNC_SEARCH', 'ANNOTATE_ORFS',
                      'GENBANK_FILE', 'CREATE_ANNOT_REPORTS', 'PATHOLOGIC_INPUT', 'COMPUTE_TPM']
        for step in steps_list:
            setattr(args, step, 'redo')    
    params = parsemod.populate_dict(args)

    # initialize the input directory or file
    input_fp = params['Minimum Required Arguments']['input_file']
    output_dir = path.abspath(params['Minimum Required Arguments']['output_dir'])
    verbose = params['Miscellaneous Arguments']['verbose']

    # Subset inputs if specified
    sample_subset = removeSuffix(params['Miscellaneous Arguments']['samples'])
    run_type = 'safe'

    # Create output directory if it doesn't exist
    if not path.exists(output_dir):
        makedirs(output_dir)

    # Set output type if verbose or not
    if verbose:
        status_update_callback = gutils.print_to_stdout
    else:
        status_update_callback = gutils.no_status_updates

    # Initialize the commandline params dictionary
    command_line_params = {}
    command_line_params['verbose'] = verbose
    command_line_params['ptools-compact-mode'] = 'yes' # this is a legacy setting, to be removed someday

    """ load the sample inputs  it expects either a fasta
        file or  a directory containing fasta and yaml file pairs
    """
    print("output dir", output_dir)
    globalerrorlogger = mputils.WorkflowLogger(mputils.generate_log_fp(output_dir, basefile_name = 'global_errors_warnings'), open_mode='w')

    input_output_list = {}
    if path.isfile(input_fp):
       """ check if it is a file """
       input_output_list = create_an_input_output_pair(input_fp, output_dir,  globalerrorlogger=globalerrorlogger)
    else:
       if path.exists(input_fp):
          """ check if dir exists """
          input_output_list = create_input_output_pairs(input_fp, output_dir, globalerrorlogger=globalerrorlogger)
       else:
          """ must be an error """
          gutils.eprintf("ERROR\tNo valid input sample file or directory containing samples exists .!")
          gutils.eprintf("ERROR\tAs provided as arguments in the -in option.!\n")
          mputils.exit_process("ERROR\tAs provided as arguments in the -in option.!\n")

    """ these are the subset of sample to process if specified
        in case of an empty subset process all the sample """

    # remove all samples that are not specifed unless sample_subset is empty
    remove_unspecified_samples(input_output_list, sample_subset, globalerrorlogger = globalerrorlogger)

    # add check the config parameters
    sorted_input_output_list = sorted(input_output_list.keys())
    filetypes = gutils.check_file_types(sorted_input_output_list)

    #stop on invalid samples
    if not halt_on_invalid_input(input_output_list, filetypes, sample_subset):
       globalerrorlogger.printf("ERROR\tInvalid inputs found. Check for file with bad format or characters!\n")

    # make sure the sample files are found
    report_missing_filenames(input_output_list, sample_subset, logger=globalerrorlogger)
    parameter =  paramsmod.Parameters()

    config = gutils.someclass()
    config.refdb_dir = params['Minimum Required Arguments']['refdb_dir']
    
    if not diagnoze.staticDiagnose(params, config,  logger = globalerrorlogger):
        gutils.eprintf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        globalerrorlogger.printf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        return

    samplesData = {}
    # PART1 before the blast

    configs = {
        "PREPROCESS_INPUT"     : "MetaPathways_filter_input",
        "PREPROCESS_AMINOS"    : "MetaPathways_preprocess_amino_input",
        "ORF_PREDICTION"       : "MetaPathways_orf_prediction",
        "ORF_TO_AMINO"         : "MetaPathways_create_amino_sequences",
        "FUNC_SEARCH"          : "MetaPathways_func_search",
        "PARSE_FUNC_SEARCH"    : "MetaPathways_parse_blast",
        "COMPUTE_REFSCORES"    : "MetaPathways_refscore",
        "ANNOTATE_ORFS"        : "MetaPathways_annotate_fast",
        "CREATE_ANNOT_REPORTS" : "MetaPathways_create_reports_fast",
        "GENBANK_FILE"         : "MetaPathways_create_genbank_ptinput",
        "SCAN_rRNA"            : "MetaPathways_rRNA_stats_calculator",
        "SCAN_tRNA"            : "MetaPathways_tRNA_scan",
        "RPKM_CALCULATION"     : "MetaPathways_tpm",
        # executables
        "BLASTP_EXECUTABLE"    : 'blastp',
        "BLASTN_EXECUTABLE"    : 'blastn',
        "BWA_EXECUTABLE"       : 'bwa',
        "FASTDB_EXECUTABLE"    : 'fastdb',
        "FAST_EXECUTABLE"      : 'fastal',
        "PRODIGAL_EXECUTABLE"  : 'pprodigal',
        "SCAN_tRNA_EXECUTABLE" : 'tRNAscan-SE',
        "RPKM_EXECUTABLE"      : 'coverm',
        "NUM_CPUS"             : params['Miscellaneous Arguments']['threads'],
        "REFDBS"               : params['Minimum Required Arguments']['refdb_dir']
    }

    block_mode = True

    try:
         # load the sample information
         print(f"RUNNING MetaPathways: v{__version__}")
         if len(input_output_list):
              for input_file in sorted_input_output_list:
                sample_output_dir = input_output_list[input_file]
                algorithm = mpsteps.get_parameter(params,
                                                  'Functional Annotation Arguments',
                                                  'annotation_algorithm',
                                                  default='FAST').upper()

                s = sampledata.SampleData()
                readsdir = params['Read Mapping Arguments']['readsdir']
                if readsdir:
                  s.setInputOutput(inputFile = input_file, sample_output_dir = sample_output_dir, readsDir=readsdir)  
                else:
                  s.setInputOutput(inputFile = input_file, sample_output_dir = sample_output_dir)
                s.setParameter('algorithm', algorithm)
                s.setParameter('FILE_TYPE', filetypes[input_file][0])
                s.setParameter('SEQ_TYPE', filetypes[input_file][1])
                s.clearJobs()

                if not  path.exists(sample_output_dir):
                   makedirs(sample_output_dir)
                s.prepareToRun()
                samplesData[input_file] = s

              # load the sample information
              mpsteps.run_metapathways(
                   samplesData,
                   sample_output_dir,
                   output_dir,
                   globallogger = globalerrorlogger,
                   command_line_params = command_line_params,
                   params = params,
                   status_update_callback = status_update_callback,
                   config_settings = configs,
                   run_type = run_type,
                   block_mode = block_mode
              )
         else:
              gutils.eprintf("ERROR\tNo valid input files/Or no files specified  to process in folder %s!\n",sQuote(input_fp) )
              globalerrorlogger.printf("ERROR\tNo valid input files to process in folder %s!\n",sQuote(input_fp) )

    except:
       mputils.exit_process(str(traceback.format_exc(10)), logger= globalerrorlogger )



    gutils.eprintf("            ***********                \n")
    gutils.eprintf("INFO : FINISHED PROCESSING THE SAMPLES \n")
    gutils.eprintf("             THE END                   \n")
    gutils.eprintf("            ***********                \n")


def build_db():
    argv = sys.argv
    gutils.eprintf("%-10s:%s\n" % ('COMMAND', ' '.join(argv)))
    parser = dbParser()
    args = parser.parse_args()

    gutils.eprintf("Building Refenence DB:")
    mc_list = ['stage_fast_full', 'stage_fast_lite',
               'stage_blast_full', 'stage_blast_lite'
               ]
    if args.db_type in mc_list:
        gutils.eprintf("Your selected database type requires a MetaCyc download...\n")
        # Get the username and password securely
        metacyc_user = input("Enter your MetaCyc username: ")
        metacyc_pswd = getpass.getpass("Enter your MetaCyc password: ")
    else:
        metacyc_user = ''
        metacyc_pswd = ''

    db_dir = args.output_dir
    isExist = path.exists(db_dir)
    if not isExist:
        makedirs(db_dir)

    add_netrc_entry("brg-files.ai.sri.com", f"{metacyc_user}", f"{metacyc_pswd}")

    gutils.eprintf(f' metapathways-data-install.sh {args.output_dir} {args.db_type} {args.num_cpus} \n')
    system(f'metapathways-data-install.sh {args.output_dir} {args.db_type} {args.num_cpus}')


def add_netrc_entry(machine, login, password):
    # Get the user's home directory
    home_directory = path.expanduser("~")

    # Define the path to the .netrc file
    netrc_path = path.join(home_directory, ".netrc")

    # Check if the .netrc file already exists
    if path.exists(netrc_path):
        # Read the existing .netrc file
        with open(netrc_path, 'r') as netrc_file:
            existing_entries = netrc_file.read()
            # Check if an entry for the specified machine already exists
            if f"machine {machine}" in existing_entries:
                print(f"An entry for '{machine}' already exists in .netrc. Not adding a duplicate entry.")
                return

    # If it doesn't exist or the entry doesn't exist, open the .netrc file in append mode
    with open(netrc_path, 'a') as netrc_file:
        # Write the new entry to the file
        netrc_file.write(f"machine {machine}\n")
        netrc_file.write(f"login {login}\n")
        netrc_file.write(f"password {password}\n")


def mag_split():
    argv = sys.argv
    gutils.eprintf("%-10s:%s\n" % ('COMMAND', ' '.join(argv)))
    parser = msParser()
    args = parser.parse_args()

    gutils.eprintf("Mapping Ptools inputs to MAGs:")
    pf_file = path.join(args.output_dir, 'ptools/0.pf')
    orf_map = path.join(args.output_dir, 'ptools/orf_map.txt')
    orf_contig_map = glob.glob(path.join(args.output_dir, 'results/annotation_table/*.ORF_annotation_table.txt'))[0]
    contig_map = glob.glob(path.join(args.output_dir, 'preprocessed/*.mapping.txt'))[0]
    mag_map = args.mag_map
    ms_outdir = path.join(args.output_dir, 'magsplitter')

    isExist = path.exists(ms_outdir)
    if not isExist:
       makedirs(ms_outdir)

    gutils.eprintf(f' magsplitter -p {pf_file} -r {orf_map} -c {orf_contig_map} -m {mag_map} -i {contig_map} -o {ms_outdir} \n')
    system(f'magsplitter -p {pf_file} -r {orf_map} -c {orf_contig_map} -m {mag_map} -i {contig_map} -o {ms_outdir}')


def ptools():
    argv = sys.argv
    gutils.eprintf("%-10s:%s\n" % ('COMMAND', ' '.join(argv)))
    parser = ptParser()
    args = parser.parse_args()
    if args.tag:
        tag = args.tag
    else:
        tag = path.basename(args.output_dir.rstrip('/'))
    gutils.eprintf("Building ePGDBs:")
    gutils.eprintf(f' pgdb_build_wf.py --mp_out {args.output_dir} --tag {tag} \n')
    system(f'pgdb_build_wf.py --mp_out {args.output_dir} --tag {tag}')


def help():
    print(f"""\
        MetaPathways: v{__version__}
        https://metapathways.readthedocs.io
        https://bitbucket.org/BCB2/metapathways

        Syntax: MetaPathways COMMAND [OPTIONS]

        Where COMMAND is one of :
            help
            version
            build_db
            run
            mag_split
            ptools

        for addional help, use:
            MetaPathways COMMAND -h
        """)

def version():
    print(f"MetaPathways v{__version__}")


def main():
    if len(sys.argv) <= 1:
        help()
        return
    {
        "help": help,
        "version": version,
        "build_db" : build_db,
        "run": run,
        "mag_split": mag_split,
        "ptools": ptools
    }.get(
        sys.argv[1],
        help # default
    )()

# the main function of metapaths
if __name__ == "__main__":
    main()
