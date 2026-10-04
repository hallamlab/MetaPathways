"""The main script that calls the pipeline """

try:
    import sys
    import traceback
    import re
    import inspect
    import shutil
    import glob
    import getpass
    import pathlib

    import subprocess
    import threading
    import os
    import signal
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
    from metapathways import nextflow
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



def runParser(command='run'):
    parser = argparse.ArgumentParser(description='MetaPathways command-line tool for annotation of contigs.')
    subparsers = parser.add_subparsers(dest="command")
    run_parser = subparsers.add_parser(
        command,
        description='Minimum REQUIRED Command:\n'
                    'MetaPathways run -i INPUT_FILE -o OUTPUT_DIR -d REFDB_DIR\n\n',
        usage=f'Metapathways {command} [options]', formatter_class=argparse.RawTextHelpFormatter)
    
    # Minimum required args
    req_args = run_parser.add_argument_group('Minimum Required Arguments')
    req_args.add_argument("-i", "--input_file", required=False, default=None,
                        help='path to the input fasta file/input dir [REQUIRED]')
    req_args.add_argument("-o", "--output_dir", required=False, default=None,
                        help='path to the output directory [REQUIRED]')
    req_args.add_argument("-d", "--refdb_dir", required=False, default=None,
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
    func_args.add_argument('--annotation_dbs', nargs='+', type=str, default=['swissprot'],
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

    # Read mapping
    reads_args = run_parser.add_argument_group('Read Mapping Arguments (single sample support only)')
    reads_args.add_argument("-1", "--fastq", dest="fwd_fastq", 
                        help="location of the raw fastq file, either forward or interleaved")
    reads_args.add_argument("-2", "--rev_fastq", dest="rev_fastq",
                        help="location of the raw reverse fastq file, if separate paired-end")
    reads_args.add_argument("--interleaved", action="store_true", default=False,
                        help="if paired-end is interleaved [False]")

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

    nextflow.add_resources(run_parser)
    run_parser.add_argument('--dryrun', action='store_true', help='show the execution plan without running tasks')

    # Other arguments
    misc_args = run_parser.add_argument_group('Miscellaneous Arguments')
    misc_args.add_argument("-s", "--samples", nargs='+', action="append", default=[],
                        help="process only specific samples, space-separated list")
    misc_args.add_argument("-t", "--threads", default=8, type=nextflow.positive,
                        help="threads per capable tool [8], capped by --max_cpus; serial stages use one CPU")
    misc_args.add_argument("-v", "--verbose", action="store_true", default=False,
                        help="print more information on the stdout")
    misc_args.add_argument("--test", action="store_true", help="use test values for all arguments")


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

    nextflow.add_resources(mag_parser)
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
    ptools_parser.add_argument('--taxprune', action="store_true", dest="taxprune", default=False,
                             help='Set taxonomic pruning in pathway tools to True')
    from metapathways.pt_taxonomy import add_taxonomy_options
    add_taxonomy_options(ptools_parser)
    ptools_parser.add_argument('--no_transport_inference', action='store_true', help='Disable TIP transport inference (SIF only)')
    ptools_parser.add_argument('--entity', help='Build only community or the specified MAG ID')
    ptools_parser.add_argument('--image', help='Pathway Tools SIF [registered by build_pt]')
    from metapathways.nextflow import add_resources
    add_resources(ptools_parser, task_memory='4 GB')
    return parser


def blParser(DBS_FUNC, DBS_FUNC_DEFAULT, ALIGNERS):
    parser = argparse.ArgumentParser(description='automated database install')
    db = parser.add_argument_group(title="database arguments")
    db.add_argument("-d", "--refdb_dir", metavar="PATH", required=False, default=None,
                    help="path to save the reference DB, [DEFAULT \"./\"]")
    db.add_argument("--func", metavar="CATEGORICAL", nargs='*', required=False, default=DBS_FUNC_DEFAULT,
                    help=f"functional references, select any combination from {DBS_FUNC}, [DEFAULT {DBS_FUNC_DEFAULT}]")
    db.add_argument("-a", "--aligner", required=False, default="fast",
                    help=f"local aligner to index for, select one of {ALIGNERS}, [DEFAULT fast]")
    db.add_argument('--metacyc_source', help='licensed MetaCyc data directory or Pathway Tools SIF [registered SIF when --func includes metacyc]')

    # "options" group
    parser.add_argument("-t", "--threads", metavar="INT", type=nextflow.positive, required=False, default=None,
                        help="total database-build CPU budget [available CPUs]")
    parser.add_argument("--dryrun", action="store_true", default=False, required=False,
                        help="show the database execution plan without running tasks")
    parser.add_argument("--snakemake", nargs='*', required=False, default=[],
                        help="legacy compatibility flags; use the resource flags for new runs")

    parser.add_argument("--test", action="store_true", help="build reviewer SwissProt/SILVA references; use -d for the prepared MPDB directory")

    nextflow.add_resources(parser)
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
        raise ValueError(errmessage)
    return False


def create_an_input_output_pair(input_file, output_dir,  globalerrorlogger=None):
    """ creates an input output pair if input is just an input file """

    input_output = {}

    if not re.search(r'.(fasta|fas|fna|faa|fa|gbk|gff|fasta.gz|fas.gz|fna.gz|faa.gz|fa.gz)$',input_file, re.IGNORECASE):
       return input_output

    shortname = None
    shortname = re.sub('[.](fasta|fas|fna|faa|fa|gbk|gff|fasta.gz|fas.gz|fna.gz|faa.gz|fa.gz)$','',input_file, re.IGNORECASE)
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
    parser = runParser()
    args = parser.parse_args()
    tasks, output_dir = prepare_annotation(args, parser)
    nextflow.launch(tasks, output_dir, args, 'run', dryrun=args.dryrun)
    print('MetaPathways processing complete.')


def prepare_annotation(args, parser):
    # Shared planner: creates tasks without starting any tools.
    if args.test:
        args.refdb_dir = str(pathlib.Path(path.abspath(__file__)).parent.joinpath("regtests/test_db"))
        args.rRNA_refdbs = ["SILVA_SSU_test", "SILVA_LSU_test"]
        args.annotation_dbs = ["swissprot_test"]
        args.input_file = str(pathlib.Path(path.abspath(__file__)).parent.joinpath("regtests/input/k12_test.fasta.gz"))
        args.fwd_fastq = str(pathlib.Path(path.abspath(__file__)).parent.joinpath("regtests/input/k12_small_R1.fastq.gz"))
        args.rev_fastq = str(pathlib.Path(path.abspath(__file__)).parent.joinpath("regtests/input/k12_small_R2.fastq.gz"))
        args.output_dir = "./test"
        args.threads = 1

        # Set required arguments to False when --test is used
        for action in parser._actions:
            if action.required:
                action.required = False

    # Validate required arguments if --test is not used
    else:
        if args.input_file == None:
            parser.error(f"Input is a required argument: \"-i\", \"--input_file\"")
        elif args.output_dir == None:
            parser.error(f"Output path is a required argument: \"-o\", \"--output_dir\"")
        elif args.refdb_dir == None:
            parser.error(f"Reference DB is a required argument: \"-d\", \"--refdb_dir\"")

    for key in ('input_file', 'output_dir', 'refdb_dir', 'fwd_fastq', 'rev_fastq'):
        value = getattr(args, key, None)
        if value and value != 'None':
            setattr(args, key, str(pathlib.Path(value).expanduser().absolute()))
    budget = args.max_cpus or (args.threads if args.executor == 'slurm' else nextflow.local_capacity()[0])
    args.threads = min(args.threads, budget)
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
    if not getattr(args, '_analysis_planning', False):
        print(f"Output directory: {output_dir}", flush=True)
    globalerrorlogger = mputils.WorkflowLogger(mputils.generate_log_fp(output_dir, basefile_name = 'global_errors_warnings'), open_mode='a')
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
          parser.error(f'Input path does not exist: {input_fp}')

    """ these are the subset of sample to process if specified
        in case of an empty subset process all the sample """

    # remove all samples that are not specifed unless sample_subset is empty
    remove_unspecified_samples(input_output_list, sample_subset, globalerrorlogger = globalerrorlogger)

    # add check the config parameters
    sorted_input_output_list = sorted(input_output_list.keys())
    filetypes = gutils.check_file_types(sorted_input_output_list)

    #stop on invalid samples
    if not halt_on_invalid_input(input_output_list, filetypes, sample_subset):
       parser.error('Invalid input sequence format; see the error log.')

    # make sure the sample files are found
    report_missing_filenames(input_output_list, sample_subset, logger=globalerrorlogger)
    parameter =  paramsmod.Parameters()

    config = gutils.someclass()
    config.refdb_dir = params['Minimum Required Arguments']['refdb_dir']
    
    if not diagnoze.staticDiagnose(params, config,  logger = globalerrorlogger):
        gutils.eprintf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        globalerrorlogger.printf("ERROR\tFailed to pass the test for required scripts and inputs before run\n")
        parser.exit(1, "Input or dependency checks failed; see the error log.\n")

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
        "SCAN_tRNA_EXECUTABLE" : 'ptRNAscan.py',
        "RPKM_EXECUTABLE"      : 'coverm',
        "NUM_CPUS"             : params['Miscellaneous Arguments']['threads'],
        "REFDBS"               : params['Minimum Required Arguments']['refdb_dir']
    }

    block_mode = True

    try:
         # load the sample information
         if not getattr(args, '_analysis_planning', False):
              print(f"RUNNING MetaPathways: v{__version__}", flush=True)
         if len(input_output_list):
              for input_file in sorted_input_output_list:
                sample_output_dir = input_output_list[input_file]
                algorithm = mpsteps.get_parameter(params,
                                                  'Functional Annotation Arguments',
                                                  'annotation_algorithm',
                                                  default='FAST').upper()

                s = sampledata.SampleData()
                
                fwd_fq = params['Read Mapping Arguments']['fwd_fastq']
                rev_fq = params['Read Mapping Arguments']['rev_fastq']
                interleaved = params['Read Mapping Arguments']['interleaved']
                fq_files = [[fwd_fq, rev_fq], interleaved]
                if fwd_fq:
                  s.setInputOutput(inputFile = input_file, sample_output_dir = sample_output_dir, fq_files = fq_files)  
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

              tasks = nextflow.annotation_tasks(samplesData, params, configs, args.memory)
              return tasks, output_dir
         else:
              gutils.eprintf("ERROR\tNo valid input files/Or no files specified  to process in folder %s!\n",gutils.sQuote(input_fp) )
              raise ValueError(f'No valid inputs in {input_fp}')

    except Exception as exc:
       parser.exit(1, f"MetaPathways run failed: {exc}\n")



def build_db():
    argv = sys.argv
    DBS_FUNC = "swissprot cazy eggnog uniref50 uniref90 metacyc".split(" ")
    DBS_FUNC_DEFAULT = "swissprot".split(" ")
    ALIGNERS = "fast blast".split(" ")
    parser = blParser(DBS_FUNC, DBS_FUNC_DEFAULT, ALIGNERS)
    args = parser.parse_args(argv[2:])

    # Check for the --test flag and set test values if present
    if args.test:
        args.refdb_dir = args.refdb_dir or pathlib.Path(path.abspath(__file__)).parent.joinpath("regtests/test_db")
        args.func = ["swissprot_test"]
        args.aligner = "fast"
        args.threads = 1


        # Set required arguments to False when --test is used
        for action in parser._actions:
            if action.required:
                action.required = False

    # Validate required arguments if --test is not used
    else:

        missing_required = []
        for action in parser._actions:
            if action.required and getattr(args, action.dest, None) is None:
                missing_required.append(action.dest)

        if missing_required:
            parser.error(f"The following arguments are required: {', '.join(missing_required)}")

    args.refdb_dir = args.refdb_dir or "./"
    input_error = False
    help_printed = False

    def _error(message: str):
        nonlocal input_error, help_printed
        if not help_printed:
            parser.print_help()
            print()
            help_printed = True
        gutils.eprintf(f"Invalid input: {message}")
        input_error = True

    selected_dbs_functional = args.func
    if len(selected_dbs_functional) == 0: _error("no functional references selected")
    for db in selected_dbs_functional:
        if ((db not in DBS_FUNC) & (db != 'swissprot_test')):
            _error(f"unknown functional reference: {db}")
        elif db == 'swissprot_test':
            print("Running test reference DB")
    alinger = args.aligner
    if alinger not in ALIGNERS: _error(f"unknown aligner: {alinger}")

    if input_error: sys.exit(1)

    from metapathways.nf_databases import plan
    # Preserve the old build_db -t total-core limit, while run -t controls searches.
    if args.max_cpus is None and args.threads is not None:
        args.max_cpus = args.threads
    force = False
    for item in args.snakemake:
        key, _, value = item.partition('=')
        if key in ('cores', 'jobs') and value:
            args.max_cpus = nextflow.positive(value)
        elif key in ('dryrun', 'dry-run'):
            args.dryrun = True
        elif key == 'forceall':
            force = True
        elif key not in ('keep-going', 'rerun-incomplete', 'printshellcmds', 'latency-wait'):
            parser.error(f'Legacy --snakemake option {key!r} has no supported Nextflow equivalent; use resource flags')
    try:
        if args.metacyc_source and 'metacyc' not in args.func:
            parser.error('--metacyc_source requires --func metacyc (optionally alongside other databases)')
        tasks = plan(args.refdb_dir, args.func, args.aligner, test=args.test, memory=args.memory,
                     metacyc_source=args.metacyc_source)
        if force:
            for t in tasks:
                t['status'] = 'redo'
        nextflow.launch(tasks, args.refdb_dir, args, 'build_db', dryrun=args.dryrun)
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'build_db: {exc}\n')


def mag_split():
    argv = sys.argv
    parser = msParser()
    args = parser.parse_args()

    gutils.eprintf("Mapping Ptools inputs to MAGs:")
    pf_file = path.join(args.output_dir, 'ptools/0.pf')
    orf_map = path.join(args.output_dir, 'ptools/orf_map.txt')
    orf_contig_map = glob.glob(path.join(args.output_dir, 'results/annotation_table/*.ORF_annotation_table.txt'))[0]
    feature_table = orf_contig_map.removesuffix('.ORF_annotation_table.txt') + '.ptinput.tsv'
    contig_map = glob.glob(path.join(args.output_dir, 'preprocessed/*.mapping.txt'))[0]
    mag_map = args.mag_map
    ms_outdir = path.join(args.output_dir, 'magsplitter')

    isExist = path.exists(ms_outdir)
    if not isExist:
       makedirs(ms_outdir)
    cmd = ['magsplitter', '-p', pf_file, '-r', orf_map, '-c', orf_contig_map, '-m', mag_map, '-i', contig_map, '-o', ms_outdir]
    cmd_str = f' magsplitter -p {pf_file} -r {orf_map} -c {orf_contig_map} -m {mag_map} -i {contig_map} -o {ms_outdir}'
    gutils.eprintf(cmd_str + '\n')

    import shlex
    cmd = [str(pathlib.Path(v).resolve()) if i in (2, 4, 6, 8, 10, 12) else v for i, v in enumerate(cmd)]
    saved_map = str(pathlib.Path(ms_outdir).resolve() / 'contig_to_mag.tsv')
    copy_map = shlex.join([sys.executable, '-c', 'import shutil,sys; from pathlib import Path; a,b=map(Path,sys.argv[1:]); shutil.copyfile(a,b) if a.resolve()!=b.resolve() else None', str(pathlib.Path(mag_map).resolve()), saved_map])
    tasks = [nextflow.task('mag-split', 'Split annotations into MAGs', [shlex.join(cmd), copy_map],
                          [str(pathlib.Path(p).resolve()) for p in (pf_file, orf_map, orf_contig_map, mag_map, contig_map, feature_table)],
                          [str(pathlib.Path(ms_outdir).resolve() / 'results'), saved_map],
                          cpus=1, memory=args.memory, adopt_existing=False,
                          cache_version='authoritative-mag-coordinates-v1')]
    try:
        nextflow.launch(tasks, args.output_dir, args, 'mag_split')
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'mag_split: {exc}\n')

    

def ptools():
    argv = sys.argv
    parser = ptParser()
    args = parser.parse_args()
    from metapathways.pt_container import registered_image
    from metapathways import nextflow
    image = args.image or (None if args.container else registered_image())
    if args.executor == 'slurm' and not image:
        parser.error('Slurm Pathway Tools tasks require a SIF; run build_pt or pass --image')
    if image:
        image = str(pathlib.Path(image).expanduser().resolve())
    if image and not pathlib.Path(image).is_file():
        parser.error(f'Pathway Tools image does not exist: {image}')
    output = pathlib.Path(args.output_dir).resolve()
    tag = args.tag or output.name
    script = shutil.which('pgdb_build_wf.py')
    if not script:
        candidate = pathlib.Path(__file__).resolve().parent.parent / 'dev/pgdb_build_wf.py'
        if not candidate.is_file():
            parser.error('pgdb_build_wf.py is missing; reinstall MetaPathways')
        script = str(candidate)
    import shlex
    tasks = []
    entities = [('community', tag, output / 'ptools', output / 'results/pgdb/community')]
    entities += [(p.name, p.name, p, output / 'results/pgdb/MAGs' / p.name)
                 for p in sorted((output / 'magsplitter/results').glob('*'))
                 if p.is_dir() and 'non_binned' not in p.name]
    if args.no_transport_inference and not image:
        parser.error('--no_transport_inference requires a Pathway Tools SIF')
    if args.entity:
        entities = [entry for entry in entities if entry[0] == args.entity]
        if not entities:
            parser.error('Unknown PGDB entity: ' + args.entity)
    for entity, entity_tag, inputs, results in entities:
        cmd = [sys.executable, script, '--mp_out', str(output), '--tag', tag,
               '--entity', entity]
        if image:
            cmd += ['--image', image]
        elif args.container:
            cmd.append('--container')
        if args.taxprune:
            cmd.append('--taxprune')
        if args.no_transport_inference:
            cmd.append('--no_transport_inference')
        from metapathways.pt_taxonomy import resolve_taxon
        taxon_id = resolve_taxon(args)
        if taxon_id is not None:
            cmd += ['--taxon_id', str(taxon_id)]
        tasks.append(nextflow.task('pgdb-' + entity, entity, [shlex.join(cmd)],
            ([image] if image else []) + [str(inputs), str(output / 'results/annotation_table'),
                str(output / 'preprocessed' / (output.name + '.fasta'))],
            [str(results / (entity_tag + suffix)) for suffix in ('cyc.tar.bz2', '_pwy.tsv', '_pwy2orf.tsv')],
            cpus=1, memory=args.memory, allow_failure=entity != 'community', adopt_existing=False,
            cache_version='sequence-backed-pgdb-trna-names-v3', host_serial=not bool(image)))
        tasks[-1]['fingerprint_inputs'] = tasks[-1]['inputs'] + [str(output / 'orf_prediction' / (output.name + '.cds.gff'))]
    if not image:
        args.max_tasks = 1
        print('Native Pathway Tools tasks are serialized; build_pt enables isolated parallel runs.')
    try:
        nextflow.launch(tasks, output, args, 'ptools')
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'ptools: {exc}\n')
    return

def help():
    print(f"""\
        MetaPathways: v{__version__}
        https://github.com/hallamlab/MetaPathways

        Syntax: MetaPathways COMMAND [OPTIONS]

        Where COMMAND is one of :
            help
            version
            prepare_test
            build_db
            build_pt
            run
            analysis_wf
            mag_split
            ptools
            report

        for addional help, use:
            MetaPathways COMMAND -h
        """)

def version():
    print(f"MetaPathways v{__version__}")


def build_pt():
    from metapathways.pt_container import main as build_image
    try:
        build_image(sys.argv[2:])
    except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
        print(f'build_pt: {exc}', file=sys.stderr)
        sys.exit(1)


def report():
    from metapathways.report_server import main as report_main
    report_main(sys.argv[2:])


def prepare_test():
    from metapathways.reviewer import main as reviewer_main
    reviewer_main(sys.argv[2:])


def analysis_wf():
    from metapathways.analysis_workflow import main as analysis_main
    analysis_main(sys.argv[2:])


def main():
    if len(sys.argv) <= 1:
        help()
        return
    command = sys.argv[1]
    fn = {'prepare_test': prepare_test, 'help': help, 'version': version, 'build_db': build_db, 'build_pt': build_pt,
          'run': run, 'analysis_wf': analysis_wf, 'mag_split': mag_split, 'ptools': ptools, 'report': report}.get(command, help)
    log_dir = None
    if command in ('run', 'analysis_wf', 'mag_split', 'ptools', 'build_db', 'build_pt', 'report') and not any(x in sys.argv for x in ('-h', '--help')):
        probe = argparse.ArgumentParser(add_help=False)
        probe.add_argument('-o', '--output_dir')
        probe.add_argument('-d', '--refdb_dir')
        probe.add_argument('--test', action='store_true')
        options, _ = probe.parse_known_args(sys.argv[2:])
        log_dir = options.refdb_dir or '.' if command == 'build_db' else options.output_dir
        if options.test and command in ('run', 'build_db'):
            log_dir = './test' if command == 'run' else options.refdb_dir or pathlib.Path(__file__).parent / 'regtests/test_db'
        if command == 'build_pt' and not log_dir:
            from metapathways.pt_container import parser as pt_parser
            log_dir = pt_parser().get_default('output_dir')
    def invoke():
        try:
            fn()
            if command in ('run', 'analysis_wf', 'mag_split', 'ptools') and log_dir and '--dryrun' not in sys.argv:
                from metapathways.reporting import build_report
                report_root = pathlib.Path(log_dir).resolve()
                parent_report = report_root.parent / 'reports/schema.json'
                if parent_report.is_file():
                    import json
                    if str(report_root) in json.loads(parent_report.read_text()).get('sample_paths', []):
                        report_root = report_root.parent
                build_report(report_root)
        except (ValueError, RuntimeError, OSError, subprocess.CalledProcessError) as exc:
            print(f'MetaPathways: {exc}', file=sys.stderr)
            sys.exit(1)
        except KeyboardInterrupt:
            print('MetaPathways interrupted; see the retained logs and work directory.', file=sys.stderr)
            sys.exit(130)
        except Exception:
            traceback.print_exc()
            sys.exit(1)
    if log_dir:
        from metapathways.cli_logging import transcript
        with transcript(log_dir, sys.argv):
            invoke()
    else:
        invoke()


# the main function of metapaths
if __name__ == "__main__":
    main()
