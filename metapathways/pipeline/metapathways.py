#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


try:
    import traceback
    from optparse import make_option
    from os import makedirs, path, listdir, remove, rename, _exit
    import os, sys, errno, shutil, re
    from glob import glob
    from datetime import date

    # from metapaths_utils  import pars[s._command_line_parameters

    from metapathways.utils.sysutil import getstatusoutput, pathDelim

    # from metapathways.utils.utils import *, hasInput, createFolderIfNotFound
    from metapathways.utils.utils import *
    from metapathways.parsers.parse import parse_metapaths_parameters
    from metapathways.pipeline.execution import (
        print_commands,
        execute_tasks,
    )
    from metapathways.pipeline.MetaPathways_gather_run_stats import (
        MetaPathways_gather_run_stats,
    )
    from metapathways.utils.metapathways_utils import (
        fprintf,
        printf,
        eprintf,
        remove_existing_pgdb,
        exit_process,
        WorkflowLogger,
        generate_log_fp,
    )

    from metapathways.pipeline.sampledata import *
    from metapathways.pipeline.jobscreator import *
    from metapathways.pipeline.commands import *
    import metapathways.scripts


except:
    print(""" Could not load some user defined  module functions""")
    sys.exit(3)


PATHDELIM = pathDelim()


def copyFile(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc:  # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else:
            raise


def dry_run_status(commands):
    for command in commands:
        printf("%s", command[0])
        if command[4] == True:
            printf("%s", " Required")
        else:
            printf("%s", " Not Required")
    printf("\n")


def get_refdb_name(dbstring):
    dbstring = dbstring.rstrip()
    dbstring = dbstring.lstrip()
    dbstring = dbstring.lower()
    return dbstring


def format_db(formatdb_executable, seqType, raw_sequence_file, formatted_db, algorithm):
    _temp_formatted_db = formatted_db + "__temp__"

    """ format with 4GB file size """
    if algorithm == "BLAST":
        cmd = "%s -dbtype %s --max_file_sz 4294967296  -in %s -out %s" % (
            formatdb_executable,
            seqType,
            raw_sequence_file,
            _temp_formatted_db,
        )

    if algorithm == "LAST":
        # dirname = os.path.dirname(raw_sequence_file)
        cmd = "%s -s 4G -p -c %s  %s" % (
            formatdb_executable,
            _temp_formatted_db,
            raw_sequence_file,
        )

    result = getstatusoutput(cmd)
    temp_fileList = glob(_temp_formatted_db + "*")
    try:
        for tempFile in temp_fileList:
            file = re.sub("__temp__", "", tempFile)
            rename(tempFile, file)

    except:
        return False

    if result[0] == 0:
        return True
    else:
        return False


# convert an input gbk file to fna faa and gff file
def convert_gbk_to_fna_faa_gff(
    input_gbk, output_fna, output_faa, output_gff, config_settings
):
    cmd = "%s  -g %s --output-fna %s --output-faa %s --output-gff %s" % (
        (config_settings["METAPATHWAYS_PATH"] + config_settings["GBK_TO_FNA_FAA_GFF"]),
        input_gbk,
        output_fna,
        output_faa,
        output_gff,
    )
    return cmd


# convert an input gff file to fna faa and gff file
def convert_gff_to_fna_faa_gff(inputs, outputs, config_settings):
    cmd = "%s " % (
        config_settings["METAPATHWAYS_PATH"] + config_settings["GFF_TO_FNA_FAA_GFF"]
    )
    for source, target in zip(inputs, outputs):
        cmd += " --source " + source + " --target " + target
    return cmd


def make_sure_map_file_exists(config_settings, dbname, globallogger=None):
    dbmapFile = (
        config_settings["REFDBS"]
        + PATHDELIM
        + "functional"
        + PATHDELIM
        + "formatted"
        + PATHDELIM
        + dbname
        + "-names.txt"
    )
    seqFilePath = (
        config_settings["REFDBS"] + PATHDELIM + "functional" + PATHDELIM + dbname
    )
    if not doFilesExist([dbmapFile]):
        eprintf("WARNING: Trying to create database map file for %s\n", dbname)
        if globallogger != None:
            globallogger.write(
                "WARNING: Trying to create database map file for %s\n" % (dbname)
            )

        if not doFilesExist([seqFilePath]):
            eprintf(
                "ERROR : You do not even have the raw sequence for Database  %s to format!\n",
                dbname,
            )
            eprintf("      : Make sure you have the file %s\n", seqFilePath)

            if globallogger != None:
                globallogger.write(
                    "ERROR \t You do not even have the raw sequence for Database  %s to format!\n"
                    % (dbname)
                )
                globallogger.write("Make sure you have the file %s\n" % (seqFilePath))

            exit_process()

        mapfile = open(dbmapFile, "w")
        seqFile = open(seqFilePath, "r")
        for line in seqFile:
            if re.match(r">", line):
                fprintf(mapfile, "%s\n", line.strip())
        seqFile.close()
        mapfile.close()

    return dbmapFile


# gets the parameter value from a category as.ecified in the
# parameter file
def get_parameter(params, category, field, default=None):
    if params == None:
        return default

    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default


# parameter file
def get_make_parameter(params, category, field, default=False):
    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default


def get_pipeline_steps(steps_log_file):
    try:
        logfile = open(steps_log_file, "r")
    except IOError:
        eprintf("Did not find %s!\n", logfile)
        eprintf("Try running in 'complete' run-type\n")
    else:
        lines = logfile.readlines()

    pipeline_steps = None
    return pipeline_steps


def write_run_parameters_file(fileName, parameters):
    try:
        paramFile = open(fileName, "w")
    except IOError:
        eprintf("Cannot write run parameters to file %s!\n", fileName)
        exit_process("Cannot write run parameters to file %s" % (fileName))

    #       16s_rRNA      {'min_identity': '40', 'max_evalue': '0.000001', 'min_bitscore': '06', 'refdbs': 'silva_104_rep_set,greengenes_db_DW'}
    paramFile.write("\nRun Date : " + str(date.today()) + " \n")

    paramFile.write("\n\nNucleotide Quality Control parameters[s.n")
    paramFile.write(
        "  min length" + "\t" + str(parameters["quality_control"]["min_length"]) + "\n"
    )

    paramFile.write("\n\nORF prediction parameters[s.n")
    paramFile.write(
        "  min length" + "\t" + str(parameters["orf_prediction"]["min_length"]) + "\n"
    )
    paramFile.write(
        "  algorithm" + "\t" + str(parameters["orf_prediction"]["algorithm"]) + "\n"
    )

    paramFile.write("\n\nAmino acid quality control and annotation parameters[s.n")
    paramFile.write(
        "  min bit score" + "\t" + str(parameters["annotation"]["min_score"]) + "\n"
    )
    paramFile.write(
        "  min seq length" + "\t" + str(parameters["annotation"]["min_length"]) + "\n"
    )
    paramFile.write(
        "  annotation reference dbs"
        + "\t"
        + str(parameters["annotation"]["dbs"])
        + "\n"
    )
    paramFile.write(
        "  min BSR" + "\t" + str(parameters["annotation"]["min_bsr"]) + "\n"
    )
    paramFile.write(
        "  max evalue" + "\t" + str(parameters["annotation"]["max_evalue"]) + "\n"
    )

    paramFile.write("\n\nPathway Tools parameters[s.n")
    paramFile.write(
        "  taxonomic pruning "
        + "\t"
        + str(parameters["ptools_settings"]["taxonomic_pruning"])
        + "\n"
    )

    paramFile.write("\n\nrRNA search/match parameters[s.n")
    paramFile.write(
        "  min identity" + "\t" + str(parameters["rRNA"]["min_identity"]) + "\n"
    )
    paramFile.write(
        "  max evalue" + "\t" + str(parameters["rRNA"]["max_evalue"]) + "\n"
    )
    paramFile.write(
        "  rRNA reference dbs" + "\t" + str(parameters["rRNA"]["refdbs"]) + "\n"
    )

    paramFile.close()


# check for empty values in parameter settings
def checkMissingParam_values(params, choices, logger=None):
    reqdCategoryParams = {
        "annotation": {"dbs": False},
        "orf_prediction": {},
        "rRNA": {},
        "metapaths_steps": {},
    }

    success = True
    for category in choices:
        for parameter in choices[category]:
            if (not params[category][parameter]) and (
                (category in reqdCategoryParams)
                and (parameter in reqdCategoryParams[category])
                and reqdCategoryParams[category][parameter]
            ):
                print(category, parameter)
                print(reqdCategoryParams)
                print(reqdCategoryParams[category])
                eprintf(
                    "ERROR: Empty parameter %s of type %s\n" % (parameter, category)
                )
                eprintf("Please select at least one database for %s\n" % (category))
                if logger != None:
                    logger.write(
                        "ERROR\tEmpty parameter %s of type %s\n" % (parameter, category)
                    )
                    logger.write(
                        "Please select at least one database for %s\n" % (category)
                    )
                success = False

    return success


# check if all of the metapaths_steps have
# settings from the valid list [ yes, skip stop, redo]


def checkParam_values(allcategorychoices, parameters, runlogger=None):
    for category in allcategorychoices:
        for choice in allcategorychoices[category]:
            if choice in parameters:

                if not parameters[choice] in allcategorychoices[category][choice]:
                    logger.write("ERROR\tIncorrect setting in your parameter file")
                    logger.write("for step %s as %s" % (choice, parameters[choices]))
                    eprintf(
                        "ERROR: Incorrect setting in your parameter file"
                        + "       for step %s as %s",
                        choice,
                        parameters[choices],
                    )
                    exit_process()


def checkMetapathsteps(params, runlogger=None):
    choices = {"metapaths_steps": {}, "annotation": {}, "INPUT": {}}

    choices["INPUT"]["format"] = [
        "fasta",
        "gbk_unannotated",
        "gbk_annotated",
        "gff_unannotated",
        "gff_annotated",
    ]

    choices["annotation"]["algorithm"] = ["last", "blast"]

    choices["metapaths_steps"]["PREPROCESS_FASTA"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["ORF_PREDICTION"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["GFF_TO_AMINO"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["FILTERED_FASTA"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["COMPUTE_REFSCORE"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["BLAST_REFDB"] = ["yes", "skip", "stop", "redo", "grid"]
    choices["metapaths_steps"]["PARSE._BLAST"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["SCAN_rRNA"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["STATS_rRNA"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["ANNOTATE"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["PATHOLOGIC_INPUT"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["GENBANK_FILE"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["CREATE_SEQUIN_FILE"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["CREATE_REPORT_FILES"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["SCAN_tRNA"] = ["yes", "skip", "stop", "redo"]
    choices["metapaths_steps"]["PATHOLOGIC"] = ["yes", "skip", "stop", "redo"]

    if params["metapaths_steps"]:
        checkParam_values(choices, params["metapaths_steps"], runlogger)

    checkparams = {}
    checkparams["annotation"] = []
    checkparams["annotation"].append("dbs")

    if not checkMissingParam_values(params, checkparams, runlogger):
        exit_process("Missing parameters")


def copy_fna_faa_gff_orf_prediction(source_files, target_files, config_settings):

    for source, target in zip(source_files, target_files):

        sourcefile = open(source, "r")
        targetfile = open(target, "w")
        sourcelines = sourcefile.readlines()
        for line in sourcelines:
            fprintf(targetfile, "%s\n", line.strip())

        sourcefile.close()
        targetfile.close()


#################################################################################
########################### BEFORE BLAST ########################################
#################################################################################
def run_metapathways(
    samplesData,
    output_dir,
    all_samples_output_dir,
    globallogger,
    command_line_params,
    params,
    status_update_callback,
    run_type,
    config_settings=None,
    block_mode=False,
):
    runid = 'random'
    jobcreator = JobCreator(params, config_settings)

    sorted_samplesData_keys = sorted(samplesData.keys())
    for input_file in sorted_samplesData_keys:
        s = samplesData[input_file]
        jobcreator.addJobs(s, block_mode=block_mode)

    _params = Singleton(Params)(params)

    if block_mode:
        eprintf("==============  RUNNING STEPS IN BLOCK 0 ================\n")
        for input_file in sorted_samplesData_keys:
            s = samplesData[input_file]
            s.stepslogger.printf(
                "\n\n==============  BEGIN RUN "
                + s.sample_name
                + " "
                + runid
                + " BLOCK0 ================\n"
            )
            sample_name_banner = "PROCESSING INPUT " + input_file
            eprintf("\n" + "#" * len(sample_name_banner) + "\n")
            eprintf("\n" + sample_name_banner + " [STEPS BLOCK 0] " + "\n")
            s.writeParamsToRunLogs(_params)
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=0)
            except:
                print(traceback.print_exc(10))
                pass

        for input_file in sorted_samplesData_keys:
            s = samplesData[input_file]
            s.stepslogger.printf(
                "\n\n==============  BEGIN RUN "
                + s.sample_name
                + " "
                + runid
                + " BLOCK1 ================\n"
            )
            sample_name_banner = "PROCESSING INPUT " + input_file
            eprintf("\n" + "#" * len(sample_name_banner) + "\n")
            eprintf("\n" + sample_name_banner + " [STEPS BLOCK 1] " + "\n")
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=1)
            except:
                pass

        for input_file in sorted_samplesData_keys:
            s = samplesData[input_file]
            s.stepslogger.printf(
                "\n\n==============  BEGIN RUN "
                + s.sample_name
                + " "
                + runid
                + " BLOCK2 ================\n"
            )
            sample_name_banner = "PROCESSING INPUT " + input_file
            eprintf("\n" + "#" * len(sample_name_banner) + "\n")
            eprintf("\n" + sample_name_banner + " [STEPS BLOCK 2] " + "\n")
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=2)
            except:
                pass

    else:
        for input_file in sorted_samplesData_keys:
            s = samplesData[input_file]
            s.stepslogger.printf(
                "\n\n==============  BEGIN RUN "
                + s.sample_name
                + " "
                + runid
                + "  ==================\n"
            )
            sample_name_banner = "PROCESSING INPUT " + input_file
            eprintf("#" * len(sample_name_banner) + "\n")
            eprintf("\n" + sample_name_banner + "\n")
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=0)
            except:
                pass
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=1)
            except:
                pass
            try:
                execute_tasks(s, verbose=command_line_params["verbose"], block=2)
            except:
                pass

    return
