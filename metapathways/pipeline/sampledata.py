#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"


"""Contains general utility code for the metapaths project"""

try:
    import os, re

    from shutil import rmtree
    from optparse import make_option
    from os import path, _exit

    from metapathways.utils.metapathways_utils import *
    from metapathways.utils.sysutil import pathDelim
    from metapathways.utils.utils import *
except:
    print("Cannot load some modules")
    sys.exit(0)

PATHDELIM = pathDelim()


class SampleData:
    """Contains the sample related data """

    runlogger = None
    stepslogger = None
    errorlogger = None
    runstatslogger = None

    input_file = None
    output_dir = None
    sample_name = None

    preprocessed_dir = None
    orf_prediction_dir = None
    genbank_dir = None
    output_run_statistics_dir = None
    blast_results_dir = None
    output_results = None
    output_results_annotation_table_dir = None
    output_results_megan_dir = None
    output_results_sequin_dir = None
    output_results_rpkm_dir = None
    output_fasta_pf_dir = None
    output_results_pgdb_dir = None
    output_results_rRNA_dir = None
    output_results_tRNA_dir = None
    ncbi_params_file = None
    ncbi_sequin_sbt = None
    bwa_folder = None

    stages = []
    stages_context = {}

    def __init__(self):
        pass

    def getContexts(self):
        contexts = []
        for name in self.stages:
            if name in self.stages_context:
                contexts.append(self.stages_context[name])
        return contexts

    def getContextBlocks(self):
        contextBlocks = []

        for block in self.stages:
            contexts = []
            for name in block:
                if name in self.stages_context:
                    contexts.append(self.stages_context[name])
            contextBlocks.append(contexts)

        return contextBlocks

    def numJobs(self):
        return len(self.stages)

    def clearJobs(self):
        self.stages = []
        self.stages_context = {}

    def setInputOutput(self, inputFile=None, sample_output_dir=None):
        if inputFile == None and sample_output_dir == None:
            return False

        self.input_file = inputFile
        self.sample_name = re.sub(r"[.][a-zA-Z]*$", "", self.input_file)
        self.sample_name = path.basename(self.sample_name)
        self.sample_name = re.sub("[.]", "_", self.sample_name)

        self.rpkm_input_dir = path.dirname(inputFile) + PATHDELIM + "reads"

        self.output_dir = sample_output_dir

        self.preprocessed_dir = self.output_dir + PATHDELIM + "preprocessed" + PATHDELIM
        self.orf_prediction_dir = (
            self.output_dir + PATHDELIM + "orf_prediction" + PATHDELIM
        )
        self.genbank_dir = self.output_dir + PATHDELIM + "genbank" + PATHDELIM
        self.output_run_statistics_dir = (
            self.output_dir + PATHDELIM + "run_statistics" + PATHDELIM
        )
        self.blast_results_dir = (
            self.output_dir + PATHDELIM + "blast_results" + PATHDELIM
        )
        self.bwa_folder = self.output_dir + PATHDELIM + "bwa" + PATHDELIM
        self.output_results = self.output_dir + PATHDELIM + "results" + PATHDELIM
        self.output_results_annotation_table_dir = (
            self.output_results + PATHDELIM + "annotation_table" + PATHDELIM
        )

        self.output_results_megan_dir = (
            self.output_results + PATHDELIM + "megan" + PATHDELIM
        )
        self.output_results_rpkm_dir = self.output_results + PATHDELIM + "rpkm"
        self.output_fasta_pf_dir = self.output_dir + PATHDELIM + "ptools" + PATHDELIM
        self.output_results_pgdb_dir = (
            self.output_results + PATHDELIM + "pgdb" + PATHDELIM
        )
        self.output_results_rRNA_dir = (
            self.output_results + PATHDELIM + "rRNA" + PATHDELIM
        )
        self.output_results_tRNA_dir = (
            self.output_results + PATHDELIM + "tRNA" + PATHDELIM
        )
        self.run_stats_file = (
            self.output_run_statistics_dir
            + PATHDELIM
            + self.sample_name
            + ".run.stats.txt"
        )

    def setParameter(self, parameter, value):
        setattr(self, parameter, value)

    def prepareToRun(self):
        self._createFolders()
        self._createLogFiles()

    def getType(self):
        if not hasattr(self, "SEQ_TYPE"):
            return None

        if not hasattr(self, "FILE_TYPE"):
            return None

        if self.FILE_TYPE == "FASTA":
            if self.SEQ_TYPE == "NUCL":
                return self.SEQ_TYPE + "-" + self.FILE_TYPE

            if self.SEQ_TYPE == "AMINO":
                return self.SEQ_TYPE + "-" + self.FILE_TYPE

        if self.FILE_TYPE == "GENBANK":
            if self.SEQ_TYPE == "NOT-USED":
                return "AMINO" "-" + self.FILE_TYPE + "-" + "UNANNOT"

        return "UNKNOWN"

    def _createLogFiles(self):
        self.runlogger = WorkflowLogger(
            generate_log_fp(self.output_dir, basefile_name="metapathways_run_log"),
            open_mode="w",
        )
        self.stepslogger = WorkflowLogger(
            generate_log_fp(self.output_dir, basefile_name="metapathways_steps_log"),
            open_mode="a",
        )
        self.errorlogger = WorkflowLogger(
            generate_log_fp(self.output_dir, basefile_name="errors_warnings_log"),
            open_mode="a",
        )
        self.runstatslogger = WorkflowLogger(
            generate_log_fp(
                self.output_run_statistics_dir,
                basefile_name=self.sample_name + ".run.stats",
            ),
            open_mode="a",
        )

    def writeParamsToRunLogs(self, params):
        param_values = [
            ["quality_control", "min_length", "180"],
            ["quality_control", "delete_replicates", "yes"],
            ["orf_prediction", "strand", "both"],
            ["orf_prediction", "algorithm", "prodigal"],
            ["orf_prediction", "min_length", "60"],
            ["orf_prediction", "translation_table", "11"],
            ["orf_prediction", "mode", "meta"],
            ["annotation", "algorithm", "BLAST"],
            ["annotation", "dbs_high", ""],
            ["annotation", "dbs_custom", ""],
            ["annotation", "dbs", ""],
            ["annotation", "dbtype", "high"],
            ["annotation", "min_bsr", "0.4"],
            ["annotation", "max_evalue", "0.000001"],
            ["annotation", "min_score", "20"],
            ["annotation", "min_length", "45"],
            ["annotation", "max_hits", "5"],
            [
                "rRNA",
                "refdbs",
                "SILVA_128_SSURef_tax_silva,SILVA_128_LSURef_tax_silva, GREENGENES_gg16S_13_5",
            ],
            ["rRNA", "max_evalue", "0.000001"],
            ["rRNA", "min_identity", "20"],
            ["rRNA", "min_bitscore", "50"],
            ["ptools_settings", "taxonomic_pruning", "no"],
            ["ptools_input", "compact_mode", "yes"],
        ]
        for parameter in param_values:
            value = params.get(parameter[0], parameter[1], default=parameter[2])
            self.runlogger.printf("%s:%s\t%s\n", parameter[0], parameter[1], value)

    def _createFolders(self):
        checkOrCreateFolder(self.preprocessed_dir)
        checkOrCreateFolder(self.orf_prediction_dir)
        checkOrCreateFolder(self.genbank_dir)
        checkOrCreateFolder(self.output_run_statistics_dir)
        checkOrCreateFolder(self.blast_results_dir)
        checkOrCreateFolder(self.bwa_folder)
        checkOrCreateFolder(self.output_results)
        checkOrCreateFolder(self.output_results_annotation_table_dir)
        checkOrCreateFolder(self.output_results_megan_dir)
        checkOrCreateFolder(self.output_results_rpkm_dir)
        eprintf("Creating ========== %s\n", self.output_fasta_pf_dir)
        checkOrCreateFolder(self.output_fasta_pf_dir)
        eprintf("Creating ======= ============  ===== 3.1 %s\n", path.exists("/home/kishori/mp_output/testsample1/ptools")) 
        checkOrCreateFolder(self.output_results_pgdb_dir)
        checkOrCreateFolder(self.output_results_rRNA_dir)
        checkOrCreateFolder(self.output_results_tRNA_dir)

    def addPipeLineStage(self, stepName, inputs=[], outputs=[], status="yes"):

        stagecontext = Context()
        stagecontext.inputs = inputs
        stagecontext.outputs = outputs
        stagecontext.name = stepName
        stagecontext.status = context

        stages_context[stepName] = stagecontext

    def addContextBlock(self, contextBlock):
        self.contextBlock.append(contextBlock)

    def addContexts(self, contextBlock):
        stages = []

        for contexts in contextBlock:
            for context in contexts:
                stages.append(context.name)
                self.stages_context[context.name] = context

        self.stages.append(stages)

    def hasPToolsInput(self):
        """ checks if the ptools folder has the right inputs"""
        files = ["0.pf", "0.fasta", "genetic-elements.dat", "organism-params.dat"]

        for file in files:
            if not doesFileExist(self.output_fasta_pf_dir + PATHDELIM + file):
                return False

        return True

    def hasGenbankFile(self):
        """ checks if genbank file has already been created"""
        output_annot_gbk = self.genbank_dir + PATHDELIM + self.sample_name + ".gbk"
        return doesFileExist(output_annot_gbk)

    def hasSequinFile(self):
        """ checks if sequin file has already been created"""
        output_annot_sequin = (
            self.output_results_sequin_dir + PATHDELIM + self.sample_name + ".tbl"
        )
        return doesFileExist(output_annot_sequin)
