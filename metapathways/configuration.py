"""Contains general utility code for the metapaths project"""

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    from shutil import rmtree
    from optparse import make_option
except:
    print("Cannot load some modules")
    print(traceback.print_exc(10))
    sys.exit(0)


class Configuration:
    """Contains the acceptable configurations of the tool
    If you want to add any new module, you will have to add the new algorithm
    this configuration file."""

    acceptableConfiguration = {}
    missingErrors = {}


    def __init__(self):
        try:
            self.file = self.initializeConfiguration()
        except IOError:
            print('ERROR : Cannot  initialize "Configuration" in file ' + sys.argv[0])

    def isValid1(self, key1):
        if key1 in self.acceptableConfiguration:
            return True
        return False

    def getValidValues1(self, key1):
        if self.isValid1(key1):
            return self.acceptableConfigurations.keys()
        return []

    def isValid2(self, key1, key2):
        if not self.isValid1(key1):
            return False

        if not key2 in self.acceptableConfigurations[key1]:
            return False

        return True

    def getValidValues2(self, key1, key2):
        if self.isValid2(key1, key2):
            return self.acceptableConfigurations[key1][key2].keys()
        return []

    def getConfiguration(self):
        return self.acceptableConfiguration

    def initializeConfiguration(self):
        self.acceptableConfiguration = {
            "PYTHON_EXECUTABLE": True,
            "PGDB_FOLDER": False,
            "METAPATHWAYS_PATH": True,
            "PATHOLOGIC_EXECUTABLE": False,
            "REFDBS": False,
            "FORMATDB_EXECUTABLE": False,
            "BLASTP_EXECUTABLE": False,
            "BLASTN_EXECUTABLE": False,
            "EXECUTABLES_DIR": True,
            "RESOURCES_DIR": True,
            "LASTDB_EXECUTABLE": False,
            "LAST_EXECUTABLE": False,
            "ORF_PREDICTION": False,
            "SCAN_tRNA": False,
            "GBK_TO_FNA_FAA_GFF": True,
            "GFF_TO_FNA_FAA_GFF": True,
            "PREPROCESS_INPUT": True,
            "ORF_TO_AMINO": True,
            "COMPUTE_REFSCORES": True,
            "PARSE_FUNC_SEARCH": True,
            "ANNOTATE_ORFS": True,
            "GENBANK_FILE": True,
            "CREATE_ANNOT_REPORTS": True,
            "STATS_rRNA": True,
            "SCAN_tRNA": False,
            "COMPUTE_RPKM": True,
            "NUM_CPUS": False,
        }

        self.missingErrors = {
            "PYTHON_EXECUTABLE": True,
            "PGDB_FOLDER": False,
            "METAPATHWAYS_PATH": True,
            "PATHOLOGIC_EXECUTABLE": False,
            "REFDBS": False,
            "FORMATDB_EXECUTABLE": False,
            "BLASTP_EXECUTABLE": False,
            "BLASTN_EXECUTABLE": False,
            "EXECUTABLES_DIR": True,
            "RESOURCES_DIR": True,
            "LASTDB_EXECUTABLE": False,
            "LAST_EXECUTABLE": False,
            "ORF_PREDICTION": False,
            "SCAN_tRNA": False,
            "GBK_TO_FNA_FAA_GFF": True,
            "GFF_TO_FNA_FAA_GFF": True,
            "PREPROCESS_INPUT": True,
            "ORF_TO_AMINO": True,
            "COMPUTE_REFSCOREs": True,
            "PARSE_BLAST": True,
            "ANNOTATE_ORFS": True,
            "GENBANK_FILE": True,
            "CREATE_ANNOT_REPORTS": True,
            "STATS_rRNA": True,
            "SCAN_tRNA": False,
            "COMPUTE_RPKM": True,
            "NUM_CPUS": True,
        }


if __name__ == "__main__":
    v = Configuration()
