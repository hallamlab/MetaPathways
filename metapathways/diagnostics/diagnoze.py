#!/usr/bin/env python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2013, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

"""Contains general utility code for the metapaths project"""

try:
    import sys
    from shutil import rmtree
    from optparse import make_option
    from metapathways.diagnostics.parameters import *
    from metapathways.diagnostics.configuration import *
    from metapathways.utils.sysutil import pathDelim, getstatusoutput
    from metapathways.utils.utils import *
    from metapathways.utils.errorcodes import *
    from metapathways.utils.metapathways_utils import fprintf
    from os import path, _exit, rename, remove
except:
    print("Cannot load some modules")
    sys.exit(0)

PATHDELIM = pathDelim()


def staticDiagnose(params, config = None, logger=None):
    """
    Diagnozes the pipeline basedon the configs and params for
    binaries, scripts and resources
    """

    """ makes sure that the choices in  parameter file are valid """
    errors = checkParams(params, logger=logger)
    if errors:
        return False

    """ Get the configurations for the executables/scripts and databasess """
    _configuration = Configuration()
    configuration = _configuration.getConfiguration()

    """ the place holders for the tools required to make the run """
    parameters = Parameters()

    """ check if the required standard databases exists """

    #    print  parameters.getRunSteps( activeOnly = True)
    if not checkForRequiredDatabases(
        params, config, "functional", logger=logger):
        insert_error(17)
        return False

    if not checkForRequiredDatabases(
        params, config, "taxonomic", logger=logger):
        insert_error(17)
        return False

    """ make sure all the executables exist """
    executables = [ 'fastal', 'fastdb', 'rpkm' ]
    message, ok = checkbinaries(executables)

    if message:
        print(message)
        return False

    return True


def checkbinaries(executables):
    message = ""
    ok = True
    cmd_exists = lambda x: shutil.which(x) is not None
    for executable in executables:
        if not cmd_exists(executable):
            message = "'" + executable + "'" + " is missing.\n"
            ok = False

    return message, ok

def checkForRequiredDatabases(params, config,  dbType, logger=None):
    """checks the
    -- database folder structure
    -- checks for raw sequences
    -- checks for formatted sequences
    -- formats if necessary
    """

    if dbType == "functional":
        dbtype = get_parameter(params, "annotation", "dbtype", default="high")
        dbstring = ""
        _algorithm = get_parameter(params, "annotation", "algorithm", default=None)

    if dbType == "taxonomic":
        dbstring = get_parameter(params, "rRNA", "refdbs", default=None)
        _algorithm = get_parameter(params, "annotation", "algorithm", default=None)

    if dbstring == None:
        eprintf(
            "WARNING\tReference databases to annotate with is unspecified, please add it in the params file\n"
        )
        return False

    dbs = [x.strip() for x in dbstring.split(",") if len(x) != 0]

    if not dbs:
        return True

    """ checks refdb path """
    if not check_if_refDB_path_valid(config.refdb_dir, logger=logger):
        return False

    """ checks raw sequences for dbtype functional/taxonimic """
    if isRefDBNecessary(params, dbType):
        if not check_for_raw_sequences(dbs, config.refdb_dir, dbType, logger=logger):
            insert_error(17)
            return False

        for db in dbs:
            algorithm = ""
            if dbType == "taxonomic":
                algorithm = _algorithm
                seqType = "nucl"
            elif dbType == "functional":
                algorithm = _algorithm
                seqType = "prot"
            else:
                algorithm = None

            """ is db formatted ? """
            if not isDBformatted(db, config.refdb_dir, dbType, seqType, \
                algorithm, logger=logger
                ):
                """ if note formatted then format it """
                eprintf(
                    "WARNING\tTrying to format %s  database %s\n", seqType, sQuote(db)
                )
                logger.printf(
                    "WARNING\tTrying to format %s database %s \n", seqType, sQuote(db)
                )

                return False


            seqFilePath = configs["REFDBS"] + PATHDELIM + dbType + PATHDELIM + db
            """ check for dbmapfile """

            if not doesFileExist(dbMapFile):
                eprintf(
                    "WARNING\tDoes not have map file %s for %s\n",
                    sQuote(dbMapFile),
                    sQuote(db),
                )
                logger.printf(
                    "WARNING\tDoes not have map file %s for %s\n",
                    sQuote(dbMapFile),
                    sQuote(db),
                )
                if not createMapFile(seqFilePath, dbMapFile):
                    eprintf(
                        "ERROR\tFailed to create map file %s for %s\n",
                        sQuote(dbMapFile),
                        sQuote(db),
                    )
                    logger.printf(
                        "ERROR\tFailed to create map file %s for %s\n",
                        sQuote(dbMapFile),
                        sQuote(db),
                    )
                    return False
                eprintf(
                    "INFO\tSuccessfully created  map file %s for %s\n",
                    sQuote(dbMapFile),
                    sQuote(db),
                )
                logger.printf(
                    "INFO\tSuccessfully created map file %s for %s\n",
                    sQuote(dbMapFile),
                    sQuote(db),
                )

    return True


def isRefDBNecessary(params, dbType):
    """ decide yes or no based on the params settings yes or redo """
    if dbType == "functional":
        status = get_parameter(params, "metapaths_steps", "FUNC_SEARCH", default=None)
        if status in ["yes", "redo"]:
            return True

    if dbType == "taxonomic":
        status = get_parameter(params, "metapaths_steps", "SCAN_rRNA", default=None)
        if status in ["yes", "redo"]:
            return True

    return False


def isDBformatted(db, refdbspath, dbType, seqType, algorithm, logger=None):
    """ check if the DB is formatted """
    """Checks if the formatted database for the specified algorithm exits """
    dbPath = refdbspath + PATHDELIM + dbType + PATHDELIM + "formatted"
    dbname = dbPath + PATHDELIM + db
    suffixes = getSuffixes(algorithm, seqType)

    # print algorithm, suffixes
    if not suffixes:
        return False

    status = False
    for suffix in suffixes:
        allfileList = glob(dbname + "*." + suffix)

        fileList = []
        tempFilePattern = re.compile(r"" + dbname + "[.\d]*." + suffix + "$")

        for aFile in allfileList:
            searchResult = tempFilePattern.search(aFile)
            if searchResult:
                fileList.append(aFile)

        if len(fileList) == 0:
            eprintf("WARNING\tsequence for db  %s not formatted\n", dbname)
            logger.printf("WARNING\tsequence for db  %s not formatted\n", dbname)
            return False

        status = True

    return status


def check_if_refDB_path_valid(refdbspath, logger=None):
    """it checks for the validity of the refdbs path structure
    refdbpath  /functional
                   /formatted
               /tanxonomic
                   /formatted
    """

    status = True
    if not doesFolderExist(refdbspath):
        eprintf("ERROR\treference sequence folder %s not found\n", sQuote(refdbspath))
        logger.printf(
            "ERROR\treference sequence folder %s not found\n", sQuote(refdbspath)
        )
        return False

    dbTypes = ["functional", "taxonomic"]
    """ now check if respective dbtype folders are available """
    status = True
    for dbType in dbTypes:
        if not doesFolderExist(refdbspath + PATHDELIM + dbType):
            eprintf(
                "ERROR\tfolder %s for reference type %s not found\n",
                sQuote(refdbspath + PATHDELIM + dbType),
                dbType,
            )
            logger.printf(
                "ERROR\tfolder %s for reference type %s not found\n",
                sQuote(refdbspath + PATHDELIM + dbType),
                dbType,
            )
            status = False

    if status == False:
        return status

    """ now check if path to drop the formatted dbs are available """
    for dbType in dbTypes:
        if not doesFolderExist(
            refdbspath + PATHDELIM + dbType + PATHDELIM + "formatted"
        ):
            eprintf(
                "ERROR\tsubfolder %s not found under the folder %s\n",
                sQuote("formatted"),
                sQuote(refdbspath + PATHDELIM + dbType + PATHDELIM),
            )
            logger.printf(
                "ERROR\tsubfolder %s not found under the folder %s\n",
                sQuote("formatted"),
                sQuote(refdbspath + PATHDELIM + dbType + PATHDELIM),
            )
            status = False

    return status


def check_for_raw_sequences(dbs, refdbspath, dbType, logger=None):
    """ check for the raw sequence file """
    status = True
    for db in dbs:
        fullPath = refdbspath + PATHDELIM + dbType + PATHDELIM + db
        if not does_plain_or_gz_FileExist(fullPath):
            eprintf(
                "ERROR\tRaw sequences %s expected for %s references\n", fullPath, dbType
            )
            logger.printf(
                "ERROR\tRaw sequences %s expected for %s references\n", fullPath, dbType
            )
            status = False
    return status


def get_parameter(params, category, field, default=None):
    """gets the parameter value from a category
    as specified in the  parameter file"""

    if params == None:
        return default

    if category in params:
        if field in params[category]:
            return params[category][field]
        else:
            return default
    return default


def _checkParams(params, paramsAccept, logger=None, errors=None):

    """make sure that every parameter in the params is valid recursively
    This is initialed by the checkParams() function
    store the erros in the erros dictionary
    """
    """ if not level to go deeper  then the leaves of the dict are reached"""

    if not type(params) is dict and type(paramsAccept) is dict:
        # print  'type ',  params, paramsAccept,  (not params in paramsAccept), (len(paramsAccept.keys())!=0)
        try:
            if (not params in paramsAccept) and len(paramsAccept.keys()) != 0:
                errors[params] = False
                choices = ", ".join(paramsAccept.keys())
                eprintf(
                    "ERROR\tValue for key %s, in param file,  is not set propertly must be one of %s \t %s\n",
                    sQuote(params),
                    sQuote(choices),
                    __name__,
                )
                logger.printf(
                    "ERROR\tValue for key %s, in param file,  is not set propertly must be one of %s\t%s\n",
                    sQuote(params),
                    sQuote(choices),
                    __name__,
                )
        except:
            pass
        return

    """  make sure that every parameter in the params is valid recursively """
    for key, value in params.items():
        if type(paramsAccept) is dict:
            if len(key) and key in paramsAccept:
                _checkParams(
                    params[key], paramsAccept[key], logger=logger, errors=errors
                )


def checkParams(params, logger=None):
    """ makes sure that all the params provides are valid or acceptable """
    """ when the choices are not any of the acceptable 
    values then it is considered erroneous"""

    _paramsAccept = Parameters()
    paramsAccept = _paramsAccept.getAcceptableParameters()
    errors = {}

    for key, value in params.items():
        if key in paramsAccept:
            _checkParams(params[key], paramsAccept[key], logger=logger, errors=errors)

    return errors


def getSuffixes(algorithm, seqType):
    """Get the suffixes for the right algorithm with the right
    sequence type
    """

    suffixes = {}
    suffixes["LAST"] = {}
    suffixes["BLAST"] = {}
    suffixes["BLAST"]["nucl"] = ["nhr", "nsq", "nin"]
    suffixes["BLAST"]["prot"] = ["phr", "psq", "pin"]

    suffixes["LAST"]["nucl"] = ["des", "sds", "suf", "bck", "prj", "ssp", "tis"]
    suffixes["LAST"]["prot"] = ["des", "sds", "suf", "bck", "prj", "ssp", "tis"]

    if not algorithm in suffixes:
        return None

    if not seqType in suffixes[algorithm]:
        return None

    return suffixes[algorithm][seqType]
