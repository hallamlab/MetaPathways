"""Scans the contigs for tRNA genes """

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import sys
    import re
    from os import path, _exit

    from optparse import OptionParser, OptionGroup

    from metapathways import metapathways_utils as mputils
    from metapathways import general_utils as gutils
    from metapathways import sysutil as sysutils
    from metapathways import errorcodes as errormod
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)

PATHDELIM = sysutils.pathDelim()
errorcode = 7

help = sys.argv[0] + """ -i input -o output [algorithm dependent options]"""

def createParser():
    epilog = """This script is used for scanning for tRNA,  using tRNA-Scan 1.4,
              on the set of metagenomics sample sequences """
    epilog = re.sub(r"\s+", " ", epilog)
    parser = OptionParser(usage=help, epilog=epilog)

    # Input options
    parser.add_option(
        "--mode",
        dest="trna_mode",
        default='B',
        help="search mode for tRNAs",
    )

    parser.add_option(
        "-i",
        dest="trna_i",
        default=None,
        help="input fasta file"
    )

    parser.add_option(
        "-o",
        dest="trna_o",
        default=None,
        help="save final results in <file>"
    )

    parser.add_option(
        "-m",
        dest="trna_m",
        default=None,
        help="save statistics summary for run in <file>"
    )

    parser.add_option(
        "-j",
        dest="trna_j",
        default=None,
        help="save results in GFF3 file format of <file>"
    )

    parser.add_option(
        "-a",
        dest="trna_a",
        default=None,
        help="save predicted tRNA sequences in FASTA file format of <file>",
    )

    parser.add_option(
        "-p",
        dest="trna_p",
        default=None,
        help="use <label> prefix for all default output file names",
    )

    parser.add_option(
        "-q",
        dest="trna_q",
        default=True,
        help="quiet mode (credits & run option selections suppressed)",
    )

    parser.add_option(
        "-t",
        dest="trna_t",
        default=str(1),
        help="number of threads used for running infernal (default 1)",
    )

    parser.add_option(
        "--executable",
        dest="trna_executable",
        default=None,
        help="The tRNA-SCAN SE executable",
    )
    return parser


def main(argv, errorlogger=None, runcommand=None, runstatslogger=None):
    parser = createParser()
    options, args = parser.parse_args(argv)
    return _execute_tRNA_Scan(options)


def _execute_tRNA_Scan(options):
    global errorcode
    args = []

    if options.trna_executable:
        args.append(options.trna_executable)

    if options.trna_mode:
        args += ['-' + options.trna_mode]

    if options.trna_o:
        args += ["-o", options.trna_o]

    if options.trna_m:
        args += ["-m", options.trna_m]

    if options.trna_j:
        args += ["-j", options.trna_j]

    if options.trna_a:
        args += ["-a", options.trna_a]

    if options.trna_p:
        args += ["-p", options.trna_p]

    if options.trna_q:
        args += ["--quiet"]

    if options.trna_t:
        args += ["--thread", options.trna_t]
       
    if options.trna_i:
        args += [options.trna_i]
    
    result = sysutils.getstatusoutput(" ".join(args))

    if result[0] != 0:
        errormod.insert_error(errorcode)
    return result


def MetaPathways_tRNA_scan(
    argv, extra_command=None, errorlogger=None, runstatslogger=None):
    global errorcode
    if errorlogger != None:
        errorlogger.write("#STEP\ttRNA_SCAN\n")

    try:
        main(
            argv,
            errorlogger=errorlogger,
            runcommand=extra_command,
            runstatslogger=runstatslogger,
        )
    except:
        errormod.insert_error(errorcode)
        return (1, traceback.print_exc(10))

    return (0, "")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        result = main(sys.argv[1:])
