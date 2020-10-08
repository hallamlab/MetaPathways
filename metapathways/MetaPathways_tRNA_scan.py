"""Scans the contigs for tRNA genes """

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
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
        "-o",
        dest="trna_o",
        default=None,
        help="Output from the tRNA-Scan 1.4 into <outfile>",
    )

    parser.add_option(
        "-i", dest="trna_i", default=None, help="reads the sequences from <input> file"
    )

    parser.add_option(
        "-T", dest="trna_T", default="6", help="reads the Tsignal from <TPCsignal>"
    )

    parser.add_option(
        "-D", dest="trna_D", default=None, help="reads the Dsignal from <Dsignal>"
    )

    parser.add_option(
        "-F",
        dest="trna_F",
        default=None,
        help="write predicted tRNA genes in fasta format<outfile>",
    )

    parser.add_option(
        "--executable",
        dest="trna_executable",
        default=None,
        help="The tRNA-SCAN 1.4 executable",
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

    if options.trna_i:
        args += ["-i", options.trna_i]

    if options.trna_o:
        args += ["-o", options.trna_o]

    if options.trna_D:
        args += ["-D", options.trna_D]

    if options.trna_T:
        args += ["-T", options.trna_T]

    if options.trna_F:
        args += ["-F", options.trna_F]

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
