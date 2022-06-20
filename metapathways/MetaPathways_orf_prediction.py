"""This script run the orf prediction """

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import sys
    import re

    from os import path, _exit, rename, remove
    from optparse import OptionParser, OptionGroup

    from metapathways import metapathways_utils as mputils
    from metapathways import general_utils as gutils
    from metapathways import sysutil as sysutils
    from metapathways import errorcodes as errormod
    from metapathways import fastareader as fastreadrmod
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)

PATHDELIM = sysutils.pathDelim()
errorcode = 2

usage = sys.argv[0] + """ --algorithm <algorithm> [algorithm dependent options]"""

def createParser():

    epilog = """
    The preprocessed nucleotide sequences (contigs) are used as
    input to a gene prediction algorithm, currently prodigal, to detect the
    gene coding regions.  The output of the prodigal run is a set of
    untranslated ORFs and the same ORFs translated (into amino acid sequences).
    The resulting files are available in the 'orf_prediction' folder.
    The translation is done based on the translation table id provided by
    the user, by default it 11
    """

    epilog = re.sub(r"\s+", " ", epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    # Input options

    parser.add_option(
        "--algorithm",
        dest="algorithm",
        default="prodigal",
        choices=["prodigal"],
        help="default : prodigal ORF prediction algorithm [prodigal]",
    )

    parser.add_option(
        "--nthreads", dest="nthreads", default="1", help="number of threads default : 1"
    )

    prodigal_group = OptionGroup(parser, "Prodigal parameters")

    prodigal_group.add_option(
        "--prod_input",
        dest="prod_input",
        default=None,
        help="the input sequences  <inputfile>",
    )

    prodigal_group.add_option(
        "--prod_output", dest="prod_output", default=None, help="the output <outfile>"
    )

    prodigal_group.add_option(
        "--prod_p",
        dest="prod_p",
        default=None,
        help="Select procedure (single or meta).  Default is single",
    )

    prodigal_group.add_option(
        "--prod_f",
        dest="prod_f",
        default="gff",
        help="Select output format (gbk, gff, or sco).  default is gff",
    )

    prodigal_group.add_option(
        "--prod_g",
        dest="prod_g",
        default="11",
        help="Specify a translation table to use (default 11)",
    )

    prodigal_group.add_option(
        "--prod_m",
        dest="prod_m",
        action="store_true",
        default=False,
        help="Treat runs of n's as masked sequence and do not build genes across them",
    )

    prodigal_group.add_option(
        "--prod_exec", dest="prod_exec", default=None, help="prodigal executable"
    )

    prodigal_group.add_option(
        "--strand",
        dest="strand",
        default="both",
        choices=["both", "pos", "neg"],
        help="strands to use in case of transcriptomic sample",
    )

    parser.add_option_group(prodigal_group)
    return parser


def main(argv, errorlogger=None, runcommand=None, runstatslogger=None):
    parser = createParser()
    options, args = parser.parse_args(argv)

    if options.algorithm == "prodigal":
        _execute_prodigal(options)


def _execute_prodigal(options):
    args = []

    if options.prod_exec:
        args.append(options.prod_exec)
    if options.prod_m:
        args.append("-m")
    if options.prod_p:
        args += ["-p", options.prod_p]
    if options.prod_f:
        args += ["-f", options.prod_f]
    if options.prod_g:
        args += ["-g", options.prod_g]
    if options.prod_input:
        args += ["-i", options.prod_input]
    if options.prod_output:
        args += ["-o", options.prod_output + ".tmp"]

    result = sysutils.getstatusoutput(" ".join(args))
    rename(options.prod_output + ".tmp", options.prod_output)
    
    return result[0]


def MetaPathways_orf_prediction(
    argv, extra_command=None, errorlogger=None, runstatslogger=None):
    global errorcode
    if errorlogger != None:
        errorlogger.write("#STEP\tORF_PREDICTION\n")
    try:
        main(
            argv,
            errorlogger=errorlogger,
            runcommand=extra_command,
            runstatslogger=runstatslogger,
        )
    except:
        errormod.insert_error(errrocode)
        return (1, traceback.print_exc(10))

    return (0, "")

if __name__ == "__main__":
    if len(sys.argv) > 1: 
       main(sys.argv[1:])
