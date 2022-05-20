"""This script runs the functional search via BLAST or FAST like
homology search tools on a set of ORFs
"""
__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import sys
    import re

    from os import path, _exit, rename
    from optparse import OptionParser, OptionGroup
    from collections import namedtuple

    from metapathways import sysutil as sysutils
    from metapathways import general_utils as gutils
    from metapathways import metapathways_utils as mputils
    from metapathways import sysutil as sysutils
    from metapathways import errorcodes as errormod
    from metapathways import MetaPathways_parse_hmmer
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)

PATHDELIM = sysutils.pathDelim()

usage = sys.argv[0] + """ -i input -o output [algorithm dependent options]"""

def createParser():
    epilog = """This script is used for running a homology search algorithm such as BLAST or FAST
              on a set of query amino acid sequences against a target of  reference protein sequences.
              Currently it supports the BLASTP and FAST algorithm. Any other homology search algorithm
              can be added by first adding the new algorithm name in upper caseusing in to the
              choices parameter in the algorithm option of this script.
              The results are put in a tabular form in the folder blast_results, with individual files
               for each of the databases. The files are named as "<samplename>.<dbname>.<algorithm>out"
              In the case of large number of amino acid sequences, this step of the computation can be
               also done using multiple grids (to use batch processing system) """

    epilog = re.sub(r"\s+", " ", epilog)
    parser = OptionParser(usage=usage, epilog=epilog)

    # Input options
    parser.add_option(
        "--algorithm",
        dest="algorithm",
        default="BLAST",
        choices=["BLAST", "FAST", "HMMER"],
        help="the homology search algorithm",
    )

    blast_group = OptionGroup(parser, "BLAST parameters")

    blast_group.add_option(
        "--blast_query",
        dest="blast_query",
        default=None,
        help="Query amino acid sequences for BLASTP",
    )

    blast_group.add_option(
        "--blast_db",
        dest="blast_db",
        default=None,
        help="Target reference database sequenes for BLASTP",
    )

    blast_group.add_option(
        "--blast_out", dest="blast_out", default=None, help="BLAST output file"
    )

    blast_group.add_option(
        "--blast_outfmt",
        dest="blast_outfmt",
        default="6",
        help="BLASTP output format [default 6, tabular]",
    )

    blast_group.add_option(
        "--blast_evalue",
        dest="blast_evalue",
        default=None,
        help="The e-value cutoff for the BLASTP",
    )

    blast_group.add_option(
        "--num_threads",
        dest="num_threads",
        default="1",
        type="str",
        help="Number of BLAST threads",
    )

    blast_group.add_option(
        "--blast_max_target_seqs",
        dest="blast_max_target_seqs",
        default=None,
        help="Maximum number of target hits per query",
    )

    blast_group.add_option(
        "--blast_executable",
        dest="blast_executable",
        default=None,
        help="The BLASTP executable",
    )

    blast_group.add_option(
        "--num_hits",
        dest="num_hits",
        default="10",
        type="str",
        help="The BLASTP executable",
    )

    parser.add_option_group(blast_group)

    last_group = OptionGroup(parser, "FAST parameters")

    last_group.add_option(
        "--last_query",
        dest="last_query",
        default=None,
        help="Query amino acid sequences for FAST",
    )

    last_group.add_option(
        "--last_db",
        dest="last_db",
        default=None,
        help="Target reference database sequenes for FAST",
    )

    last_group.add_option(
        "--last_f",
        dest="last_f",
        default="0",
        help="FAST output format [default 0, tabular]",
    )

    last_group.add_option(
        "--last_o", dest="last_o", default=None, help="FAST output file"
    )

    last_group.add_option(
        "--last_executable",
        dest="last_executable",
        default=None,
        help="The FAST executable",
    )

    parser.add_option_group(last_group)


    return parser


def main(argv, errorlogger=None, runcommand=None, runstatslogger=None):
    parser = createParser()
    options, args = parser.parse_args(argv)

    if options.algorithm == "BLAST":
        (code, message) = _execute_BLAST(options, logger=errorlogger)

    elif options.algorithm == "FAST":
        (code, message) = _execute_FAST(options, logger=errorlogger)
    else:
        gutils.eprintf("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")
        if errorlogger:
            errorlogger.printf("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")
        # exit_process("ERROR\tUnrecognized algorithm name for FUNC_SEARCH\n")
        return -1

    if code != 0:
        a = "\nERROR\tCannot successfully execute the %s for FUNC_SEARCH\n" % (
            options.algorithm
        )
        b = "ERROR\t%s\n" % (message)
        c = "INFO\tDatabase you are searching against may not be formatted correctly (if it was formatted for an earlier version) \n"
        code = -1
        d = "INFO\tTry removing the files for that database in 'formatted' subfolder for MetaPathways to trigger reformatting \n"
        if options.algorithm == "BLAST":
            e = "INFO\tYou can remove as 'rm %s.*','\n" % (options.blast_db)
        if options.algorithm == "FAST":
            e = "INFO\tYou can remove as 'rm %s.*','\n" % (options.last_db)

        (code, message) = _execute_FAST(options, logger=errorlogger)
        f = "INFO\tIf removing the files did not work then format it manually (see manual)"
        outputStr = a + b + c + d + e + f

        gutils.eprintf(outputStr + "\n")

        if errorlogger:
            errorlogger.printf(outputStr + "\n")
        return code

    return 0


def _execute_FAST(options, logger=None):
    args = []

    if options.last_executable:
        args.append(options.last_executable)

    if options.last_f:
        args += ["-f", options.last_f]

    if options.last_o:
        args += ["-o", options.last_o + ".tmp"]

    if options.num_threads:
        args += ["-P", options.num_threads]

    args += [" -K", options.num_hits]

    if options.last_db:
        args += [options.last_db]

    if options.last_query:
        args += [options.last_query]

    result = None
    try:
        result = sysutils.getstatusoutput(" ".join(args))
        rename(options.last_o + ".tmp", options.last_o)
    except:
        message = "Could not run FAST correctly"
        if result and len(result) > 1:
            message = result[1]
        if logger:
            logger.printf("ERROR\t%s\n", message)
        return (1, message)

    return (result[0], result[1])


def _execute_BLAST(options, logger=None):
    args = []

    if options.blast_executable:
        args.append(options.blast_executable)

    if options.blast_max_target_seqs:
        args += ["-max_target_seqs", options.blast_max_target_seqs]

    if options.num_threads:
        args += ["-num_threads", options.num_threads]

    if options.blast_outfmt:
        args += ["-outfmt", options.blast_outfmt]

    if options.blast_db:
        args += ["-db", options.blast_db]

    if options.blast_query:
        args += ["-query", options.blast_query]

    if options.blast_evalue:
        args += ["-evalue", options.blast_evalue]

    if options.blast_out:
        args += ["-out", options.blast_out + ".tmp"]

    result = sysutils.getstatusoutput(" ".join(args))
    rename(options.blast_out + ".tmp", options.blast_out)
    return (result[0], result[1])

def MetaPathways_func_search(
    argv, extra_command=None, errorlogger=None, runstatslogger=None):

    if errorlogger != None:
        errorlogger.write("#STEP\tFUNC_SEARCH\n")
    try:
         main(
             argv,
             errorlogger = errorlogger,
             runcommand = extra_command,
             runstatslogger = runstatslogger,
            )
    except:
        errormod.insert_error(4)
        return (1, traceback.print_exc(10))

    return (0, "")


def define_hmm_domtbl_thresholds(stringency: str, hmm_cov: int, query_cov: int) -> namedtuple:
    thresholds_nt = namedtuple("thresholds", ["perc_aligned", "query_aligned",
                                              "min_acc", "max_e", "max_ie", "min_score",
                                              "profile_match"])

    for opt, value in {"hmm_coverage": hmm_cov, "query_coverage": query_cov}.items():
        if not 1 <= value <= 100:
            ts_logger.error("Option '{}' needs to be between 1 and 100 percent (currently {}).\n"
                                 .format(opt, value))
            sys.exit(3)

    # Parameterizing the hmmsearch output parsing:
    if stringency == "relaxed":
        domtbl_thresholds = thresholds_nt(perc_aligned=hmm_cov, query_aligned=query_cov,
                                          min_acc=0.7, max_e=1E-3, max_ie=1E-1, min_score=15, profile_match=False)
    elif stringency == "strict":
        domtbl_thresholds = thresholds_nt(perc_aligned=hmm_cov, query_aligned=query_cov,
                                          min_acc=0.7, max_e=1E-5, max_ie=1E-3, min_score=30, profile_match=False)
    else:
        ts_logger.error("Unknown HMM-parsing stringency option '{}'.\n".format(stringency))
        sys.exit(3)
    return domtbl_thresholds


def best_discrete_matches(matches: list) -> list:
    """
    Function for finding the best alignment in a list of HmmMatch() objects
    The best match is based off of the full sequence score

    :param matches: A list of HmmMatch() objects
    :return: List of the best HmmMatch's
    """
    # Code currently only permits multi-domains of the same gene
    dropped_annotations = list()
    len_sorted_matches = sorted(matches, key=lambda x: x.end - x.start)
    i = 0
    orf = len_sorted_matches[0].orf
    while i + 1 < len(len_sorted_matches):
        j = i + 1
        a_match = len_sorted_matches[i]  # type HmmMatch
        while j < len(len_sorted_matches):
            b_match = len_sorted_matches[j]  # type HmmMatch
            if a_match.target_hmm != b_match.target_hmm:
                if MetaPathways_parse_hmmer.detect_orientation(a_match.start, a_match.end,
                                                       b_match.start, b_match.end) != "satellite":
                    if a_match.full_score > b_match.full_score:
                        dropped_annotations.append(len_sorted_matches.pop(j))
                        j -= 1
                    else:
                        dropped_annotations.append(len_sorted_matches.pop(i))
                        j = len(len_sorted_matches)
                        i -= 1
            j += 1
        i += 1

    if len(len_sorted_matches) == 0:
        LOGGER.error("All alignments were discarded while deciding the best discrete HMM-match.\n")
        sys.exit(3)

    LOGGER.debug("HMM search annotations for " + orf +
                 ":\n\tRetained\t" +
                 ', '.join([match.target_hmm +
                            " (%d-%d)" % (match.start, match.end) for match in len_sorted_matches]) +
                 "\n\tDropped\t\t" +
                 ', '.join([match.target_hmm +
                            " (%d-%d)" % (match.start, match.end) for match in dropped_annotations]) + "\n")
    return len_sorted_matches


def parse_domain_tables(thresholds, hmm_domtbl_files: dict) -> dict:
    """
    Parses HMMER domain tables using predetermined thresholds

    :param thresholds: A namedtuple instance: namedtuple("thresholds", "max_e max_ie min_acc min_score perc_aligned")
    :param hmm_domtbl_files: A list of domain table files written by hmmsearch
    :return: Dictionary of HmmMatch objects indexed by their reference package and/or HMM name
    """
    # Check if the HMM filtering thresholds have been set
    LOGGER.info("Parsing HMMER domain tables for high-quality matches... ")

    search_stats = MetaPathways_parse_hmmer.HmmSearchStats()
    hmm_matches = dict()
    orf_gene_map = dict()
    optional_matches = list()

    # TODO: Capture multimatches across multiple domain table files
    for r_q, domtbl_file in hmm_domtbl_files.items():
        _prefix, reference = r_q
        domain_table = MetaPathways_parse_hmmer.DomainTableParser(domtbl_file)
        domain_table.read_domtbl_lines()
        distinct_hits = MetaPathways_parse_hmmer.format_split_alignments(domain_table, search_stats)
        purified_hits = MetaPathways_parse_hmmer.filter_poor_hits(thresholds, distinct_hits, search_stats)
        complete_hits = MetaPathways_parse_hmmer.filter_incomplete_hits(thresholds, purified_hits, search_stats)
        MetaPathways_parse_hmmer.renumber_multi_matches(complete_hits)

        for match in complete_hits:
            match.genome = reference
            if match.orf not in orf_gene_map:
                orf_gene_map[match.orf] = dict()
            try:
                orf_gene_map[match.orf][match.target_hmm].append(match)
            except KeyError:
                orf_gene_map[match.orf][match.target_hmm] = [match]
            if match.target_hmm not in hmm_matches.keys():
                hmm_matches[match.target_hmm] = list()
    search_stats.num_dropped()
    for orf in orf_gene_map:
        if len(orf_gene_map[orf]) == 1:
            for target_hmm in orf_gene_map[orf]:
                for match in orf_gene_map[orf][target_hmm]:
                    hmm_matches[target_hmm].append(match)
                    search_stats.seqs_identified += 1
        else:
            search_stats.multi_alignments += 1
            # Remove all the overlapping domains - there can only be one highlander
            for target_hmm in orf_gene_map[orf]:
                optional_matches += orf_gene_map[orf][target_hmm]
            retained = 0
            for discrete_match in best_discrete_matches(optional_matches):
                hmm_matches[discrete_match.target_hmm].append(discrete_match)
                retained += 1
            search_stats.dropped += (len(optional_matches) - retained)
            search_stats.seqs_identified += retained
            optional_matches.clear()

    LOGGER.info("done.\n")

    alignment_stat_string = search_stats.summarize()

    if search_stats.seqs_identified == 0 and search_stats.dropped == 0:
        LOGGER.warning("No alignments found.\n")
        sys.exit(0)
    if search_stats.seqs_identified == 0 and search_stats.dropped > 0:
        LOGGER.warning("No alignments (" + str(search_stats.seqs_identified) + '/' + str(search_stats.dropped) +
                       ") met the quality cut-offs!\n")
        alignment_stat_string += "\tPoor quality alignments:\t" + str(search_stats.bad) + "\n"
        alignment_stat_string += "\tShort alignments:\t" + str(search_stats.short) + "\n"
        LOGGER.debug(alignment_stat_string)
        sys.exit(0)

    alignment_stat_string += "\n\tNumber of markers identified:\n"
    for marker in sorted(hmm_matches):
        alignment_stat_string += "\t\t" + marker + "\t" + str(len(hmm_matches[marker])) + "\n"

    LOGGER.debug(alignment_stat_string)
    return hmm_matches


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
