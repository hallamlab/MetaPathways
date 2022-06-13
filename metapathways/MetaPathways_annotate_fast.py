#!/usr/bin/python

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import re
    import glob

    from os import makedirs, sys, remove, rename
    from sys import path
    from optparse import OptionParser, OptionGroup

    from metapathways import general_utils as gutils
    from metapathways import metapathways_utils as mputils
    from metapathways import sysutil as sysutils
    from metapathways import errorcodes as errormod
    from metapathways import gff_compare as gffc
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)


usage = (
    sys.argv[0]
    + """ -d dbname1 -b parsed_blastout_for_database1 -w weight_for_database1 
         [-d dbname2 -b parsed_blastout_for_database2 -w weight_for_database2 ]
         [ --rRNA_16S  16SrRNA-stats-table ] [ --tRNA tRNA-stats-table ]
         [ --compact_output ]
         """
)

errrocode = 8
def createParser():
    epilog = """Based on the parsed homology search results against reference protein databases,
               tRNA scan results and matches against the LSU and SSU rRNA sequence databases
               annotations for individual ORFs are created.
               The functional annotation is done by synthesizing  the annotations from the
               reference databases, after removing redundant texts and also making sure
               the annotations from different reference databases are in agreement with each other.
               The taxonomic annotation is done by applying the LCA (lowest common ancestor rule)
               on the hits from the RefSeq NR protein database."""

    epilog = re.sub(r"\s+", " ", epilog)
    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option(
        "-b",
        "--blastoutput",
        dest="input_blastout",
        action="append",
        default=[],
        help="blastout files in TSV format [at least 1 REQUIRED]",
    )
    parser.add_option(
        "-a",
        "--algorithm",
        dest="algorithm",
        default="BLAST",
        help="algorithm BLAST or LAST",
    )

    parser.add_option(
        "-m",
        "--contig-map",
        dest="contig_map_file",
        default=None,
        help="contig map file",
    )

    parser.add_option(
        "-d",
        "--dbasename",
        dest="database_name",
        action="append",
        default=[],
        help="the database names [at least 1 REQUIRED]",
    )

    parser.add_option(
        "-D",
        "--blastdir",
        dest="blastdir",
        default=None,
        help="the blast dir where all the BLAST outputs are located",
    )

    parser.add_option(
        "-s",
        "--samplename",
        dest="sample_name",
        default=None,
        help="the sample name"
    )

    parser.add_option(
        "-w",
        "--weight_for_database",
        dest="weight_db",
        action="append",
        default=[],
        type="float",
        help="the map file for the database  [at least 1 REQUIRED]",
    )

    parser.add_option(
        "--rRNA",
        dest="rRNA",
        action="append",
        default=[],
        help="the SSU rRNA stats file [OPTIONAL]",
    )

    parser.add_option(
        "--tRNA",
        dest="tRNA",
        action="append",
        default=[],
        help="the tRNA stats file [OPTIONAL]",
    )

    parser.add_option(
        "--rRNA_gff",
        dest="rRNA_gff",
        default=None,
        help="the SSU rRNA GFF file [OPTIONAL]",
    )

    parser.add_option(
        "--tRNA_gff",
        dest="tRNA_gff",
        default=None,
        help="the tRNA GFF file [OPTIONAL]",
    )

    cutoffs_group = OptionGroup(parser, "Cuttoff Related Options")

    cutoffs_group.add_option(
        "--min_score",
        dest="min_score",
        type="float",
        default=20,
        help="the minimum bit score cutoff [default = 20 ] ",
    )

    cutoffs_group.add_option(
        "--max_evalue",
        dest="max_evalue",
        type="float",
        default=1e-6,
        help="the maximum E-value cutoff [ default = 1e-6 ] ",
    )
    cutoffs_group.add_option(
        "--min_length",
        dest="min_length",
        type="float",
        default=30,
        help="the minimum length of query cutoff [default = 30 ] ",
    )
    cutoffs_group.add_option(
        "--max_length",
        dest="max_length",
        type="float",
        default=10000,
        help="the maximum length of query cutoff [default = 10000 ] ",
    )

    cutoffs_group.add_option(
        "--min_identity",
        dest="min_identity",
        type="float",
        default=20,
        help="the minimum identity of query cutoff [default 30 ] ",
    )
    cutoffs_group.add_option(
        "--max_identity",
        dest="max_identity",
        type="float",
        default=100,
        help="the maximum identity of query cutoff [default = 100 ] ",
    )

    cutoffs_group.add_option(
        "--limit",
        dest="limit",
        type="float",
        default=5,
        help="max number of hits per query cutoff [default = 5 ] ",
    )

    cutoffs_group.add_option(
        "--min_bsr",
        dest="min_bsr",
        type="float",
        default=0.00,
        help="minimum BIT SCORE RATIO [default = 0.00 ] ",
    )
    parser.add_option_group(cutoffs_group)

    output_options_group = OptionGroup(parser, "Output Options")
    output_options_group.add_option(
        "--tax",
        dest="taxonomy",
        action="store_true",
        default=False,
        help="add the taxonomy info [useful for refseq] ",
    )
    parser.add_option_group(output_options_group)

    parser.add_option(
        "-o",
        "--output_gff",
        dest="output_gff",
        help="the output gff file [REQUIRED]"
    )

    parser.add_option(
        "--output-comparative-annotation",
        dest="output_comparative_annotation",
        help="the comparative output table [REQUIRED]",
    )

    parser.add_option(
        "--compact_output",
        dest="compact_output",
        action="store_true",
        default=False,
        help="compact output [OPTIONAL]",
    )

    parser.add_option(
        "--input_gff",
        dest="input_gff",
        metavar="INPUT",
        help="Unannotated gff file [REQUIRED]",
    )
    return parser


def check_arguments(opts, args):

    return True

    if len(opts.input_blastout) == 0:
        gutils.eprintf("There should be at least one blastoutput file\n")
        return False

    if len(opts.database_name) == 0:
        gutils.eprintf("There should be at least one database name\n")
        return False

    if len(opts.weight_db) == 0:
        eprint("There should be at least one weight\n")
        return False

    if len(opts.input_blastout) != len(opts.database_name) or len(
        opts.input_blastout
    ) != len(opts.weight_db):
        eprint(
            "The num of database names, blastoutputs and database map file should be equal\n"
        )
        return False

    if opts.output_gff == None:
        gutils.eprintf("Must specify the output gff file\n")
        return False

    if opts.output_comparative_annotation == None:
        gutils.eprintf("Must specify the output tables for comparative annotation\n")
        return False

    if opts.input_gff == None:
        gutils.eprintf("Must specify the input gff file\n")
        return False

    return True


def insert_attribute(attributes, attribStr):
    rawfields = re.split("=", attribStr)
    if len(rawfields) == 2:
        if rawfields[0].strip().lower() == "id":
            orfid = mputils.ShortenORFId(rawfields[1].strip())
            attributes[rawfields[0].strip().lower()] = orfid
            attributes['orf_id'] = rawfields[1].strip()
        else:
            attributes[rawfields[0].strip().lower()] = rawfields[1].strip()


def split_attributes(str, attributes):
    rawattributes = re.split(";", str)
    for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

    return attributes


def insert_orf_into_dict(line, contig_dict, shortenorfid=False):
    rawfields = re.split("\t", line)

    fields = []
    for field in rawfields:
        fields.append(field.strip())

    if len(fields) != 9:
        return

    attributes = {}

    seqname = fields[0]
    try:
        if shortenorfid:
            seqname = mputils.ShortenContigId(fields[0])
        else:
            seqname = fields[0]
    except:
        seqname = fields[0]

    attributes["seqname"] = seqname  # this is a bit of a  duplication
    attributes["source"] = fields[1]
    attributes["feature"] = fields[2]
    attributes["start"] = int(fields[3])
    attributes["end"] = int(fields[4])

    try:
        attributes["score"] = float(fields[5])
    except:
        attributes["score"] = fields[5]

    attributes["strand"] = fields[6]
    attributes["frame"] = fields[7]

    split_attributes(fields[8], attributes)

    if not seqname in contig_dict:
        contig_dict[seqname] = []

    contig_dict[seqname].append(attributes)


class GffFileParser(object):
    def __init__(self, gff_filename, shortenorfid=False):
        self.Size = 10000
        self.i = 0
        self.orf_dictionary = {}
        self.gff_beg_pattern = re.compile("^#")
        self.lines = []
        self.size = 0
        self.shortenorfid = shortenorfid
        try:
            self.gff_file = open(gff_filename, "r")
        except AttributeError:
            gutils.eprintf("Cannot read the map file for database : %s\n", dbname)
            exit_process()

    def __iter__(self):
        return self

    def refillBuffer(self):
        self.orf_dictionary = {}
        i = 0
        while i < self.Size:
            line = self.gff_file.readline()
            if not line:
                break
            if self.gff_beg_pattern.search(line):
                continue
            insert_orf_into_dict(line, self.orf_dictionary, self.shortenorfid)
            # print self.orf_dictionary
            i += 1

        self.orfs = list(self.orf_dictionary.keys())
        self.size = len(self.orfs)
        self.i = 0

    def __next__(self):
        if self.i == self.size:
            self.refillBuffer()

        if self.size == 0:
            self.gff_file.close()
            raise StopIteration()

        # print self.i
        if self.i < self.size:
            self.i = self.i + 1
            tmp = self.i - 1
            return self.orfs[tmp]


def process_gff_file(gff_file_name, orf_dictionary):
    try:
        gfffile = open(gff_file_name, "r")
    except IOError:
        gutils.eprintf("Cannot read file %s!\n", gff_file_name)

    gff_lines = gfffile.readlines()
    gff_beg_pattern = re.compile("^#")
    gfffile.close()

    count = 0
    for line in gff_lines:
        line = line.strip()
        if gff_beg_pattern.search(line):
            continue
        insert_orf_into_dict(line, orf_dictionary)
        count += 1
        # if count %10000 == 0:
        #   print count


def create_dictionary(databasemapfile, annot_map):
    seq_beg_pattern = re.compile(">")

    dbmapfile = open(databasemapfile, "r")
    lines = dbmapfile.readlines()
    dbmapfile.close()
    for line in lines:
        if seq_beg_pattern.search(line):
            words = line.rstrip().split()
            name = words[0].replace(">", "", 1)

            words.pop(0)
            annotation = " ".join(words)
            annot_map[name] = annotation


def write_annotation_for_orf(
    outputgff_file,
    candidatedbname,
    dbname_weight,
    results_dictionary,
    orf_dictionary,
    contig,
    candidate_orf_pos,
    orfid,
    sample_name,
    compact_output=True,
    ):
    global errorcode
    try:
        fields = ["source", "feature", "start", "end", "score", "strand", "frame"]

        output_line = orf_dictionary[contig][candidate_orf_pos]["seqname"]

        # if compact_output:
        # output_line = ShortenContigId(output_line)

        for field in fields:
            output_line += "\t" + str(orf_dictionary[contig][candidate_orf_pos][field])

        # if compact_output:
        try:
            if compact_output:
                attributes = "ID=" + mputils.ShortenORFId(
                    orf_dictionary[contig][candidate_orf_pos]["id"]
                )
            else:
                attributes = (
                    "ID="
                    + sample_name
                    + "_"
                    + mputils.ShortenORFId(orf_dictionary[contig][candidate_orf_pos]["id"])
                )

            attributes += (
                ";"
                + "locus_tag="
                + mputils.ShortenORFId(orf_dictionary[contig][candidate_orf_pos]["locus_tag"])
            )
        except:
            attributes = "ID=" + orf_dictionary[contig][candidate_orf_pos]["id"]
            attributes += (
                ";"
                + "locus_tag="
                + orf_dictionary[contig][candidate_orf_pos]["locus_tag"]
            )

        attributes += (
            ";"
            + "contig_length="
            + orf_dictionary[contig][candidate_orf_pos]["contig_length"]
        )
        attributes += (
            ";"
            + "orf_length="
            + orf_dictionary[contig][candidate_orf_pos]["orf_length"]
        )
        if "partial" in orf_dictionary[contig][candidate_orf_pos].keys():
            attributes += (
                ";" + "partial=" + orf_dictionary[contig][candidate_orf_pos]["partial"]
            )
        attributes += ";" + "sourcedb=" + candidatedbname

        if candidatedbname in results_dictionary:
            attributes += ";" + "target=" + results_dictionary[candidatedbname][orfid]["target"]
            attributes += (
                ";"
                + "annotvalue="
                + str(results_dictionary[candidatedbname][orfid]["value"])
            )
            attributes += (
                ";" + "ec=" + str(results_dictionary[candidatedbname][orfid]["ec"])
            )
            attributes += (
                ";" + "product=" + results_dictionary[candidatedbname][orfid]["product"]
            )
        else:
            attributes += ";" + "target=" + str("")
            attributes += ";" + "annotvalue=" + str("0")
            attributes += ";" + "ec=" + str("")
            attributes += ";" + "product=" + "hypothetical protein"

        output_line += "\t" + attributes

        # print any db name for now
        # if candidatedbname in results_dictionary:
        gutils.fprintf(outputgff_file, "%s\n", output_line)
    except:
        gutils.eprintf("ERROR : Failure to annotate in contig %s\n", contig)
        # print orf_dictionary[contig]
        errormod.insert_error(errorcode)
        exit_process()


def write_16S_tRNA_gene_info(contig_id, f_rec, outputgff_file, tag):
    output_line = str(f_rec["seqname"])
    length = abs(f_rec["end"] - f_rec["start"])
    if tag == "_rRNA":
        output_line += "\t" + str(f_rec["source"])
        output_line += "\t" + str(f_rec["feature"])
        output_line += "\t" + str(f_rec["start"])
        output_line += "\t" + str(f_rec["end"])
        output_line += "\t" + str(f_rec["score"])
        output_line += "\t" + str(f_rec["strand"])
        output_line += "\t" + str(f_rec["frame"])
        attributes = "ID=" + str(f_rec["name"])
        attributes += ";" + "locus_tag=" + str(f_rec["seqname"]) + tag
        attributes += ";" + "orf_length=" + str(length)
        #attributes += ";" + "contig_length=" + str(orf_rec["contig_length"])
        attributes += ";" + "product=" + f_rec["product"]
        #attributes += ";" + "target=" + r_dictionary[orf_id]["taxonomy"]
        #attributes += ";" + "ec="
        output_line += "\t" + attributes

    elif tag == "_tRNA":
        output_line += "\t" + str(f_rec["source"])
        output_line += "\t" + str(f_rec["feature"])
        output_line += "\t" + str(f_rec["start"])
        output_line += "\t" + str(f_rec["end"])
        output_line += "\t" + str(f_rec["score"])
        output_line += "\t" + str(f_rec["strand"])
        output_line += "\t" + str(f_rec["frame"])
        attributes = "ID=" + str(f_rec["name"])
        attributes += ";" + "locus_tag=" + str(f_rec["seqname"]) + tag
        attributes += ";" + "orf_length="  + str(length)
        #attributes += ";" + "contig_length=" + str(r_dictionary[orf_id]["contig_length"])
        #attributes += ";" + "product=" + f_rec["product"]
        #attributes += ";" + "target=" + f_rec["product"]
        #attributes += ";" + "ec="
        output_line += "\t" + attributes
    gutils.fprintf(outputgff_file, "%s\n", output_line)


def process_rRNA_16S_stats(rRNA_16S_file, rRNA_16S_dictionary, shortenorfid=False):
    counter_rRNA = {}
    if not gutils.doesFileExist(rRNA_16S_file):
        return
    try:
        taxonomy_file = open(rRNA_16S_file, "r")
    except IOError:
        gutils.eprintf("Cannot read file %s!\n", rRNA_16S_file)
        exit_process()

    tax_lines = taxonomy_file.readlines()
    similarity_pattern = re.compile("similarity")
    evalue_pattern = re.compile("evalue")
    bitscore_pattern = re.compile("bitscore")
    taxonomy_pattern = re.compile("taxonomy")
    headerScanned = False
    for line in tax_lines:
        if headerScanned == False:
            if (
                similarity_pattern.search(line)
                and evalue_pattern.search(line)
                and bitscore_pattern.search(line)
                and taxonomy_pattern.search(line)
            ):
                headerScanned = True
            continue
        fields = [x.strip() for x in line.split("\t")]
        if len(fields) >= 6:

            if shortenorfid:
                name = get_sequence_number(fields[0])
            else:
                name = fields[0]

            if not name in counter_rRNA:
                counter_rRNA[name] = 0

            _name = name #+ "_" + str(counter_rRNA[name])
            counter_rRNA[name] = counter_rRNA[name] + 1

            if fields[1] != "-":
                rRNA_16S_dictionary[_name] = {
                                              'sequence': fields[0],
                                              'start': fields[1],
                                              'end': fields[2],
                                              'similarity': fields[3],
                                              'evalue': fields[4],
                                              'bitscore': fields[5],
                                              'taxonomy': fields[6]
                                              }

    taxonomy_file.close()


def get_sequence_number(line):
    seqnamePATT = re.compile(r"[\S]+_(\d+)$")
    result = seqnamePATT.search(line.strip())
    if result:
        return result.group(1)
    return line


def process_tRNA_stats(tRNA_stats_file, tRNA_dictionary, shortenorfid=False):
    counter_tRNA = {}
    try:
        tRNA_file = open(tRNA_stats_file, "r")
    except IOError:
        gutils.eprintf("Cannot read file %s!\n", tRNA_stats_file)
        exit_process()
    tRNA_lines = tRNA_file.readlines()

    sequence_name_pattern = re.compile("#seq_name", re.I)
    number_pattern = re.compile("number", re.I)

    headerScanned = False
    for line in tRNA_lines:
        if number_pattern.search(line):
            continue
        if headerScanned == False:
            if sequence_name_pattern.search(line):
                headerScanned = True
            continue
        fields = [x.strip() for x in line.split("\t")]
        if len(fields) >= 6:

            if shortenorfid:
                name = get_sequence_number(fields[0])
            else:
                name = fields[0]
            if not name in counter_tRNA:
                counter_tRNA[name] = 0

            _name = name #  + "_" + str(counter_tRNA[name])
            counter_tRNA[name] = counter_tRNA[name] + 1

            tRNA_dictionary[_name] = [fields[3], fields[4], fields[5], fields[1]]


# this adds the features and attributes to  be added to the gff file format for the tRNA genes
def add_tRNA_genes(tRNA_dictionary, tRNA_gff_dictionary, contig_lengths):
    for tRNA in tRNA_dictionary:
        start = int(tRNA_dictionary[tRNA][0])
        end = int(tRNA_dictionary[tRNA][1])
        if start > end:
            start = int(tRNA_dictionary[tRNA][1])
            end = int(tRNA_dictionary[tRNA][0])

        try:
            orf_length = end - start
        except:
            orf_length = 0

        contig_name = re.sub(r"_\d+$", "", tRNA)

        if contig_name in contig_lengths:
            contig_length = contig_lengths[contig_name]
        else:
            contig_length = 0

        if start > end or contig_length < end:
            gutils.eprintf(
                "trna {}   {}  {}  {} {}\n".format(
                    tRNA, start, end, end - start, contig_length
                )
            )
            end = contig_length
            gutils.eprintf(
                "trna {}   {}  {}  {} {}\n".format(
                    tRNA, start, end, end - start, contig_length
                )
            )

        dict = {
            "id": mputils.ContigID(tRNA),
            "seqname": tRNA,
            "start": start,
            "end": end,
            "strand": tRNA_dictionary[tRNA][2],
            "score": " ",
            "orf_length": str(orf_length),
            "contig_length": str(contig_length),
            "feature": "tRNA",
            "source": "trnaScan-1.4",
            "frame": 0,
            "product": "tRNA-" + tRNA_dictionary[tRNA][3],
            "ec": "",
        }
        tRNA_gff_dictionary[tRNA] = dict.copy()


# this adds the features and attributes to  be added to the gff file format for the 16S rRNA genes
def add_16S_genes(rRNA_16S_dictionary, rRNA_dictionary, contig_lengths):
    for rRNA in rRNA_16S_dictionary:
        start = rRNA_16S_dictionary[rRNA]['start']
        end = rRNA_16S_dictionary[rRNA]['end']
        if start > end:
            strand = '-'
        else:
            strand = '+'
        score = rRNA_16S_dictionary[rRNA]['bitscore']
        taxonomy = rRNA_16S_dictionary[rRNA]['taxonomy']
        similarity = rRNA_16S_dictionary[rRNA]['similarity']
        evalue = rRNA_16S_dictionary[rRNA]['evalue']


        dict = {
            "id": mputils.ContigID(rRNA),
            "seqname": rRNA,
            "start": str(start),
            "end": str(end),
            "strand": strand,
            "score": str(score),
            "feature": "rRNA",
            "taxonomy": str(taxonomy),
            "similarity": str(similarity),
            "evalue": str(evalue),
            "frame": 0,
            "ec": ""
            }
        rRNA_dictionary[rRNA] = dict.copy()


def create_annotation(
    dbname_weight,
    results_dictionary,
    input_gff,
    rRNA_16S_stats_files,
    rRNA_gff,
    tRNA_stats_files,
    tRNA_gff,
    output_gff,
    output_comparative_annotation,
    contig_lengths,
    sample_name,
    compact_output=False,
):
    orf_dictionary = {}
    #    process_gff_file(input_gff, orf_dictionary)
    gffreader = GffFileParser(input_gff)

    output_gff_tmp = output_gff + ".tmp"
    outputgff_file = open(output_gff_tmp, "w")

    output_comp_annot_file1 = open(output_comparative_annotation + ".1.txt", "w")
    output_comp_annot_file2 = open(output_comparative_annotation + ".2.txt", "w")

    output_comp_annot_file1_Str = "orf_id\tref dbname\ttarget\tproduct\tvalue"
    gutils.fprintf(output_comp_annot_file1, "%s\n", output_comp_annot_file1_Str)

    output_comp_annot_file2_Str = "orf_id"
    dbnames = dbname_weight.keys()
    for dbname in dbnames:
        weight = dbname_weight[dbname]
        output_comp_annot_file2_Str += "\t{0}(target)\t{0}(product)\t{0}(value)".format(
            dbname
        )
    gutils.fprintf(output_comp_annot_file2, "%s\n", output_comp_annot_file2_Str)

    
    # Deal with the rRNA sequences if there is rRNA stats file
    if len(rRNA_16S_stats_files) > 0 and contig_lengths:
        rRNA_reader = GffFileParser(rRNA_gff)
        '''
        rRNA_16S_dictionary = {}
        for rRNA_16S_stats_file in rRNA_16S_stats_files:
            process_rRNA_16S_stats(rRNA_16S_stats_file, rRNA_16S_dictionary)
        '''
        rRNA_dictionary = {}
        for contig in rRNA_reader:
            rec = rRNA_reader.orf_dictionary[contig]
            rRNA_dictionary[contig] = rec
        #add_16S_genes(rRNA_16S_dictionary, rRNA_dictionary, contig_lengths)
        
    # now deal with the tRNA sequences  if there is tRNA stats file
    if len(tRNA_stats_files) > 0 and contig_lengths:
        tRNA_reader = GffFileParser(tRNA_gff)
        '''
        tRNA_dictionary = {}
        for tRNA_stats_file in tRNA_stats_files:
            process_tRNA_stats(tRNA_stats_file, tRNA_dictionary)

        tRNA_gff_dictionary = {}
        add_tRNA_genes(tRNA_dictionary, tRNA_gff_dictionary, contig_lengths)
        '''
        tRNA_dictionary = {}
        for contig in tRNA_reader:
            rec = tRNA_reader.orf_dictionary[contig]
            tRNA_dictionary[contig] = rec

    values = {}
    i = 0
    for contig in gffreader:
        count = 0
        for orf in gffreader.orf_dictionary[contig]:
            orf_id = orf['orf_id']
            for dbname in dbnames:
                values[dbname] = 0.0001
            success = False
            output_comp_annot_file1_Str = ""
            output_comp_annot_file2_Str = ""
            orf_id = orf["id"]
            best_anno_dict = {'DB': 'None', 'weighted_value': 0}
            # check the annotation of the orf by dbname
            for dbname in dbnames:
                weight = dbname_weight[dbname]
                if orf_id in results_dictionary[dbname]:
                    if float(values[dbname]) < float(results_dictionary[dbname][orf_id]["value"]):
                        values[dbname] = float(results_dictionary[dbname][orf_id]["value"])
                        success = True
                        candidate_orf_pos = count
                        wt_val = results_dictionary[dbname][orf_id]["value"] * float(weight)
                        if wt_val > best_anno_dict['weighted_value']:
                            best_anno_dict['weighted_value'] = wt_val
                            best_anno_dict['DB'] = dbname
                        if output_comp_annot_file1_Str:
                            output_comp_annot_file1_Str += (
                                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                                    "",
                                    dbname,
                                    results_dictionary[dbname][orf_id]["target"],
                                    results_dictionary[dbname][orf_id]["product"],
                                    str(wt_val),
                                )
                            )
                        else:
                            output_comp_annot_file1_Str += (
                                "{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                                    orf_id,
                                    dbname,
                                    results_dictionary[dbname][orf_id]["target"],
                                    results_dictionary[dbname][orf_id]["product"],
                                    str(wt_val),
                                )
                            )

                        if output_comp_annot_file2_Str:
                            output_comp_annot_file2_Str += "\t{0}\t{1}\t{2}".format(
                                results_dictionary[dbname][orf_id]["target"],
                                results_dictionary[dbname][orf_id]["product"],
                                str(wt_val),
                            )
                        else:
                            output_comp_annot_file2_Str += "{0}\t{1}\t{2}\t{3}".format(
                                orf_id,
                                results_dictionary[dbname][orf_id]["target"],
                                results_dictionary[dbname][orf_id]["product"],
                                str(wt_val),
                            )

                else:
                    if not output_comp_annot_file1_Str:
                        output_comp_annot_file1_Str += (
                            "{0}\t{1}\t{2}\t{3}\t{4}\n".format(orf_id, "", "", "", "")
                        )

                    if output_comp_annot_file2_Str:
                        output_comp_annot_file2_Str += "\t{0}\t{1}\t{2}".format(
                            "", "", ""
                        )
                    else:
                        output_comp_annot_file2_Str += "{0}\t{1}\t{2}\t{3}".format(
                            orf_id, "", "", "", ""
                        )

            # end of for all dbname

            if success:  # there was a database hit
                gutils.fprintf(output_comp_annot_file1, "%s\n", output_comp_annot_file1_Str)
                gutils.fprintf(output_comp_annot_file2, "%s\n", output_comp_annot_file2_Str)
                write_annotation_for_orf(
                    outputgff_file,
                    best_anno_dict['DB'],
                    dbname_weight[best_anno_dict['DB']],
                    results_dictionary,
                    gffreader.orf_dictionary,
                    contig,
                    candidate_orf_pos,
                    orf_id,
                    sample_name,
                    compact_output=compact_output,
                )
            else:  # if it was not a hit then it is a hypothetical protein
                # print gffreader.orf_dictionary
                write_annotation_for_orf(
                    outputgff_file,
                    "None",
                    "0",
                    results_dictionary,
                    gffreader.orf_dictionary,
                    contig,
                    count,
                    orf_id,
                    sample_name,
                    compact_output=compact_output,
                )

            count += 1  # move to the next orf
    output_comp_annot_file1.close()
    output_comp_annot_file2.close()

    # Add rRNA and tRNA records to output gff
    for contig in rRNA_dictionary:
        rrna_rec = rRNA_dictionary[contig][0]
        write_16S_tRNA_gene_info(contig, rrna_rec, outputgff_file, "_rRNA")
    for contig in tRNA_dictionary:
        trna_rec = tRNA_dictionary[contig][0]
        write_16S_tRNA_gene_info(contig, trna_rec, outputgff_file, "_tRNA")

    outputgff_file.close()
    rename(output_gff_tmp, output_gff)


def process_product(product, database, similarity_threshold=0.9):
    """Returns the best set of products from the list of (*database*,
    *product*) tuples *products*.

    Each product in the set is first trimmed down, removing database-specific
    information.

    The set is then determined by first sorting the products by length
    (ascending), and then, for each product, sequentially applying the longest
    common substring algorithm to determine the similarity between the product
    and already determined products. If this similarity is greater than the
    specified *similarity_threshold*, the longer of the two products is chosen
    to be a determined product.
    """

    processed_product = ""
    if 'MULTISPECIES: ' in product:
        product = product.replace('MULTISPECIES: ', '')
    try:
        if 'cazy' in database:
            processed_product = product
            comment = ''

        elif (('cog' in database) or ('refseq' in database)):
            if '[' in product: # helps with edge cases
                processed_product = product.rsplit(' [', 1)[0]
                comment = '[' + product.rsplit(' [', 1)[1]
            else:
                processed_product = product
                comment = ''

        elif 'metacyc' in database:
            processed_product = product.rsplit(' (', 1)[0]
            comment = '(' + product.rsplit(' (', 1)[1]

        elif 'sprot' in database:
            processed_product = product.rsplit(' OS=', 1)[0]
            comment = ' OS=' + product.rsplit(' OS=', 1)[1]

        elif 'uniref' in database:
            processed_product = product.rsplit(' n=', 1)[0]
            comment = ' n=' + product.rsplit(' n=', 1)[1]
    except:
        gutils.eprintf("\n\tERROR:\tDatabase:%s\n", database)
        gutils.eprintf("\n\t\t\tProduct:%s\n", product)
        gutils.eprintf("\n\t\t\tHas unknown format\n")

    words = [x.strip() for x in processed_product.split()]
    filtered_words = []
    underscore_pattern = re.compile("_")
    arrow_pattern = re.compile(">")
    for word in words:
        if not underscore_pattern.search(word) and not arrow_pattern.search(word):
            filtered_words.append(word)

    processed_product = remove_repeats(filtered_words)
    processed_product = re.sub(";", "", processed_product)

    return processed_product, comment


def remove_repeats(filtered_words):
    word_dict = {}
    newlist = []
    for word in filtered_words:
        if not word in word_dict:
            if not word in [
                "",
                "is",
                "have",
                "has",
                "will",
                "can",
                "should",
                "in",
                "at",
                "upon",
                "the",
                "a",
                "an",
                "on",
                "for",
                "of",
                "by",
                "with",
                "and",
                ">",
            ]:
                if not word == "#":  # for metacyc this is important
                    word_dict[word] = 1

                newlist.append(word)
    return " ".join(newlist)


class BlastOutputTsvParser(object):
    def __init__(self, dbname, blastoutput, shortenorfid=False):
        self.dbname = dbname
        self.blastoutput = blastoutput
        self.i = 1
        self.data = {}
        self.fieldmap = {}
        self.shortenorfid = shortenorfid
        self.seq_beg_pattern = re.compile("#")

        try:
            self.blastoutputfile = open(blastoutput, "r")
            self.lines = self.blastoutputfile.readlines()
            self.blastoutputfile.close()
            self.size = len(self.lines)
            if not self.seq_beg_pattern.search(self.lines[0]):
                exit_process(
                    'First line must have field header names and begin with "#"'
                )
            header = self.lines[0].replace("#", "", 1)
            fields = [x.strip() for x in header.rstrip().split("\t")]
            k = 0
            for x in fields:
                self.fieldmap[x] = k
                k += 1
            gutils.eprintf("\n\tINFO:\tProcessing database : %s\n", dbname)

        except AttributeError:
            gutils.eprintf("\n\tERROR:\tCannot read the map file for database :%s\n", dbname)
            exit_process()

    def __iter__(self):
        return self

    count = 0

    def __next__(self):
        if self.i < self.size:

            try:
                fields = [x.strip() for x in self.lines[self.i].split("\t")]
                # print self.fieldmap['ec'], fields, self.i,  self.blastoutput
                if self.shortenorfid:
                    self.data["query"] = mputils.ShortenORFId(fields[self.fieldmap["query"]])
                else:
                    self.data["query"] = fields[self.fieldmap["query"]]
                self.data["target"] = fields[self.fieldmap["target"]]
                self.data["q_length"] = int(fields[self.fieldmap["q_length"]])
                self.data["bitscore"] = float(fields[self.fieldmap["bitscore"]])
                self.data["bsr"] = float(fields[self.fieldmap["bsr"]])
                self.data["expect"] = float(fields[self.fieldmap["expect"]])
                self.data["identity"] = float(fields[self.fieldmap["identity"]])
                self.data["ec"] = fields[self.fieldmap["ec"]]
                self.data["product"] = fields[self.fieldmap["product"]]

                self.i = self.i + 1
                return self.data
            except:
                sys.exit(0)
                return None
        else:
            raise StopIteration()


def isWithinCutoffs(data, cutoffs):
    if data["q_length"] < cutoffs.min_length:
        return False

    if data["bitscore"] < cutoffs.min_score:
        return False

    if data["expect"] > cutoffs.max_evalue:
        return False

    if data["identity"] < cutoffs.min_identity:
        return False

    if data["bsr"] < cutoffs.min_bsr:
        return False

    return True


def word_information(string_of_words):
    words = [x.strip() for x in string_of_words.split()]

    information = 0
    wordlist = {}
    underscore_pattern = re.compile("_")
    for word in words:
        if not word in [
            "",
            "is",
            "have",
            "has",
            "will",
            "can",
            "should",
            "in",
            "at",
            "upon",
            "the",
            "a",
            "an",
            "on",
            "for",
            "of",
            "by",
            "with",
            "and",
            ">",
            "predicted",
            "protein",
            "conserved",
        ]:
            if not underscore_pattern.search(word):
                wordlist[word] = 1

    # print string_of_words
    # print wordlist
    # print len(wordlist)
    return len(wordlist)


def compute_annotation_value(data):
    score = 0
    if len(data["ec"]) > 0:
        score += 10

    if not re.search(r"hypothetical protein", data["product"]):
        score += word_information(data["product"])

    return score


# compute the refscores
def process_parsed_blastoutput(
    dbname, weight, blastoutput, cutoffs, annotation_results
):
    blastparser = BlastOutputTsvParser(dbname, blastoutput, shortenorfid=False)

    fields = ["target", "q_length", "bitscore", "bsr", "expect", "aln_length", "identity", "ec"]
    if cutoffs.taxonomy:
        fields.append("taxonomy")
    fields.append("product")

    annotation = {}
    for data in blastparser:
        if isWithinCutoffs(data, cutoffs):
            annotation["target"] = data["target"]
            annotation["bsr"] = data["bsr"]
            annotation["ec"] = data["ec"]
            pprod, comm = process_product(data["product"], dbname)
            annotation["product"] = pprod #mputils.strip_taxonomy()
            annotation["comment"] = comm
            annotation["value"] = compute_annotation_value(annotation) * weight
            annotation["bitscore"] = data["bitscore"]

            if not data["query"] in annotation_results:
                annotation_results[data["query"]] = {"value": 0, "bitscore": 0.0}
            if annotation_results[data["query"]]["bitscore"] <= annotation["bitscore"]:
                if annotation_results[data["query"]]["value"] < annotation["value"]:
                    annotation_results[data["query"]] = annotation.copy()

    count = len(annotation_results.keys())
    return count


def read_contig_lengths(contig_map_file, contig_lengths):
    try:
        mapfile = open(contig_map_file, "r")
    except IOError:
        errormod.insert_error(errorcode)
        return

    mapfile_lines = mapfile.readlines()
    mapfile.close()

    for line in mapfile_lines:
        line = line.strip()
        fields = [x.strip() for x in line.split("\t")]
        if len(fields) != 3:
            contig_lengths = {}
            return
        contig_lengths[fields[0]] = int(fields[2])


def getBlastFileNames(opts):

    database_names = []
    parsed_blastouts = []
    weight_dbs = []

    dbnamePATT = re.compile(
        r""
        + opts.blastdir
        + "*"
        + opts.sample_name
        + "*[.](.*)[.]"
        + opts.algorithm.upper()
        + "out.parsed.txt"
    )

    blastOutNames = glob.glob(
        opts.blastdir + "*" + opts.algorithm.upper() + "out.parsed.txt"
    )
    for blastoutname in blastOutNames:
        result = dbnamePATT.search(blastoutname)
        if result:
            dbname = result.group(1)
            database_names.append(dbname)
            parsed_blastouts.append(blastoutname)

            if re.search(r"metacyc", blastoutname, re.I):
                weight_dbs.append(2)
            else:
                weight_dbs.append(1)

    return database_names, parsed_blastouts, weight_dbs


# the main function
def main(argv, errorlogger=None, runstatslogger=None):
    parser = createParser()
    (opts, args) = parser.parse_args(argv)

    if not check_arguments(opts, args):
        print(usage)
        sys.exit(0)

    results_dictionary = {}
    dbname_weight = {}

    contig_lengths = {}
    read_contig_lengths(opts.contig_map_file, contig_lengths)

    if opts.blastdir != None and opts.sample_name != None:
        try:
            database_names, input_blastouts, weight_dbs = getBlastFileNames(opts)
        except:
            errormod.insert_error(errorcode)
            pass
    else:
        database_names = opts.database_name
        input_blastouts = opts.input_blastout
        weight_dbs = opts.weight_db

    priority = 6000
    count_annotations = {}

    for dbname, blastoutput, weight in zip(database_names, input_blastouts, weight_dbs):
        results_dictionary[dbname] = {}
        dbname_weight[dbname] = weight
        count = process_parsed_blastoutput(
            dbname, weight, blastoutput, opts, results_dictionary[dbname]
        )
        if runstatslogger != None:
            runstatslogger.write(
                "%s\tProtein Annotations from %s\t%s\n"
                % (str(priority), dbname, str(count))
            )
        priority += 1

    for dbname in results_dictionary:
        gutils.eprintf("\n\tINFO:\tNumber of annotations in DB  {}: {}".format(dbname, len(results_dictionary[dbname].keys())))
        for seqname in results_dictionary[dbname]:
            count_annotations[seqname] = True
    count = len(count_annotations)
    if runstatslogger != None:
        runstatslogger.write(
            "%s\tTotal Protein Annotations\t%s\n" % (str(priority), str(count))
        )

    # create the annotations from he results
    create_annotation(
        dbname_weight,
        results_dictionary,
        opts.input_gff,
        opts.rRNA,
        opts.rRNA_gff,
        opts.tRNA,
        opts.tRNA_gff,
        opts.output_gff,
        opts.output_comparative_annotation,
        contig_lengths,
        sample_name=opts.sample_name,
        compact_output=opts.compact_output,
    )

    # Lastly, find overlaps of annotated ORFs against rRNAs and tRNAs
    if opts.rRNA_gff != None: # rRNAs?
        gffc.compare_gffs(opts.output_gff, opts.rRNA_gff, tag='rRNA')
    if opts.tRNA_gff != None: # tRNAs?
        gffc.compare_gffs(opts.output_gff, opts.tRNA_gff, tag='tRNA')



def MetaPathways_annotate_fast(argv, errorlogger=None, runstatslogger=None):
    errorlogger.write("#STEP\tANNOTATE_ORFS\n")
    try:
        main(argv, errorlogger=errorlogger, runstatslogger=runstatslogger)
    except:
        errormod.insert_error(errorcode)
        return (0, "")

    return (0, "")


# the main function of metapaths
if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
