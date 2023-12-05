"""Computes rpkm by aligning the read to the contigs """

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import copy
    import optparse
    import sys
    import re
    import glob
    import multiprocessing
    import shutil
    import subprocess
    
    from os import path, _exit, rename, system

    from metapathways import errorcodes as errormod
    from metapathways import metapathways_utils as mputils
    from metapathways import general_utils as gutils
    from metapathways import sysutil as sysutils
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)

PATHDELIM = sysutils.pathDelim()

epilog = (
    """\n"""
    + """
                This script computes the RPKM values for each ORF, from the BWA recruits.  The input reads
                (in the form of fastq files)  for this step must be added to the subdirectory reads in the
                input folder (where the input fasta  files are located). The read file sare identified by
                the name format of the files: For examples, if the sample name is "abcd" then the following
                read files in the "reads" folders associated with the samples abcd:

                1.   abcd.fastq : this means non-paired reads

                2.   abcd.b1.fastq  : means only unpaired read from batch b1

                3.   abcd_R1.fastq  and abcd_R2.fastq: this means paired reads for sample

                4.   abcd_R1.fastq or abcd_R2.fastq: this means only one end of a paired read

                5.   abcd_R1.b2.fastq and  abcd_R2.b2.fastq: this means paried reads from batch b2, note that batches are idenfied as bn, where n is a number

                6.   abcd_R1.b1.fastq or abcd_R2.b1.fastq: this means only one of a paried read from batch b1
             """
)


usage = (
    sys.argv[0]
    + """ -c <contigs> -o <output> -r <reads>  -g <orfgff> --rpkmExec <rpkmexec> """
    + epilog
)

def createParser():
    parser = optparse.OptionParser(usage=usage)

    # Input options
    parser.add_option(
        "-c", "--contigs", dest="contigs",
        default=None, help="the contigs file"
        )
    parser.add_option(
        "-o", "--output", dest="output",
        default=None, help="orfwise read count file"
        )
    parser.add_option(
        "--stats", dest="stats", default=None,
        help="output stats for ORFs into file"
        )
    parser.add_option(
        "-1", "--fwd_fastq", dest="fwd_fastq",
        default=None, help=""
        )
    parser.add_option(
        "-2", "--rev_fastq", dest="rev_fastq",
        default=None, help=""
        )
    parser.add_option(
        "--interleaved", dest="interleaved",
        action="store_true", default=False,
        help=""
        )
    #parser.add_option(
    #    "-r", "--readsdir", dest="readsdir",
    #    default=None, help="the directory that should have the read files"
    #    )
    parser.add_option(
        "-g", "--orfgff", dest="orfgff",
        default=None, help="annotated GFF file"
        )
    parser.add_option(
        "-s", "--sample_name", dest="sample_name",
        default=None, help="name of the sample"
        )
    parser.add_option(
        "--rpkmExec", dest="rpkmExec",
        default=None, help="RPKM Executable"
        )
    parser.add_option("--bwaExec", dest="bwaExec",
        default=None, help="BWA Executable"
        )
    parser.add_option("--bwaFolder", dest="bwaFolder",
        default=None, help="BWA Folder"
        )
    parser.add_option(
        "-n", "--num_threads", dest="num_threads",
        default=1, help="number of threads",
        )

    return parser


def getSamFiles(readdir, sample_name):
    """This function finds the set of SAM files that has the BWA recruitment information"""

    samFiles = glob.glob(readdir + PATHDELIM + sample_name + "*.sam")

    return samFiles


def getBamFiles(readdir, sample_name):
    """This function finds the set of SAM files that has the BWA recruitment information"""

    bamFiles = glob.glob(readdir + PATHDELIM + sample_name + "*.bam")

    return bamFiles


def indexForBWA(bwaExec, contigs, indexfile):
    cmd = "%s index -p %s %s" % (
        bwaExec,
        indexfile,
        contigs,
    )

    result = sysutils.getstatusoutput(cmd)

    if result[0] == 0:
        return True

    return False


def runUsingBWA(bwaExec, sample_name, indexFile, readgroup, readFiles, bwaFolder, num_threads):

    status = True

    bwaOutput = bwaFolder + PATHDELIM + readgroup + ".sam"
    stOutput = bwaFolder + PATHDELIM + readgroup + ".bam"

    stOutputTmp = stOutput + ".tmp"
    stInterTmp = bwaFolder + PATHDELIM + readgroup + ".tmp"

    cmd = "command not prepared"

    if len(readFiles) == 2:
        bwa_cmd = "%s mem -t %d -o %s %s %s %s 2> /dev/null" % (
            bwaExec,
            num_threads,
            bwaOutput,
            indexFile,
            readFiles[0],
            readFiles[1],
        )

    if len(readFiles) == 1:
        bwa_cmd = "%s mem -t %d -o %s %s %s 2> /dev/null" % (
            bwaExec,
            num_threads,
            bwaOutput,
            indexFile,
            readFiles[0],
        )

    st_cmd = "samtools sort -O bam -o %s -T %s %s" % (
        stOutputTmp,
        stInterTmp,
        bwaOutput,
    )
    
    bwaResult = sysutils.getstatusoutput(bwa_cmd)
    stResult = sysutils.getstatusoutput(st_cmd)

    clean_cmd = st_cmd = "rm %s" % (bwaOutput)

    if stResult[0] == 0:
        rename(stOutputTmp, stOutput)
        cleanup = sysutils.getstatusoutput(clean_cmd)

    else:
        gutils.eprintf("ERROR:\tError in file processing read files %s\n", readFiles)
        status = False

    return status

def read_genome_equivalent(microbecensusoutput):
    gen_equiv_patt = re.compile(r"genome_equivalents:\s+(.*)$")

    with open(microbecensusoutput, "r") as inputfile:
        lines = inputfile.readlines()
        for line in lines:
            result = gen_equiv_patt.search(line)
            if result:
                genome_equivalent = result.group(1)
                try:
                    return float(genome_equivalent)
                except:
                    return 1

    return 1

def getReadFiles(readdir, sample_name):
    """This function finds the set of fastq files that has the reads"""
    _fastqfiles = glob.glob(
       readdir + PATHDELIM + sample_name + "_*"
    )

    fastqgroups = {}
    for _fastqfile in _fastqfiles:
       fastqfile = re.sub(r'.gz$','', _fastqfile, flags=re.IGNORECASE) 
       fastqfile = re.sub(r'.fastq$','', fastqfile, flags=re.IGNORECASE) 
       fastqfile = re.sub(r'.fq$','', fastqfile, flags=re.IGNORECASE) 
       trimmedfastq = path.basename(re.sub(r'_[12]$', '', fastqfile))
       if len(trimmedfastq.rsplit('.', 1)) > 1: # goofy hack to get batches to work, need to improve
           trimmedfastq = re.sub(r'_[12]', '', trimmedfastq.rsplit('.', 1)[0]) + '.' + trimmedfastq.rsplit('.', 1)[1]

       if trimmedfastq not in fastqgroups:
           fastqgroups[trimmedfastq] = []
  
       fastqgroups[trimmedfastq].append(_fastqfile)

    keys = fastqgroups.keys()
    for key in keys:
        fastqgroups[key].sort(key = lambda x: x)

    return fastqgroups


def main(argv, errorlogger=None, runcommand=None, runstatslogger=None):
    parser = createParser()

    options, args = parser.parse_args(argv)
    if not (options.contigs != None and path.exists(options.contigs)):
        parser.error("ERROR\tThe contigs file is missing")
        errormod.insert_error(10)
        return 1

    if not (options.orfgff != None and path.exists(options.orfgff)):
        parser.error("ERROR\tThe GFF annotation file is missing")
        errormod.insert_error(10)
        return 1

    cmd_exists = lambda x: shutil.which(x) is not None
    if not (options.rpkmExec != None and cmd_exists(options.rpkmExec)):
        parser.error("ERROR\tThe RPKM executable is missing")
        errormod.insert_error(10)
        return 1

    if not (options.bwaExec != None and cmd_exists(options.bwaExec)):
        parser.error("ERROR\tThe BWA executable is missing")
        errormod.insert_error(10)
        return 1

    if ((options.fwd_fastq != None and options.fwd_fastq != 'None') and not path.exists(options.fwd_fastq)):
        parser.error("ERROR\tThe expected FASTQ \'{}\' is missing.".format(options.fwd_fastq))
        errormod.insert_error(10)
        return 1

    if ((options.rev_fastq != None and options.rev_fastq != 'None') and not path.exists(options.rev_fastq)):
        parser.error("ERROR\tThe expected FASTQ \'{}\' is missing.".format(options.rev_fastq))
        errormod.insert_error(10)
        return 1

    #if not (options.readsdir != None and path.exists(options.readsdir)):
    #    parser.error("ERROR\tThe expected RPKM directory \'{}\' is missing.".format(options.readsdir))
    #    errormod.insert_error(10)
    #    return 1

    if not (options.bwaFolder != None and path.exists(options.bwaFolder)):
        parser.error("ERROR\tThe BWA directory is missing.")
        errormod.insert_error(10)
        return 1

    if options.sample_name == None:
        parser.error("ERROR\tThe sample name is missing.")
        errormod.insert_error(10)
        return 1

    #readFiles = getReadFiles(options.readsdir, options.sample_name)
    command = [options.rpkmExec]
    command.append('contig')
    command.append('-m length')
    command.append('-m count')
    command.append('-m mean')
    command.append('-m variance')
    command.append('-m trimmed_mean')
    command.append('-m rpkm')
    command.append('-m tpm')
    # set up fastqs
    if options.fwd_fastq and options.rev_fastq:
        command.append('-1')
        command.append(options.fwd_fastq)
        command.append('-2')
        command.append(options.fwd_fastq)
    elif options.fwd_fastq and options.interleaved:
        command.append('--interleaved')
        command.append(options.fwd_fastq)
    elif options.fwd_fastq and not options.interleaved:
        command.append('--single')
        command.append(options.fwd_fastq)
    else:
        parser.error("ERROR\tFASTQs not specified correctly.")
        errormod.insert_error(10)
        return 1
    command.append('-r')
    command.append(options.contigs)
    command.append('--output-file')
    command.append(options.output)
    command.append('-t')
    command.append(options.num_threads)
    command.append('--bam-file-cache-directory')
    command.append(options.bwaFolder)
    command.append('--min-read-percent-identity 97')
    command.append('--min-read-aligned-percent 97')
    command.append('--exclude-supplementary')

    rpkmstatus = 0
    rpkmtext = ''
    try:
        rpkmstatus, rpkmtext = runRPKMCommand(runcommand=" ".join(command))
    except:
        rpkmstatus = 1
        pass

    with open(options.stats, 'w') as stat_out:
        stat_out.write(' '.join(command))
        stat_out.write('\n\n')
        stat_out.write(rpkmtext)
        stat_out.write('\n\n')

        #Start custom abundance analysis
        bam_in = glob.glob(options.bwaFolder + PATHDELIM + '*.bam')[0]
        sort_out = options.bwaFolder + PATHDELIM + options.sample_name + '.sorted.bam'
        tmp_bam = options.bwaFolder + PATHDELIM + options.sample_name + '.tmp'
        gff_in = options.orfgff
        gtf_out = path.join(path.dirname(options.orfgff), path.basename(options.orfgff).rsplit('.', 1)[0] + '.annot.gtf')
        feat = options.bwaFolder + PATHDELIM + options.sample_name + '.annot.featurecounts.txt'
        gen_lens = options.bwaFolder + PATHDELIM + options.sample_name + '.annot.genelengths.txt'
        abund_calc = path.join(path.dirname(options.output), path.basename(options.output).rsplit('.', 2)[0] + '.orf_counts.tsv')

        samcmd = ['samtools', 'sort', '-n', '-T', f'{tmp_bam}', '-o', f'{sort_out}', f'{bam_in}']
        stat_out.write(' '.join(samcmd))
        stat_out.write('\n\n')
        sam_result = subprocess.Popen(' '.join(samcmd), stdout=stat_out, stderr=subprocess.PIPE, shell=True)    
        stat_out.write('\n')
        sam_result.wait()

        gtfcmd = ['gff2gtf.py', '--feature', 'CDS,rRNA,tRNA,pseudogene', f'{gff_in}', f'{gtf_out}']
        stat_out.write(' '.join(gtfcmd))
        stat_out.write('\n')
        gtf_result = subprocess.Popen(' '.join(gtfcmd), stdout=stat_out, stderr=subprocess.PIPE, shell=True)    
        stat_out.write('\n')
        gtf_result.wait()

        featcmd = ['featureCounts', '-t', 'CDS,rRNA,tRNA,pseudogene', '-p', '-O', '-T', f'{options.num_threads}', '-a', f'{gtf_out}', '-o', f'{feat}', f'{sort_out}']
        stat_out.write(' '.join(featcmd))
        stat_out.write('\n')
        feat_result = subprocess.Popen(' '.join(featcmd), stdout=stat_out, stderr=subprocess.PIPE, shell=True)    
        stat_out.write('\n')
        feat_result.wait()

        cutcmd = f'cut -f4,5,9 {gtf_out} | sed \'s/gene_id //g\' | gawk \'{{print $3,$2-$1+1}}\' | tr \' \' \'\t\' > {gen_lens}'
        stat_out.write(cutcmd)
        stat_out.write('\n')
        cut_result = subprocess.Popen(cutcmd, stdout=stat_out, stderr=subprocess.PIPE, shell=True)    
        stat_out.write('\n')
        cut_result.wait()    

        abuncmd = ['abund_calc.py', '--counts-file', f'{feat}', '--gene-lengths-file', f'{gen_lens}', '--gtf-file', f'{gtf_out}', '--output', f'{abund_calc}']
        stat_out.write('\n')
        stat_out.write(' '.join(abuncmd))
        stat_out.write('\n')
        abund_result = subprocess.Popen(' '.join(abuncmd), stdout=stat_out, stderr=subprocess.PIPE, shell=True)    
        stat_out.write('\n')
        abund_result.wait()

    if rpkmstatus != 0:
        gutils.eprintf("ERROR\tRPKM calculation was unsuccessful\n")
        errormod.insert_error(10)
        return 1
        # exit_process("ERROR\tFailed to run RPKM" )

    return rpkmstatus


def runRPKMCommand(runcommand=None):
    if runcommand == None:
        return False

    result = sysutils.getstatusoutput(runcommand)
    #if result[1]:
    #    print(result[1])
    return result[0], result[1]


# this is the portion of the code that fixes the name


def split_attributes(str, attributes):
    rawattributes = re.split(";", str)
    for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)

    return attributes


# this is the function that fixes the name
def fix_pgdb_input_files(pgdb_folder, pgdbs=[]):
    pgdb_list = glob.glob(pgdb_folder + "/*/input/organism.dat")

    for pgdb_organism_file in pgdb_list:
        process_organism_file(pgdb_organism_file)


def fixLine(line, id):
    fields = line.split("\t")
    if len(fields) == 2:
        return fields[0] + "\t" + id


def getID(line):
    fields = line.split("\t")
    if len(fields) == 2:
        return fields[1]


def process_organism_file(filel):
    patternsToFix = [
        re.compile(r"NAME\tunclassified sequences"),
        re.compile(r"ABBREV-NAME\tu. sequences"),
    ]
    patternID = re.compile(r"^ID\t.*")
    try:
        orgfile = open(filel, "r")
    except IOError:
        print("ERROR : Cannot open organism file" + str(filel))
        errormod.insert_error(10)
        return

    lines = orgfile.readlines()
    newlines = []
    needsFixing = False

    id = None
    for line in lines:
        line = line.strip()
        if len(line) == 0:
            continue
        flag = False

        result = patternID.search(line)
        if result:
            id = getID(line)

        for patternToFix in patternsToFix:
            result = patternToFix.search(line)
            if result and id:
                newline = fixLine(line, id)
                newlines.append(newline)
                flag = True
                needsFixing = True

        if flag == False:
            newlines.append(line)

    orgfile.close()
    if needsFixing:
        write_new_file(newlines, filel)


def write_new_file(lines, output_file):

    print("Fixing file " + output_file)
    try:
        outputfile = open(output_file, "w")
        pass
    except IOError:
        print("ERROR :Cannot open output file " + output_file)

    for line in lines:
        gutils.fprintf(outputfile, "%s\n", line)

    outputfile.close()


def MetaPathways_tpm(argv, extra_command=None, errorlogger=None, runstatslogger=None):
    if errorlogger != None:
        errorlogger.write("#STEP\tRPKM_CALCULATION\n")
    try:
        main(
            argv,
            errorlogger=errorlogger,
            runcommand=extra_command,
            runstatslogger=runstatslogger,
        )
    except:
        errormod.insert_error(10)
        return (1, traceback.print_exc(10))

    return (0, "")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
