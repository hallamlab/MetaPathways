from __future__ import division
__author__ = "Kishori M Konwar, Aria Hahn"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

import re
import shutil
import pyfastx
import sys
import os
from optparse import OptionParser, OptionGroup
from glob import glob

script_name = sys.argv[0]
usage= script_name + """ --input-folder <inputfolder> --output-folder <output folder> --contigs <contigs file>
                          Takes as the files \"0.pf\", \"genetic-elements.dat\", \"organism-params.dat\" in the <input folder> and 
                          produces inputs for pathway tools for MAGs in the output folder <output folder>"
                     """
parser = OptionParser(usage)

parser.add_option("-i", "--input-folder", dest="input_folder", default =None, 
                  help='the input folder')

parser.add_option("-o", "--output-folder", dest="output_folder", default =None, 
                  help='the output folder')

parser.add_option("--contigs", dest="contigs", default = None, 
                  help='the contigs file')

parser.add_option("-s", "--sample-name", dest="sample_name", default = None, 
                  help='sample name')


def check_arguments(opts, args):
    if opts.input_folder == None or \
       opts.output_folder == None or \
       opts.sample_name == None or \
       opts.contigs == None:
         return False

    return True

def read_pf_file(opts):
     contig_id_in_orfid_patt = re.compile(r'[^_]+_([0-9]+)_[0-9]+')
     remove_patterns = []
     remove_patterns.append(re.compile(r'^MULTISPECIES:\s*'))
 
     contig_no = None
     annot_by_contig = {}
     with open(opts.input_folder + "/0.pf",'r') as fin:
         for line in fin:
            fields = line.strip().split('\t')

            if fields[0]=="FUNCTION":
                funct = fields[1]
                for pat in remove_patterns:
                    funct = re.sub(pat, '', funct)

                if contig_no != None:
                    annot_by_contig[contig_no].append(fields[0] + '\t' +  funct)
                else:
                    raise
            else:
                if fields[0]=="ID":
                    res =  contig_id_in_orfid_patt.search(fields[1])
                    if res: 
                        contig_no = res.group(1)
                        if contig_no not in annot_by_contig: 
                           annot_by_contig[contig_no] = []
                    else:
                        eprintf("WARNING: Cannot determine contig number in the ORF id\n")
                        raise 

                if fields[0]=="//":
                   annot_by_contig[contig_no].append(fields[0])
                else:
                   annot_by_contig[contig_no].append(fields[0] + '\t' + fields[1])

     if not os.path.exists(opts.output_folder):
         os.makedirs(opts.output_folder, exist_ok=True)

     # create the pf file
     for contig_no in annot_by_contig:
        with open(opts.output_folder + "/" + opts.sample_name + "_" + contig_no + ".pf", 'w') as fout:
           fout.write('\n'.join(annot_by_contig[contig_no]))

     # create  the genetic elements file 
     with open(opts.output_folder + "/" + "genetic-elements.dat", 'w') as fout:
         for contig_no in annot_by_contig:
             contig_id = opts.sample_name + "_" + contig_no 
             annot_file = opts.sample_name + "_" + contig_no + ".pf"
             contig_file = opts.sample_name + "_" + contig_no + ".fasta"
             outputlines = ['ID' + '\t' + contig_id]
             outputlines.append('NAME' + '\t' + contig_id)
             outputlines.append('TYPE' + '\t' + ':CONTIG')
             outputlines.append('ANNOT-FILE' + '\t' + annot_file)
             outputlines.append('SEQ-FILE' + '\t' + contig_file)
             outputlines.append('//')
             fout.write('\n'.join(outputlines))
             
     # copy the organism-params.dat file
     shutil.copy2(opts.input_folder + "/organism-params.dat", opts.output_folder + "/" + "organism-params.dat")

     # write the contig files
     contig_id = re.compile(r'[^_]+_([0-9]+)')
     for seq in pyfastx.Fasta(opts.contigs):
        #contents1.append(Sequence(seq.name, seq.seq))
        res =  contig_id.search(seq.name)
        if res: 
            contig_no = res.group(1)
            # write the contig if is is already in the 0.pf file ORF
            if contig_no in annot_by_contig:
                contig_file_path = opts.output_folder + "/" + opts.sample_name + "_" + contig_no + ".fasta"
                with open(contig_file_path, 'w') as fout:
                    fout.write(">" + seq.name + '\n')
                    fout.write(seq.seq + '\n')

# the main function
def main(argv): 
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print(usage)
       sys.exit(0)

    read_pf_file(opts)


#   print  len(orfids), len(neworfs), len(equivalent)
#    print("Original # ORFs   : {}".format(len(orfids)))
#    print("Compressed # ORFs : {}".format(len(neworfs)))
#    print("Equiv # ORFS      : {}".format(len(equivalent)))
#    print("Annotations #     : {}".format(len(annotations)))

def write_new_pf_file(neworfs, annotations, filename):
    outputfile = open(filename,'w')
    for orfid in neworfs: 
       for line in annotations[orfid]:
           fprintf(outputfile, "%s\n",line)
       fprintf(outputfile, "//\n")
    outputfile.close()

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

