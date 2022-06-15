import re

from os import makedirs, sys, remove, rename
from sys import path
from optparse import OptionParser
import re

usage= sys.argv[0] + """ -a <unfiltered annot file> -A <filtered annot file> -p <0.pf> """

def createParser():
    epilog = """This script cleans up annotations in a .pf file."""

    epilog = re.sub(r'[ \t\f\v]+',' ', epilog)

    parser = OptionParser(usage=usage, epilog=epilog)

    parser.add_option("-a", "--annot-unfiltered", dest="unfiltered_annot",
                      help='functional annotation that are not filtered')

    parser.add_option("-A", "--annot-filtered", dest="filtered_annot",
                      help='functional annotation that are filtered')

    parser.add_option("-p", "--pf-file", dest="pf_file",
                      help='pf file')
    return parser


def main(argv): 
    parser = createParser()
    (opts, args) = parser.parse_args(argv)

    unfiltered_filtered = {}
    with open(opts.unfiltered_annot, 'r') as fraw, open(opts.filtered_annot, 'r') as fclean:
        for raw, clean in zip(fraw, fclean):
            unfiltered_filtered[raw.split('\t')[1].strip()] = clean.split('\t')[1].strip()

    with open(opts.pf_file, 'r') as fin:
        for _line in fin:
           line = _line.strip()
           fields = [x.strip() for x in line.split('\t') ]
           if fields[0] == "FUNCTION" and len(fields) == 2:
               if fields[1] in unfiltered_filtered:
                  print('\t'.join([fields[0], unfiltered_filtered[fields[1]]]))
               else:
                  print('\t'.join([fields[0], fields[1]]))
           else:
               print(line)



if __name__ == "__main__":
    main(sys.argv[1:])

