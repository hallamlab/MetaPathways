#!/usr/bin/python
# File created on Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re, os
     import sys
     import copy
     from optparse import OptionParser, OptionGroup
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

usage= sys.argv[0] + """ --taxa taxa_file --path pathway_file -o output.csv"""
parser = OptionParser(usage)

parser.add_option( "-i", dest="input_folders", action="append", default=[],  help='input_folder')
parser.add_option( "-s", dest="samples", action="append", default=[],  help='sample name')

parser.add_option( "-o", dest="output_file", help='the output file')
parser.add_option( "-f", dest="category", default="M",  help='the functional category "KEGG" (kegg)/"COG" (cog)/"MetaCyc" (meta) [default M]')
parser.add_option( "-c", dest="category_file", help='the functional category file')

# The C programming language is a way of life
def fprintf(file, fmt, *args):
    file.write(fmt % args)


# check the arguments
def check_arguments(opts, args):
    
    if opts.output_file == None:
         print """Must have an output file"""
         return False
         
    if opts.input_folders == None:
         print """Must have an input folder"""
         return False
    return True

def read_pathways_for_sample(pathway_file, rpkms, all_pwys):
    pathways = {}
    O_patt = re.compile(r'(\d+_\d+)')
    with open(pathway_file,'r') as pathwayfp:
       for line in pathwayfp:
          fields = [ x.strip() for x in line.strip().split('\t') ] 
          if fields[0]=="SAMPLE":
             continue

          all_pwys[fields[1]]= fields[2]
          pwy=fields[1]
          cname=fields[2]
          orfstr= re.sub("\[","",fields[9])
          orfstr= re.sub("\]","",orfstr)
          rpkm_pwy=0
          for orf in orfstr.split(','):
             res = O_patt.search(orf)
             if res:
                orfname = res.group(1)
                if orfname in rpkms:
                   rpkm_pwy += rpkms[orfname]
           
          pathways[pwy] =rpkm_pwy
    
    return pathways


def read_rpkm_for_sample(rpkm_file):
    # try to open pathway file
    rpkms={}
    with open(rpkm_file,'r') as rpkmfp:
       for line in rpkmfp:
          fields = [ x.strip() for x in line.strip().split('\t') ]

          if len(fields) ==2:
            try:
               count = float(fields[1])
               rpkms[fields[0]]= count
            except:
               pass
    return rpkms
    
def write_pathway_rxn_taxa_table(read_to_taxa,pathway_rxns,output_file):
    # try to open pathway file
    try:
        mypath_output_file = open(output_file,'w')
    except IOError:
        print "Cannot open " + str(output_file)
        
    
    
    mypath_output_file.close()


# the main function
def main(argv): 
    # grab and check arguments
    (opts, args) = parser.parse_args() 
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)

    if opts.category=="MetaCyc":    
       create_metacyc_table(opts)
    else:
       create_cogg_table(opts, type=opts.category)

def read_category_table(opts):

    if opts.category=='COG':
      categoryPatt =re.compile('^\\t[^\\t]')
      subcategoryPatt =re.compile('^\\t\\t[^\\t]')

    if opts.category=='KEGG':
      categoryPatt =re.compile('^\\t\\t[^\\t]')
      subcategoryPatt =re.compile('^\\t\\t\\t[^\\t]')

    subcategory_to_category0 = {}
    category_counts0 = {}
    category_counts1 = {}
    with open(opts.category_file, 'r') as fin:
        category = None
        subcategory = None
        for line in fin:
           line = line.rstrip() 
           if categoryPatt.search(line):
              fields = [ x.strip() for x in line.strip().split('\t') ]
              category_counts0[fields[0]] = [ fields[1], 0 ]
              category = fields[0]

           if subcategoryPatt.search(line):
              fields = [ x.strip() for x in line.strip().split('\t') ]
              subcategory_to_category0[fields[0]] = category
              category_counts1[fields[0]] = [ fields[1], 0 ]

    subcategory_to_category1 = { x:x for x in subcategory_to_category0.keys() }
    return subcategory_to_category0, subcategory_to_category1, category_counts0, category_counts1

def get_sample_names(opts):
    samples={}
    for folder in opts.input_folders:
      for file in os.listdir(folder):
        if os.path.exists(folder + "/" +  file) and os.path.exists(folder + "/" +  file + "/preprocessed" ):
           samples[file] = True


    return samples.keys()


def create_cogg_table(opts, type):

     
    if opts.samples==[]: 
       sample_names  = get_sample_names(opts)
    else:
       samples_names = opts.samples
    
    match_to_category0, match_to_category1, empty_category_counts0, empty_category_counts1= read_category_table(opts)

    samples = read_rpkm_by_category(opts.category, opts.input_folders, sample_names, empty_category_counts0, match_to_category0)
    write_cog_output(opts.output_file + "_level0.txt", sample_names, empty_category_counts0, samples)

    samples = read_rpkm_by_category(opts.category, opts.input_folders, sample_names, empty_category_counts1, match_to_category1)
    write_cog_output(opts.output_file + "_level1.txt", sample_names, empty_category_counts1, samples)

def read_rpkm_by_category(category, input_folders, sample_names, empty_category_counts, match_to_category):
    commentPatt =re.compile('^#')
    all_pathways = {}
    samples ={}
    file = 0
    for sample in sample_names:
      for folder in input_folders:
         rpkm_file = folder + "/" + sample + "/results/rpkm/" + sample  + ".orfwise"
         annotation_file = folder + "/" + sample + "/results/annotation_table/" + sample + ".ORF_annotation_table.txt" 

         if os.path.exists(rpkm_file):
           rpkms = read_rpkm_for_sample(rpkm_file)

           file += 1
           print annotation_file, file
          
           mis =0
           with open(annotation_file, 'r') as fin:
              category_counts = copy.deepcopy(empty_category_counts)
              for line in fin:
                 if commentPatt.search(line):
                    fields = [ x.strip() for x in line.strip().split('\t') ]
                    index=0
                    for field in fields:
                      if field==category:
                         break
                      index += 1
                    continue
                
                 fields = [ x.strip() for x in line.strip().split('\t') ]
               
                 if not fields[0] in rpkms:
                    continue

                 if len(fields) <= index:
                    mis += 1
                    continue

                 if not fields[index] in match_to_category:
                    continue
                 category_counts[ match_to_category[fields[index]] ][1] += rpkms[fields[0]]

           samples[sample] =copy.deepcopy(category_counts)
           if mis:
             print "# Missing lines", mis

    return samples
    
   
def write_cog_output(output_file, sample_names, empty_category_counts, samples):
    with open(output_file, "w") as fw:
     fprintf(fw, "FUNCTION\tCOMMON-NAME")
     for sample in sample_names:
       fprintf(fw, "\t%s", sample)
     fprintf(fw, "\n")

     for cat in sorted(empty_category_counts.keys()):
       fprintf(fw, "%s\t%s",cat, empty_category_counts[cat][0])
       for sample in sample_names:
          if cat in samples[sample]:
              fprintf(fw, "\t%.2f",samples[sample][cat][1])
          else:
              fprintf(fw, "\t0.0")
       fprintf(fw, "\n")


def create_metacyc_table(opts):
    if opts.samples==[]: 
       sample_names  = get_sample_names(opts)
    else:
       sample_names = opts.samples
    

    samples ={}
    all_pathways = {}

    for sample in sample_names:
      for folder in opts.input_folders:
         rpkm_file = folder + "/" + sample + "/results/rpkm/" + sample  + ".orfwise"
         pathway_file = folder + "/" + sample + "/results/pgdb/" + sample.lower() + ".pwy.txt" 

         if os.path.exists(rpkm_file) and os.path.exists(pathway_file):
           print rpkm_file, pathway_file
           rpkms = read_rpkm_for_sample(rpkm_file)
           samples[sample] = read_pathways_for_sample(pathway_file, rpkms, all_pathways)

           

    with open(opts.output_file, "w") as fw:
     fprintf(fw, "PATHWAY\tCOMMON-NAME")
     for sample in sample_names:
       fprintf(fw, "\t%s", sample)
     fprintf(fw, "\n")

     for pwy in sorted(all_pathways.keys()):
       fprintf(fw, "%s\t%s",pwy, all_pathways[pwy])
       for sample in sample_names:
          if pwy in samples[sample]:
              fprintf(fw, "\t%.2f",samples[sample][pwy])
          else:
              fprintf(fw, "\t0.0")
       fprintf(fw, "\n")

    # create a read to taxa dictionary
    #read_to_taxa = process_read_taxa_file(opts.taxa_file)
    
    # create a MetaCyc pathway shortname to longname hash
    #pathway_short_to_long =  process_path_short_long(opts.pathway_file)
    
    # collect pathway names, reads, and reactions
    #pathways = process_pathways(opts.pathway_file)
    
    #if(opts.path_rxn_file != None):
    #    pathway_rxns = process_pathway_reaction(opts.path_rxn_file)
    #    write_pathway_rxn_taxa_table(read_to_taxa,pathway_rxns,opts.output_file)
    
    
    # create table of pathways vs taxa in the output file
    # write_pathways_taxa(read_to_taxa, pathway_short_to_long, pathways, opts.output_file)
    
# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

