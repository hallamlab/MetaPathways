"""Script creates
   a) ptools input
   b) genbank file with annotations
"""

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import optparse
    import shutil
    import sys
    import re
    import gzip
    import glob
    import pandas as pd
    import numpy as np
    from tqdm import tqdm

    from os import makedirs, path, listdir, remove, rename

    from metapathways import errorcodes as errormod
    from metapathways import general_utils as gutils
    from metapathways import metapathways_utils as mputils
    from metapathways import sysutil as sysutils
    from metapathways import parse as parsemod
except:
    print(""" Could not load some user defined  module functions""")
    print(""" """)
    sys.exit(3)

PATHDELIM = sysutils.pathDelim()

def createParser():
    epilog = """
    This script has three functions:
    (i) The functional and taxonomic annotations created for the individual
    ORFs are used to create the inputs required by the Pathway-Tools's
    Pathologic algorithm to build the ePGDBs.
    The input consists of 4 files that contains functional annotations and
    sequences with relevant information, this information is used by
    Pathologic to create the ePGDBs.
    (ii) It can create a genbank file for the ORFs and their annotations
    """
    epilog = re.sub(r'\s+', ' ', epilog)

    parser = optparse.OptionParser(usage=usage, epilog = epilog)

    # Input options

    input_group = optparse.OptionGroup(parser, 'input options')

    input_group.add_option('-g', '--gff', dest='gff_file',
                           metavar='GFF_FILE',
                           help='GFF files, with annotations,  to convert to genbank format')

    input_group.add_option('-n', '--nucleotide', dest='nucleotide_sequences',
                           metavar='NUCLEOTIDE_SEQUENCE',
                           help='Nucleotide sequences')

    input_group.add_option('-p', '--protein', dest='protein_sequences',
                           metavar='PROTEIN_SEQUENCE',
                           help='Protein sequences')


    input_group.add_option('-o', '--output', dest='output',
                           metavar='OUTPUT', help='Genbank file')

    parser.add_option_group(input_group)

    output_options_group =  optparse.OptionGroup(parser, 'Output Options')

    output_options_group.add_option("--out-gbk", dest="gbk_file", default = None,
                     help='option to create a genbank file')

    output_options_group.add_option("--ncbi-tree", dest="ncbi_taxonomy_tree",  default=None,
                       help='add the ncbi taxonomy tree ')

    output_options_group.add_option("--taxonomy-table", dest="taxonomy_table",  default=None,
                       help='table with taxonomy')

    output_options_group.add_option("--annotation-table", dest="annotation_table",  default=None,
                       help='table with all annotations')
    
    output_options_group.add_option("--ec-mapping", dest="ec_mapping",  default=None,
                       help='table with mappings for Uniprot/KEGG/EC')

    output_options_group.add_option("--rxn-mapping", dest="rxn_mapping",  default=None,
                       help='table with mappings for MetaCyc PROT IDs to RXN IDS')

    output_options_group.add_option("--out-ptinput", dest="ptinput_file",  default=None,
                     help='option and directory  to create ptools input files')

    output_options_group.add_option("--compact-output", dest="compact_output",  default=True, action='store_true',
                     help='option to create compact orfid names')

    parser.add_option_group(output_options_group)

    return parser

def insert_attribute(attributes, attribStr):
     rawfields = attribStr.split('=', 1)
     if len(rawfields) == 2:
       attributes[rawfields[0].strip().lower()] = rawfields[1].strip()

def split_attributes(str, attributes):
     rawattributes = re.split(';', str)
     for attribStr in rawattributes:
        insert_attribute(attributes, attribStr)
     return attributes


def insert_orf_into_dict(line, contig_dict):
     rawfields = re.split('\t', line)
     fields = []
     for field in rawfields:
        fields.append(field.strip());

     if( len(fields) != 9):
       return

     attributes = {}
     try:
         attributes['seqname'] =  fields[0]   # this is a bit of a  duplication
         attributes['source'] =  fields[1]
         attributes['feature'] =  fields[2]
         attributes['start'] =  int(fields[3])
         attributes['end'] =  int(fields[4])
     except:
         print(line)
         sys.exit(0)

     try:
        attributes['score'] =  float(fields[5])
     except:
        attributes['score'] =  fields[5]

     attributes['strand'] =  fields[6]
     attributes['frame'] =  fields[7]
     split_attributes(fields[8], attributes)

     if not fields[0] in contig_dict :
       contig_dict[fields[0]] = []

     contig_dict[fields[0]].append(attributes)


def get_sequence_name(line):
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     #print name
     return name


def get_sequence_number(line):
     fields = re.split(' ', line)
     name = re.sub('>','',fields[0])
     seqnamePATT = re.compile(r'[\S]+_(\d+)$')
     result = seqnamePATT.search(line.strip())
     return result.group(1)

     #print name


note = """GFF File Format
Fields

Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
"""

#def get_date():


def get_sample_name(gff_file_name):
     sample_name= re.sub('.annotated.gff', '', gff_file_name)
     sample_name= re.sub('.annot.gff', '', gff_file_name)
     sample_name= re.sub(r'.*[/\\]', '', sample_name)
     return sample_name


def process_gff_file(gff_file_name, output_filenames, nucleotide_seq_dict, \

    protein_seq_dict, input_filenames, orf_to_taxonid = {},  compact_output=True):

    gff_file_name = gutils.correct_filename_extension(gff_file_name)
    with gzip.open(gff_file_name, 'r') if gff_file_name.endswith('.gz') \
        else open(gff_file_name, 'r') as gfffile:
        sample_name=get_sample_name(gff_file_name)

        gff_lines = gfffile.readlines()
        gff_beg_pattern = re.compile("^#")

        contig_dict={}
        count = 0

        for line in gff_lines:
           line = line.strip()
           if gff_beg_pattern.search(line):
             continue
           """  Do not add tRNA """
           try:
              insert_orf_into_dict(line, contig_dict)
           except:
              gutils.eprintf("%s\n", raceback.print_exc(10))

        if "gbk" in output_filenames:
          write_gbk_file(  # TODO: this function should be parallelized
              output_filenames['gbk'], contig_dict, sample_name, \
              nucleotide_seq_dict, protein_seq_dict
              )

        if "ptinput" in output_filenames:
          write_ptinput_files(  # TODO: this function should be parallelized
              output_filenames, contig_dict, sample_name, \
              nucleotide_seq_dict, protein_seq_dict, compact_output, orf_to_taxonid=orf_to_taxonid
          )

# this function creates the pathway tools input files
def write_ptinput_files(outfiles, contig_dict, sample_name, nucleotide_seq_dict, \
        protein_seq_dict, compact_output, orf_to_taxonid={}):
     
    # check if EC available in attributes
    ec_map = {} 
    for key in contig_dict:
        for attrib in contig_dict[key]:
            if attrib['feature'] == 'CDS':
                ec = attrib['ec']
                if ec:
                    ec_map[attrib['target']] = [str(ec)]
    
    # Get all required mapping files
    output_dir_name = outfiles['ptinput']
    ec_map_files = outfiles['ec_mapping']
    rxn_map_file = outfiles['rxn_mapping']
    anno_file = outfiles['annotation_table']
    ft_file = outfiles['taxonomy_table']
    anno_dir_path = path.dirname(ft_file)

    try:
        gutils.remove_dir(output_dir_name)
        makedirs(output_dir_name)
        genetic_elements_file = open(output_dir_name + "/.tmp.genetic-elements.dat", 'w')
    except:
        print("cannot create the pathway tools files")
        print("perhaps there is already a folder " + output_dir_name)
        traceback.print_exc(file=sys.stdout)

    # load annotable
    anno_df = pd.read_csv(anno_file, sep='\t', header=0, low_memory=False)
    anno_df.loc[anno_df['orf_id'] == '','orf_id'] = np.nan
    anno_df['orf_id']  = anno_df['orf_id'].ffill()
    anno_df.dropna(subset = ['target'], inplace=True)
    
    # load funct/tax table
    ft_df = pd.read_csv(ft_file, sep='\t', header=0, low_memory=False)
    ft_df.dropna(subset = ['target'], inplace=True)
    ft_df['orf_id'] = [x for x in ft_df['ORF_ID']]
    # build product dictionary to remove duplicate ORF annotations
    cprod_dict = {k : [] for k in set(ft_df['product'])}
    rep_orf_list = list(ft_df.drop_duplicates('product')['orf_id'])

    # modify target headers for later mapping
    trim_t_list = []
    for d,t in zip(anno_df['ref dbname'], anno_df['target']):
        if t != "<no accession>":
            if 'metacyc' in d:
                trim_t = t.split('|', 2)[2]
                trim_t_list.append(trim_t)
            elif 'swissprot' in d:
                trim_t = t.split('|', 2)[1]
                trim_t_list.append(trim_t)
            elif 'uniref' in d:
                trim_t = t.split('_', 1)[1]
                trim_t_list.append(t) #trim_t)
            elif 'eggnog' in d:
                trim_t = t
                trim_t_list.append(t)
            else:
                trim_t_list.append('NONE')
        else:
                trim_t_list.append('NONE')
    anno_df['trim_target'] = trim_t_list
    # remove uninformative targets
    anno_sub_df = anno_df.query("trim_target != 'NONE'")
    # subset annotable to only rep ORFs
    anno_sub_df = anno_sub_df.copy().query('orf_id in @rep_orf_list')

    trim_set_list = list(set(trim_t_list))
    # Load all mapping files, subset to only annotated targets, create dicts
    ec_dfs = []
    if ec_map_files:
        for ec_map_file in glob.glob(ec_map_files):
            ec_map_df = pd.read_csv(ec_map_file, sep='\t', header=0, low_memory=False)
            ec_dfs.append(ec_map_df)
        ec_df = pd.concat(ec_dfs)
        ec_df.dropna(subset = ['EC'], inplace=True)
        ec_df = ec_df.copy().query('AccID in @trim_set_list')
        ec_df = ec_df[~ec_df['EC'].str.contains(".-")]
        ec_df['EC'] = ec_df['EC'].str.replace('EC:', '')
        ec_dict = ec_df.groupby('AccID', group_keys=True)['EC'].apply(list).to_dict()
    else:
        ec_dict = {}
    if path.exists(rxn_map_file):
        rxn_df = pd.read_csv(rxn_map_file, sep='\t', header=0, low_memory=False)
        rxn_df = rxn_df.query('MC in @trim_set_list')
        rxn_dict = rxn_df.groupby('MC')['RXN'].apply(list).to_dict()
    else:
        rxn_dict = {}
    # Merge the ec_map from earlier attrib search and ec_dict
    ec_m_dict = ec_dict | ec_map

    # Add EC and RXN to anno table
    anno_df['EC'] = ['|'.join(list(set(ec_m_dict[x]))) if x in ec_m_dict else np.nan for x in anno_df['trim_target']]
    anno_df['RXN'] = ['|'.join(list(set(rxn_dict[x]))) if x in rxn_dict else np.nan for x in anno_df['trim_target']]
    anno_df.to_csv(anno_dir_path + "/" + sample_name + ".EC_RXN_map.tsv", sep='\t', index=False)
    pt_attrib_dict = {}
    with open(output_dir_name + "/" + "0.pf", 'w') as pfFile:
        # iterate over every contig sequence
        if compact_output:
            prefix = '' #'O_'
        else:
            prefix = sample_name + '_'

        countError = 0
        for key in contig_dict:
            first = True
            for attrib in contig_dict[key]:
                id  =  attrib['id']
                shortid = ""
                compactid = ""
                if attrib['feature'] == 'CDS':
                    shortid = prefix + attrib['id']
                    compactid = attrib['id']
                    # clean function
                    l_func = attrib['product']
                    if l_func != '':
                        try:
                            protein_seq = protein_seq_dict[id]
                        except:
                            protein_seq = ""
                        try:
                            if ((attrib['product']=='hypothetical protein') |
                                (attrib['product']=='<unannotated protein>')):
                                continue
                        except:
                            print(attrib)
                            sys.exit(0)
                        # find all DBs and map to general DB names
                        orf_anno_df = anno_sub_df.query('orf_id == @compactid')
                        db_dict = {}
                        for dbname in set(orf_anno_df['ref dbname']):
                            if 'metacyc' in dbname:
                                if 'metacyc' in db_dict:
                                    db_dict['metacyc'].append(dbname)
                                else:
                                    db_dict['metacyc'] = [dbname]
                            elif 'swissprot' in dbname:
                                if 'swissprot' in db_dict:
                                    db_dict['swissprot'].append(dbname)
                                else:
                                    db_dict['swissprot'] = [dbname]
                            elif 'uniref' in dbname:
                                if 'uniref' in db_dict:
                                    db_dict['uniref'].append(dbname)
                                else:
                                    db_dict['uniref'] = [dbname]

                        # map all DBs to targets and add attributes if RXNs or ECs exist
                        db_targ_dict = dict(zip(orf_anno_df['ref dbname'], orf_anno_df['trim_target']))
                        if 'metacyc' in db_dict:
                            for db in db_dict['metacyc']:
                                target = db_targ_dict[db]
                                if target in rxn_dict:
                                    if 'rxn' in attrib:
                                        attrib['rxn'].extend(rxn_dict[target])
                                    else:
                                        attrib['rxn'] = rxn_dict[target]
                        if 'swissprot' in db_dict:
                            for db in db_dict['swissprot']:
                                target = db_targ_dict[db]
                                if target in ec_dict:
                                    if 'ec' in attrib:
                                        if isinstance(attrib['ec'], str):
                                            if attrib['ec'] == "":
                                                attrib['ec'] = []
                                            else:
                                                attrib['ec'] = [attrib['ec']]

                                            attrib['ec'].extend(ec_dict[target])
                                    else:
                                        attrib['ec'] = ec_dict[target]
                        if 'uniref' in db_dict:
                            for db in db_dict['uniref']:
                                target = db_targ_dict[db]
                                if target in ec_dict:
                                    if 'ec' in attrib:
                                        if isinstance(attrib['ec'], str):
                                            if attrib['ec'] == "":
                                                attrib['ec'] = []
                                            else:
                                                attrib['ec'] = [attrib['ec']]

                                            attrib['ec'].extend(ec_dict[target])
                                    else:
                                        attrib['ec'] = [ec_dict[target]]
                        if compactid in orf_to_taxonid:
                            attrib['taxon'] = orf_to_taxonid[compactid]
                        if isinstance(attrib['ec'], str):
                            if attrib['ec'] == "":
                                attrib['ec'] = []
                            else:
                                attrib['ec'] = [attrib['ec']]
                        
                        # keep all ORFs used for Ptools
                        pt_attrib_dict[shortid] = attrib
                        
                        if compactid in rep_orf_list: # only want to keep one ORF/functional annotation
                            # add ORF record to pf file
                            write_to_pf_file(output_dir_name, shortid, attrib, pfFile, compact_output=True)
                            cprod_dict[l_func].append(shortid)
                            # append to the gen elements file
                            if compact_output==False:
                                append_genetic_elements_file(genetic_elements_file, output_dir_name, shortid)
                        elif l_func in cprod_dict:
                            cprod_dict[l_func].append(shortid)
                    # Remove any duplicate ECs
                    if attrib['ec']:
                        attrib['ec'] = list(set(attrib['ec']))
                    #endfor

                elif attrib['feature'] == 'rRNA':
                    shortid = prefix + attrib['id']  #mputils.ShortenrRNAId(attrib['id'])
                    compactid = attrib['id']  #mputils.ShortenrRNAId(attrib['id'])
                    if compactid in orf_to_taxonid:
                        attrib['taxon'] = orf_to_taxonid[compactid]
                    # add ORF record to pf file
                    write_to_pf_file(output_dir_name, shortid, attrib, pfFile, compact_output=True)
                    pt_attrib_dict[shortid] = attrib

                    # append to the gen elements file
                    if compact_output==False:
                        append_genetic_elements_file(genetic_elements_file, output_dir_name, shortid)

                elif attrib['feature'] == 'tRNA':
                    shortid = prefix + attrib['id']  #mputils.ShortentRNAId(attrib['id'])
                    compactid = attrib['id']  #mputils.ShortentRNAId(attrib['id'])
                    if compactid in orf_to_taxonid:
                        attrib['taxon'] = orf_to_taxonid[compactid]
                    # add ORF record to pf file
                    write_to_pf_file(output_dir_name, shortid, attrib, pfFile, compact_output=True)
                    pt_attrib_dict[shortid] = attrib

                    # append to the gen elements file
                    if compact_output==False:
                        append_genetic_elements_file(genetic_elements_file, output_dir_name, shortid)


            '''
            #write the sequence now only once per contig
            try:
                contig_seq =  nucleotide_seq_dict[key]
            except:
                #print nucleotide_seq_dict.keys()[0]
                if countError < 10:
                    gutils.printf("ERROR: Contig %s missing file in \"preprocessed\" folder for sample\n", key)
                    countError += 1
                    if countError == 10:
                        gutils.printf("...................................................................\n")
                continue

            fastaStr=wrap("",0,62, contig_seq)

            #write_ptools_input_files(genetic_elements_file, output_dir_name, shortid, fastaStr)
            if compact_output==False:
                write_input_sequence_file(output_dir_name, shortid, fastaStr)
            '''
        #endif
    pt_attrib_df = pd.DataFrame.from_dict(pt_attrib_dict, orient='index')
    if 'ec' in pt_attrib_df.columns:
        pt_attrib_df['ec'] = ['|'.join(x) if isinstance(x, list) else '' for x in pt_attrib_df['ec']]
    if 'rxn' in pt_attrib_df.columns:
        pt_attrib_df['rxn'] = ['|'.join(x) if isinstance(x, list) else '' for x in pt_attrib_df['rxn']]
    pt_attrib_df.rename(columns={'id': 'orf_id'}, inplace=True) 
    orf = pt_attrib_df['orf_id']
    pt_attrib_df.drop(labels=['orf_id'], axis=1,inplace = True)
    pt_attrib_df.insert(0, 'orf_id', orf)
    pt_attrib_df.to_csv(anno_dir_path + "/" + sample_name + ".ptinput.tsv", sep='\t', index=False)
    # Save the complete list of ORFs that were deduped
    with open(output_dir_name + "/" + "orf_map.txt", 'w') as map_out:
        for l in cprod_dict.keys():
            map_out.write('\t'.join(cprod_dict[l]) + '\n')

    if compact_output==True:
        add_genetic_elements_file(genetic_elements_file)

    # removing annotated.gff from sample_name
    sample_name = re.sub(".annot.gff", '', sample_name)
    sample_name = re.sub('.*/', '', sample_name)
    sample_name = re.sub(r'[\\].', '', sample_name)

    # trim sample_name to less than 35 characters
    # as it causes PGDB creation to fail
    if (len(sample_name) > 35):
        sample_name = sample_name[0:35]


    if not sample_name[0].isalpha() :
        sample_name = 'E' + sample_name

    write_organisms_dat_file(output_dir_name, sample_name)

    genetic_elements_file.close()
    rename(output_dir_name + "/.tmp.genetic-elements.dat", output_dir_name + "/genetic-elements.dat")



def write_to_pf_file(output_dir_name, shortid, attrib, pfFile, compact_output):
    
    gutils.fprintf(pfFile, "ID\t%s\n", shortid)
    gutils.fprintf(pfFile, "NAME\t%s\n", shortid)

    # for orfs on the negative strand the start and the end positions
    # are reversed
    if attrib['strand'] == '+':
        gutils.fprintf(pfFile, "STARTBASE\t%s\n", attrib['start'])
        gutils.fprintf(pfFile, "ENDBASE\t%s\n", attrib['end'])
    else: # if on the negative strand
        gutils.fprintf(pfFile, "STARTBASE\t%s\n", attrib['end'])
        gutils.fprintf(pfFile, "ENDBASE\t%s\n", attrib['start'])

    if 'product' in attrib:
        gutils.fprintf(pfFile, "FUNCTION\t%s\n", attrib['product'])
    
    if 'gene-comment' in attrib:
        gutils.fprintf(pfFile, "GENE-COMMENT\t%s\n", attrib['gene-comment'])

    if 'rxn' in attrib:
        rxn_list = list(set(attrib['rxn']))
        for rxn_val in rxn_list:
            gutils.fprintf(pfFile, "METACYC\t%s\n", rxn_val)

    if 'ec' in attrib:
        ec_val = attrib['ec']
        #if ec_val:
        #    gutils.fprintf(pfFile, "EC\t%s\n", ec_val)
        ec_list = list(set(attrib['ec']))
        for ec_val in ec_list:
            gutils.fprintf(pfFile, "EC\t%s\n", ec_val)

    if 'taxon' in attrib:
        gutils.fprintf(pfFile, "TAXONOMIC-ANNOT\tTAX-%s\n", attrib['taxon'])

    if attrib['feature']=='CDS':
        gutils.fprintf(pfFile, "PRODUCT-TYPE\tP\n")

    if attrib['feature']=='tRNA':
        gutils.fprintf(pfFile, "PRODUCT-TYPE\tTRNA\n")

    if attrib['feature']=='rRNA':
        gutils.fprintf(pfFile, "PRODUCT-TYPE\trRNA\n")

    gutils.fprintf(pfFile, "//\n")


def clean_up_function_text(_product):
     _fields = [ x.strip() for x in  _product.split(' ')]
     prev_field=''
     fields = []
     for _field in _fields:
        if _field!=prev_field and _field:
            fields.append(_field)
            prev_field = _field

     return ' '.join(fields)


def write_input_sequence_file(output_dir_name, shortid, contig_sequence):

    seqfile = open(output_dir_name + "/" + shortid + ".fasta", 'w')
    gutils.fprintf(seqfile, ">%s\n%s",shortid, contig_sequence);
    seqfile.close()



def add_genetic_elements_file(genetic_elementsfile):
    gutils.fprintf(genetic_elementsfile,"ID\t0\n")
    gutils.fprintf(genetic_elementsfile,"NAME\t0\n")
    gutils.fprintf(genetic_elementsfile,"TYPE\t:CONTIG\n")
    gutils.fprintf(genetic_elementsfile,"ANNOT-FILE\t0.pf\n")
    gutils.fprintf(genetic_elementsfile,"//\n")


def append_genetic_elements_file(genetic_elementsfile, output_dir_name, shortid):
    gutils.fprintf(genetic_elementsfile,"ID\t%s\n", shortid)
    gutils.fprintf(genetic_elementsfile,"NAME\t%s\n", shortid)
    gutils.fprintf(genetic_elementsfile,"TYPE\t:CONTIG\n")
    gutils.fprintf(genetic_elementsfile,"ANNOT-FILE\t%s.pf\n", shortid)
    gutils.fprintf(genetic_elementsfile,"SEQ-FILE\t%s.fasta\n", shortid)
    gutils.fprintf(genetic_elementsfile,"//\n")


def write_organisms_dat_file(output_dir_name, sample_name):
    organism_paramsfile = open(output_dir_name + "/organism-params.dat", 'w')
    gutils.fprintf(organism_paramsfile,"ID\t%s\n",sample_name)
    gutils.fprintf(organism_paramsfile,"STORAGE\tFILE\n")
    gutils.fprintf(organism_paramsfile,"NAME\t%s\n",sample_name)
    gutils.fprintf(organism_paramsfile,"ABBREV-NAME\t%s\n",sample_name)
    gutils.fprintf(organism_paramsfile,"STRAIN\t1\n")
    gutils.fprintf(organism_paramsfile,"RANK\t|species|\n")
    gutils.fprintf(organism_paramsfile,"NCBI-TAXON-ID\t12908\n")
    organism_paramsfile.close()


def get_parameter(config_params, category, field, default = None):
     if config_params == None:
       return default

     if category in config_params:
        if field in config_params[category]:
             if config_params[category][field]:
                return config_params[category][field]
             else:
                return default
        else:
             return default
     return default


#this function creates the genbank file from the gff, protein and nucleotide sequences
def  write_gbk_file(output_file_name, contig_dict, sample_name, nucleotide_seq_dict, protein_seq_dict):
     date = sysutils.genbankDate()
     output_file_name_tmp = output_file_name + ".tmp"
     outputfile = open(output_file_name_tmp, 'w')

     count = 0
     outputStr = ""
     for key in contig_dict:
        first = True
        if count %10000 == 0:
           outputfile.write(outputStr)
           outputStr=""
        count+=1

        for attrib in contig_dict[key]:
           id  = attrib['id']
           try:
              protein_seq = protein_seq_dict[id]
           except:
              protein_seq = ""
              None

           definition = sample_name
           accession = '.'
           version = '.' +spaces(10) + "GI:."
           dblink = sample_name
           keywords = '.'
           source = sample_name
           organism = sample_name
           if first:
              first = False
              try:
                dna_seq =  nucleotide_seq_dict[key]
                dna_seq_formatted =  format_sequence_origin(dna_seq)
                dna_length = len(dna_seq)
                sourceStr = "1.." + str(dna_length)
              except:
                dna_seq = ""
                dna_seq_formatted =  ""
                dna_length = 0
                sourceStr ="0..0"

              outputStr+=("LOCUS       %-18s  %4d bp   DNA           BCT      %-11s\n" % (key, dna_length,  date))
              outputStr+=(wrap("DEFINITION  ",12,74, definition)+'\n')
              outputStr+=(wrap("ACCESSION   ", 12, 74, accession)+'\n')
              outputStr+=(wrap("VERSION     ", 12, 74, version)+'\n')
              outputStr+=(wrap("DBLINK      ", 12, 74, dblink)+'\n')
              outputStr+=(wrap("KEYWORDS    ", 12, 74,keywords)+'\n')
              outputStr+=(wrap("SOURCE    ", 12, 74, keywords)+'\n')
              outputStr+=(wrap("  ORGANISM  ",12, 74, organism)+'\n')
              outputStr+=(wrap("", 12, 74, "Metagenome")+'\n')
              outputStr+=( wrap("REFERENCE   ",12,74, "1  (bases 1 to XXXXX)")+'\n')
              outputStr+=( wrap("  AUTHORS   ",12,74, "YYYYYY,X.")+'\n')
              outputStr+=( wrap("  CONSRTM   ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("  TITLE     ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("  JOURNAL   ",12,74,"XXXXX")+'\n')
              outputStr+=( wrap("   PUBMED   ",12,74,"XXXXX")+'\n')
              outputStr+=( wrap("  REMARK   ",12,74, "XXXXX")+'\n')
              outputStr+=( wrap("COMMENT     ", 12, 74,"PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review   COMPLETENESS: XXXXX")+'\n')

              outputStr+=( wrap("FEATURES ",21,74,"Location/Qualifiers") +'\n')
              outputStr+=( wrap("     source",21,74,sourceStr) +'\n')
              outputStr+=( wrap("",21,74,"/organism=\"" + sourceStr +"\"") +'\n')
              outputStr+=( wrap("",21,74,"/strain=\"1\"")+'\n')
              outputStr+=( wrap("",21,74,"/chromosome=\"1\"") +'\n')


           if 'start' in attrib and 'end' in attrib:
               geneLoc = str(attrib['start']) +".." + str(attrib['end'])
           else:
               geneLoc = "0..0"

           if 'strand' in attrib:
              if attrib['strand']=='-':
                 geneLoc='complement' + '(' + geneLoc +')'

           outputStr+=( wrap("     gene",21,74,geneLoc) +'\n')
           if 'locus_tag' in attrib:
               locus_tag = "/locus_tag=" + "\"" + attrib['locus_tag'] + "\""
           else:
               locus_tag = "/locus_tag" + "\"\""
           outputStr+=( wrap("",21,74,locus_tag) +'\n')
           outputStr+=( wrap("     CDS",21,74, geneLoc) +'\n')
           if 'product' in attrib:
              product="/product=" + "\""+ attrib['product'] + "\""
           else:
              product="/product=\"\""
           outputStr+=( wrap("",21,74,product) +'\n')
           outputStr+=( wrap("",21,74,locus_tag) +'\n')

           codon_start="/codon_start=1"
           translation_table="/transl_table=11"
           outputStr+=( wrap("",21,74,codon_start) +'\n')
           outputStr+=( wrap("",21,74,translation_table) +'\n')

           translation= "/translation="+ protein_seq
           outputStr+=( wrap("",21,74,translation) +'\n')

        outputStr+=(wrap("ORIGIN", 21, 74, "")+'\n')
        outputStr+=(dna_seq_formatted +'\n')
        outputStr+=("//\n")

     outputfile.write(outputStr)
     outputfile.close()
     rename(output_file_name_tmp, output_file_name)


def format_sequence_origin(dna_seq):
    output=""
    Len =  len(dna_seq)
    for i in range(0, Len):
       if i==0:
          output+= '%9d' % (i+1)
       if i%10==0:
          output+=' '
       if i!=0  and i%60==0:
          output+= '\n%9d ' % (i+1)
       output +=dna_seq[i]
       i+=1
    return output


def spaces(n):
    space=''
    for  i in  range(0, n):
      space+=' '
    return space


def wrap(prefix, start, end, string):
    output=''
    prefixLen = len(prefix)
    i = prefixLen
    output+=prefix
    while i< start:
      output += ' '
      i += 1

    for c in string:
       if i > end:
          output+='\n'
          i = start
          output+=spaces(i)
       i+=1
       output+=c
    return output

def process_sequence_file(sequence_file_name,  seq_dictionary, shortorfid=False):
     sequence_file_name = gutils.correct_filename_extension(sequence_file_name)

     with gzip.open(sequence_file_name, 'rt') if sequence_file_name.endswith('.gz') \
         else open(sequence_file_name, 'r') as sequencefile:

         sequence_lines = sequencefile.readlines()
         sequencefile.close()
         fragments= []
         name=""

         seq_beg_pattern = re.compile(">(\S+)")
         for line in sequence_lines:
            line = line.strip()
            res = seq_beg_pattern.search(line)
            if res:
              if len(name) > 0:
                 sequence=''.join(fragments)
                 seq_dictionary[name]=sequence
                 fragments = []
              if shortorfid:
                 name=get_sequence_number(line)
              else:
                 name=res.group(1)
            else:
              fragments.append(line)

         # add the final sequence
         if len(name) > 0:
             sequence=''.join(fragments)
             seq_dictionary[name]=sequence

#        if count > 1000:
#           sys.exit(0)
#        count=count+1
#     fields = re.split('\t', line)
        #print table[str(fields[0].strip())]
     #print blast_file + ' ' + tax_maps + ' ' + database


usage =  sys.argv[0] + """ -g gff_files -n nucleotide_sequences -p protein_sequences [--out-gbk gbkfile --out-ptinput ptinputdir]\n"""
def read_taxons_for_orfs(ncbi_taxonomy_tree, taxonomy_table):
    num = 20
    if not path.exists(ncbi_taxonomy_tree):
        ncbi_taxonomy_tree = ncbi_taxonomy_tree + ".gz"

    name_to_taxonid = {}
    ncbi_taxonomy_tree = gutils.correct_filename_extension(ncbi_taxonomy_tree)
    with gzip.open(ncbi_taxonomy_tree, 'rt') if ncbi_taxonomy_tree.endswith(".gz") \
        else open(ncbi_taxonomy_tree, 'r') as f:

        for _line in f:
           fields = [ x.strip() for x in _line.split('\t') ]
           if len(fields) == 2:
              name_to_taxon[fields[0]] = fields[1]


    tax_PATT = re.compile(r'(\d+)')
    orf_to_taxonid = {}
    tax_col = 0

    taxonomy_table = gutils.correct_filename_extension(taxonomy_table)
    with gzip.open(taxonomy_table, 'rt') if taxonomy_table.endswith('.gz') \
        else open(taxonomy_table, 'r') as ftaxin:
        fields = [ x.strip() for x in ftaxin.readline().strip().split('\t') ]
        for i in range(0, len(fields)):
           if fields[i]=='taxonomy':
              tax_col =i

        for _line in ftaxin:
           fields = [ x.strip() for x in _line.split('\t') ]
           if len(fields) > tax_col:
              res = tax_PATT.search(fields[tax_col])
              if res:
                  orf_to_taxonid[fields[0]] = res.group(1)

    return orf_to_taxonid

def main(argv, errorlogger = None, runstatslogger = None):
    # Parse options (many!)
    # TODO: Create option groups
    # filtering options
    parser = createParser()
    options, args = parser.parse_args(argv)

    if not(options.gff_file or options.nucleotide_sequences or options.protein_sequences or options.output):
      sys.exit(0)


    if not options.gff_file:
       gutils.eprintf("ERROR\tGFF file not specified\n")
       errorlogger.printf("ERROR\tGFF file not specified\n")


    if not options.gbk_file and not options.ptinput_file:
       gutils.eprintf("ERROR:No genbank or ptools input is specified\n")
       errormod.insert_error(errorcode)
       return (1,'')


    if not path.exists(options.gff_file):
        print("gff file does not exist")
        gutils.eprintf("ERROR\tGFF file %s  not found\n", options.gff_file)
        errorlogger.printf("ERROR\tGFF file %s  not found\n", options.gff_file)
        sys.exit(0)

    output_files = {}
    input_files = {}
    if  options.gbk_file:
       output_files['gbk'] = options.gbk_file

    if  options.ptinput_file:
       output_files['ptinput'] = options.ptinput_file
       output_files['ec_mapping'] = options.ec_mapping
       output_files['rxn_mapping'] = options.rxn_mapping
       output_files['annotation_table'] = options.annotation_table
       output_files['taxonomy_table'] = options.taxonomy_table

    nucleotide_seq_dict = {}
    protein_seq_dict = {}
    if options.nucleotide_sequences  and path.exists(options.nucleotide_sequences):
       process_sequence_file(options.nucleotide_sequences, nucleotide_seq_dict)

    if gutils.plain_or_gz_file_exists(options.protein_sequences):
       process_sequence_file(options.protein_sequences, protein_seq_dict)

    orf_to_taxonid={}
    if options.ncbi_taxonomy_tree != None and options.taxonomy_table != None:
      orf_to_taxonid = read_taxons_for_orfs(options.ncbi_taxonomy_tree, options.taxonomy_table)

    process_gff_file(options.gff_file, output_files, nucleotide_seq_dict, \
        protein_seq_dict, input_files,  orf_to_taxonid=orf_to_taxonid, \
        compact_output=options.compact_output)

    sample_name = get_sample_name(options.gff_file)

    if  options.ptinput_file:
      gutils.createDummyFile(options.ptinput_file + PATHDELIM + sample_name + ".dummy.txt")

def MetaPathways_create_genbank_ptinput(argv, errorlogger = None, runstatslogger = None):
    global errorcode
    try:
       main(argv, errorlogger = errorlogger, runstatslogger = runstatslogger)
    except:
       errormod.insert_error(12)
       return (1, traceback.print_exc(10))

    return (0,'')

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main(sys.argv[1:])


