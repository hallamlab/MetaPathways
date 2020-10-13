"""Creates the jobs (corresponding to the steps of the pipeline) that
   executed the by execution managment code
"""

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2020, MetaPathways"
__version__ = "3.5.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
    import traceback
    import os
    import re
    import sys
    import shutil

    from shutil import rmtree
    from optparse import make_option
    from os import path, _exit

    from metapathways import sysutil as sysutils
    from metapathways import general_utils as gutils
    from metapathways import context as contextmod
except:
    print(""" Could not load some user defined  module functions""")
    print(traceback.print_exc(10))
    sys.exit(3)


PATHDELIM = sysutils.pathDelim()

@gutils.Singleton
class Params:
    params  = {}
    def __init__(self, params):
        self.params = params
        pass

    def get(self, key1, key2 = None, default = None):
        if not key1 in self.params:
            return default

        if key2 == None:
            return self.params[key1]

        if not key2 in self.params[key1]:
            return default

        return self.params[key1][key2]

    def print_key_values(self):
        print(self.params)


@gutils.Singleton
class Configs:
    def __init__(self, configs):
        for key, value in configs.items():
            setattr(self, key, value)

@gutils.Singleton
class ContextCreator:
    params = None
    configs = None
    factory = {}
    stageList = {}
    format = None

    def _Message(self, str):
        return '{0: <60}'.format(str)

    def create_quality_check_cmd(self, s):
        """ PREPROCESS_INPUT """
        contexts = []

        '''inputs'''
        input_file = s.input_file


        '''outputs'''
        output_fas = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        mapping_file = s.preprocessed_dir + PATHDELIM + s.sample_name + ".mapping.txt"
        nuc_stats_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".nuc.stats"
        contig_lengths_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".contig.lengths.txt"

        '''params'''
        min_length = self.params.get('quality_control','min_length', default = '180')
        type = 'nucleotide'

        context = contextmod.Context()
        context.name = 'PREPROCESS_INPUT'
        context.inputs = { 'input_file' : input_file }

        context.outputs = {
                            'output_fas': output_fas, 'mapping_file': mapping_file,\
                            'nuc_stats_file' : nuc_stats_file,\
                            'contig_lengths_file' : contig_lengths_file,\
                            'mapping_file': mapping_file
                          }


        context.status = self.params.get('metapaths_steps','PREPROCESS_INPUT')


        pyScript = self.configs.PREPROCESS_INPUT

        cmd = "%s --min_length %d --log_file %s  -i %s -o  %s -M %s -t %s -L %s"\
             %( pyScript, float(min_length), context.outputs['nuc_stats_file'], context.inputs['input_file'],
                context.outputs['output_fas'], context.outputs['mapping_file'],\
                 type, context.outputs['contig_lengths_file']\
               )

        context.message = self._Message("PREPROCESSING THE INPUT")
        context.commands = [cmd]
        contexts.append(context)
        return contexts


    def create_preprocess_input_aminos_cmd(self, s):
        """ PREPROCESS_AMINOS """
        contexts = []

        '''inputs'''
        input_file = s.input_file


        '''outputs'''
        output_fasta = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        mapping_file =  s.preprocessed_dir + PATHDELIM + s.sample_name + ".mapping.txt"
        nuc_stats_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".nuc.stats"
        contig_lengths_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".contig.lengths.txt"

        output_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".faa"
        output_fna = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".fna"
        output_gff = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".unannot.gff"

        '''params'''

        context = contextmod.Context()
        context.name = 'PREPROCESS_AMINOS'
        context.inputs = { 'input_file' : input_file }

        context.outputs = {
                            'output_fasta': output_fasta, 'output_faa': output_faa,\
                            'output_fna': output_fna, 'output_gff': output_gff,\
                            'mapping_file': mapping_file,\
                            'nuc_stats_file' : nuc_stats_file,\
                            'contig_lengths_file' : contig_lengths_file\
                          }

        context.status = self.params.get('metapaths_steps','PREPROCESS_INPUT')

        pyScript = self.configs.PREPROCESS_AMINOS

        cmd = "%s --log_file %s  -i %s  -M %s  -L %s --fasta %s --faa %s --fna %s --gff %s"\
             %( pyScript,  context.outputs['nuc_stats_file'],\
                 context.inputs['input_file'], context.outputs['mapping_file'],\
                 context.outputs['contig_lengths_file'],\
                 context.outputs['output_fasta'], context.outputs['output_faa'],\
                 context.outputs['output_fna'],context.outputs['output_gff']
               )

        context.message = self._Message("PREPROCESSING THE AMINOS INPUT")
        context.commands = [cmd]
        contexts.append(context)
        return contexts


    def  convert_gbk_to_fna_faa_gff_annotated(self, s):
        contexts = self._convert_gbk_to_fna_faa_gff(s, annotated = True,\
                                                    create_functional_table = True)
        return contexts

    def  convert_gbk_to_fna_faa_gff_unannotated(self, s):
        contexts = self._convert_gbk_to_fna_faa_gff(s, annotated = False)
        return contexts

    def  _convert_gbk_to_fna_faa_gff(self, s, annotated = False, create_functional_table = False):
        """ GBK_TO_FNA_FAA_GFF """
        contexts = []

        '''inputs'''
        input_gbk = s.input_file

        '''outputs'''
        output_fna = s.preprocessed_dir + s.sample_name + ".fasta"
        output_faa = s.orf_prediction_dir + s.sample_name + ".faa"
        output_annot_table = s.output_results_annotation_table_dir +\
                             PATHDELIM + 'functional_and_taxonomic_table.txt'
        contig_lengths_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".contig.lengths.txt"


        output_gff = ''
        if annotated == True:
            output_gff = s.genbank_dir + s.sample_name + ".annot.gff"
        else:
            output_gff = s.orf_prediction_dir + s.sample_name + ".unannot.gff"

        mapping_file =  s.preprocessed_dir + PATHDELIM + s.sample_name + ".mapping.txt"

        context = contextmod.Context()
        context.inputs = { 'input_gbk' : input_gbk }
        context.outputs = {
                             'output_gff' : output_gff,
                             'output_fna':output_fna,
                             'output_faa':output_faa,
                             'mapping_file': mapping_file,
                             'output_annot_table':output_annot_table,
                             'contig_lengths_file':contig_lengths_file
                          }


        context.name = 'PREPROCESS_INPUT'
        context.status = self.params.get('metapaths_steps','PREPROCESS_INPUT')
        pyScript = self.configs.GBK_TO_FNA_FAA_GFF


        cmd = "%s -g %s --output-fna %s --output-faa %s --output-gff %s -M %s -L %s "\
               %( pyScript, context.inputs['input_gbk'], context.outputs['output_fna'],\
                 context.outputs['output_faa'], context.outputs['output_gff'],\
                 context.outputs['mapping_file'],\
                 context.outputs['contig_lengths_file']\
                )

        if create_functional_table :
            cmd = cmd + " --create-functional-table %s" %(context.outputs['output_annot_table'])

        context.message = self._Message("PREPROCESSING THE GBK INPUT")
        context.commands = [cmd]

        contexts.append(context)

        return contexts



    def create_orf_prediction_cmd(self, s) :
        """ ORF_PREDICTION """
        contexts = []

        '''inputs'''
        input_file = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"

        '''outputs'''
        output_gff = s.orf_prediction_dir + s.sample_name + ".gff"

        context = contextmod.Context()
        context.name = 'ORF_PREDICTION'
        context.inputs = { 'input_file' : input_file }
        context.outputs = { 'output_gff' : output_gff }
        context.status = self.params.get('metapaths_steps','ORF_PREDICTION')
        translation_table = self.params.get('orf_prediction', 'translation_table')
        algorithm = self.params.get('orf_prediction', 'algorithm')
        strand = self.params.get('orf_prediction', 'strand')

        mode = self.params.get('orf_prediction', 'mode')

        pyScript = self.configs.ORF_PREDICTION

        if algorithm == "prodigal":
            executable =  self.configs.PRODIGAL_EXECUTABLE

        if algorithm == "FGS+":
            executable = self.configs.FGSPlus_EXECUTABLE


        cmd = [
                  pyScript,
                  "--prod_exec", executable,
                  "--prod_m",
                  "--prod_p", mode,
                  "--prod_f", "gff",
                  "--prod_g", translation_table,
                  "--prod_input", context.inputs['input_file'],
                  "--prod_output", context.outputs['output_gff'], #"--strand",  strand
             ]
        if algorithm == "FGS+":
            num_threads = self.configs.NUM_CPUS
            cmd += ["--algorithm", "FGS+"]
            cmd += ["--nthreads", num_threads]

        context.commands = [' '.join(cmd)]
        context.message = self._Message("ORF PREDICTION")
        contexts.append(context)
        return contexts

    def create_aa_orf_sequences_cmd(self, s):
        """ ORF_TO_AMINO """
        contexts = []

        ''' inputs '''
        input_gff = s.orf_prediction_dir + s.sample_name + ".gff"
        input_fasta = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"

        '''outputs'''
        output_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".faa"
        output_fna = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".fna"
        output_gff = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".unannot.gff"

        context = contextmod.Context()
        context.name = 'ORF_TO_AMINO'
        context.inputs = { 'input_gff' : input_gff, 'input_fasta': input_fasta }
        context.outputs = { 'output_faa': output_faa, 'output_fna': output_fna, 'output_gff' : output_gff }
        context.status = self.params.get('metapaths_steps','ORF_PREDICTION')

        pyScript = self.configs.ORF_TO_AMINO
        cmd = "%s -g  %s  -n %s --output_nuc %s --output_amino %s --output_gff %s"\
               %(pyScript, context.inputs['input_gff'], context.inputs['input_fasta'],\
                 context.outputs['output_fna'], context.outputs['output_faa'],\
                 context.outputs['output_gff'])

        context.message = self._Message("CREATING AMINO ACID SEQS FROM GFF FILE")
        context.commands = [cmd]
        contexts.append(context)
        return contexts

    def create_create_filtered_amino_acid_sequences_cmd(self, s):
        """FILTER_AMINOS"""
        contexts = []

        '''inputs'''
        input_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".faa"


        '''outputs'''
        output_filtered_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".qced.faa"
        amino_stats_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".amino.stats"
        orf_lengths_file = s.output_run_statistics_dir + PATHDELIM + s.sample_name + ".orf.lengths.txt"

        '''params'''
        min_length = self.params.get('orf_prediction', 'min_length', default = 60)
        type = 'amino'

        context = contextmod.Context()
        context.name = 'FILTER_AMINOS'
        context.inputs = { 'input_faa' : input_faa }
        context.outputs = { 'output_filtered_faa': output_filtered_faa,\
                            'amino_stats_file': amino_stats_file,\
                            'orf_lengths_file': orf_lengths_file }


        context.status = self.params.get('metapaths_steps','FILTER_AMINOS')

        pyScript = self.configs.PREPROCESS_INPUT

        cmd = "%s  --min_length %s --log_file %s  -i %s -o  %s -L %s -t %s"\
              %( pyScript, min_length, context.outputs['amino_stats_file'],\
               context.inputs['input_faa'],  context.outputs['output_filtered_faa'],\
               context.outputs['orf_lengths_file'], type)

        context.commands = [cmd]
        context.message = self._Message("FILTER AMINO ACID SEQS")
        contexts.append(context)
        return contexts



    def create_refscores_compute_cmd(self, s):
        """COMPUTE_REFSCORES"""
        contexts = []

        '''inputs'''
        input_filtered_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".qced.faa"

        '''outputs'''
        output_refscores =  s.blast_results_dir + PATHDELIM + s.sample_name + ".refscores" + "." + s.algorithm

        context = contextmod.Context()
        context.name = 'COMPUTE_REFSCORES'
        context.inputs = { 'input_filtered_faa' : input_filtered_faa }
        context.outputs = { 'output_refscores': output_refscores}

        context.status = self.params.get('metapaths_steps','PARSE_FUNC_SEARCH')

        cmd = None
        if s.algorithm == 'BLAST':
            pyScript = self.configs.COMPUTE_REFSCORES
            cmd = "%s   -o %s -i %s -a  %s"\
                 %( pyScript,  output_refscores, input_filtered_faa, s.algorithm)

        elif s.algorithm == 'LAST':
            pyScript = self.configs.COMPUTE_REFSCORES
            cmd = "%s   -o %s -i %s -a %s"\
                     %( pyScript,  output_refscores, input_filtered_faa, s.algorithm)

        context.commands = [cmd]
        context.message = self._Message("COMPUTING REFSCORES FOR BITSCORE")
        contexts.append(context)
        return contexts

    def get_dbstring(self) :
        dbstring = self.params.get('annotation', 'dbs', default = '')
        return dbstring

    def create_blastp_against_refdb_cmd(self, s):
        """FUNC_SEARCH"""
        contexts = []

        '''parameters'''

        max_evalue = self.params.get('annotation', 'max_evalue', default = 0.000001)
        max_hits = self.params.get('annotation', 'max_hits', default = 5)
        min_score = self.params.get('annotation', 'min_score', default = 20)

        num_threads = self.configs.NUM_CPUS

        dbstring = self.get_dbstring()
        dbs = [x.strip() for x in dbstring.split(",") if len(x) != 0]

        for db in dbs:
            '''inputs'''
            input_filtered_faa = s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".qced.faa"

            '''outputs'''
            blastoutput = s.blast_results_dir + PATHDELIM + s.sample_name + "." + db + "." + s.algorithm + "out"

            refDbFullName = self.configs.REFDBS + PATHDELIM + 'functional'+ PATHDELIM +\
                            'formatted' + PATHDELIM + db

            context = contextmod.Context()
            context.name = 'FUNC_SEARCH:' +db
            context.inputs = { 'input_filtered_faa' : input_filtered_faa }
            context.outputs = { 'blastoutput': blastoutput}

            cmd = None
            if s.algorithm == 'BLAST':
                pyScript =  self.configs.FUNC_SEARCH

                searchExec = self.configs.BLASTP_EXECUTABLE
                if not path.exists(self.configs.BLASTP_EXECUTABLE):
                    searchExec =  self.configs.BLASTP_EXECUTABLE

                cmd = ("%s "
                      "--algorithm %s "
                      "--blast_executable %s "
                      "--num_threads %s "
                      "--blast_max_target_seqs %s "
                      "--blast_outfmt 6 "
                      "--blast_db %s "
                      "--blast_query  %s "
                      "--blast_evalue  %s "
                      "--blast_out %s")\
                     %(pyScript,
                       s.algorithm,
                       searchExec,
                       num_threads,
                       str(max_hits),
                       refDbFullName,
                       input_filtered_faa,
                       str(max_evalue), blastoutput)

                context.message = self._Message("BLASTING AMINO SEQS AGAINST " + db)

            if s.algorithm == 'LAST':
                pyScript =  self.configs.FUNC_SEARCH
                searchExec = self.configs.LAST_EXECUTABLE
                cmd = ("%s "
                    "--algorithm %s "
                    "--last_executable %s "
                    "--num_threads %s "
                    "--last_o %s "
                    "--last_f 2 "
                    "--last_db %s "
                    "--last_query %s") \
                    %(pyScript,
                    s.algorithm ,
                    searchExec,
                    num_threads,
                    blastoutput,
                    refDbFullName,
                    input_filtered_faa)

                context.message = self._Message("LASTING AMINO SEQS AGAINST " + db)

            context.status = self.params.get('metapaths_steps','FUNC_SEARCH')
            context.commands = [cmd]

            contexts.append(context)

        return contexts


    def create_parse_blast_cmd(self, s ):
        """  Command for parsing the blast flie snd create the parse blast files
          input -- blastoutput
          output -- parseed files
          refscorefile   -- refscore file
          min_bsr   -- minimum bsr ratio for accepting into annotation
          max_evalue  -- max evalue
          min_score   -- min score
          min_length  -- min_length in amino acids, typcically 100 amino acids.ould be minimum
        """
        """PARSE_FUNC_SEARCH"""
        contexts = []

        '''parameters'''
        min_bsr = self.params.get('annotation', 'min_bsr', default = 0.4)
        min_score = self.params.get('annotation', 'min_score', default = 0.0)
        min_length = self.params.get('annotation', 'min_length', default = 100)
        max_evalue = self.params.get('annotation', 'max_evalue', default = 1000)


        dbstring = self.get_dbstring()
        dbs = [x.strip() for x in dbstring.split(",")  if len(x) != 0]

        pyScript =  self.configs.PARSE_FUNC_SEARCH
        for db in dbs:
            '''inputs'''
            input_db_blastout = s.blast_results_dir + PATHDELIM + s.sample_name + "." + db + "." + s.algorithm+"out"
            refscorefile = s.blast_results_dir + PATHDELIM + s.sample_name +\
                           ".refscores" +"." + s.algorithm
            dbmapFile = self.configs.REFDBS + PATHDELIM + 'functional' + PATHDELIM + 'formatted'\
                        + PATHDELIM + db + "-names.txt"

            '''outputs'''
            output_db_blast_parse = s.blast_results_dir + PATHDELIM + s.sample_name +\
                                    "." + db + "." + s.algorithm+"out.parsed.txt"

            context = contextmod.Context()
            context.name = 'PARSE_FUNC_SEARCH:' + db
            context.inputs = { 'input_db_blastout' : input_db_blastout,\
                               'dbmapFile': dbmapFile,\
                               'refscorefile': refscorefile
                              }

            context.outputs = { 'output_db_blast_parse':output_db_blast_parse}

            cmd = "%s -d %s  -b %s -m %s  -r  %s  --min_bsr %s  --min_score %s --min_length %s --max_evalue %s"\
                  %( pyScript, db, context.inputs['input_db_blastout'],\
                  context.inputs['dbmapFile'],  context.inputs['refscorefile'],\
                  min_bsr, min_score, min_length, max_evalue)

            if s.algorithm == 'LAST':
                cmd = cmd + ' --algorithm LAST'

            if s.algorithm == 'BLAST':
                cmd = cmd + ' --algorithm BLAST'
            context.commands = [cmd]
            context.status = self.params.get('metapaths_steps','PARSE_FUNC_SEARCH')
            context.message = self._Message("PARSING " + s.algorithm + " OUTPUT FOR " + db)
            contexts.append(context)

        return contexts


    def create_scan_rRNA_seqs_cmd(self, s):
        """SCAN_rRNA"""
        contexts = []

        '''parameters'''
        bscore_cutoff = self.params.get('rRNA', 'min_bitscore', default = 27)
        eval_cutoff = self.params.get( 'rRNA', 'max_evalue', default = 6)
        identity_cutoff = self.params.get('rRNA', 'min_identity', default = 40)
        dbstring = self.params.get('rRNA', 'refdbs', default = None)
        refrRNArefDBs = [x.strip() for x in dbstring.split(',') if len(x.strip())]

        pyScript = self.configs.PARSE_FUNC_SEARCH

        '''inputs'''
        input_fasta = s.preprocessed_dir +  PATHDELIM + s.sample_name + ".fasta"

        pyScript = self.configs.SCAN_rRNA

        algorithm  = s.algorithm.upper()
        for db in refrRNArefDBs:
            '''inputs'''
            dbpath = self.configs.REFDBS + PATHDELIM + 'taxonomic' + PATHDELIM + 'formatted' + PATHDELIM + db
            dbsequences = self.configs.REFDBS + PATHDELIM + "taxonomic" + PATHDELIM+  db

            '''outputs'''
            rRNA_blastout = s.blast_results_dir + PATHDELIM + s.sample_name + ".rRNA." + db + "." + algorithm + "out"
            rRNA_stat_results = s.output_results_rRNA_dir + s.sample_name + "." + db + ".rRNA.stats.txt"

            context = contextmod.Context()
            context.name = 'SCAN_rRNA:' + db
            context.inputs = {  'input_fasta':input_fasta, 'dbsequences':dbsequences }
            context.inputs1 = { 'dbpath' : dbpath }
            context.outputs = { 'rRNA_blastout':rRNA_blastout, 'rRNA_stat_results': rRNA_stat_results }

            """ the BLAST part first """
            blast_cmd = ""
            if True or algorithm == "BLAST":
                executable = shutil.which('blastn')
                if executable == None:
                    eprintf("ERROR\tCannot find blastn\n")
                        #logger.printf("ERROR\tCannot find blastn to format\n")

                blast_cmd = "%s -outfmt 6 -num_threads 8  -query %s -out %s -db %s -max_target_seqs 5"\
                      %(executable, context.inputs['input_fasta'], context.outputs['rRNA_blastout'], context.inputs1['dbpath'])

            if False and algorithm == "LAST":
                executable = shutil.which(self.configs.LAST_EXECUTABLE)
                if executable == None:
                    eprintf("ERROR\tCannot find blastn\n")
                        #logger.printf("ERROR\tCannot find blastn to format\n")
                blast_cmd = "%s -f 2 -o %s %s %s"\
                     %(executable, context.outputs['rRNA_blastout'], context.inputs1['dbpath'], context.inputs['input_fasta'])


            """ now the scanning part"""
            scan_cmd = "%s -o %s -b %s -e %s -s %s"  %(pyScript, context.outputs['rRNA_stat_results'],\
                  bscore_cutoff, eval_cutoff, identity_cutoff)

            scan_cmd = scan_cmd +  " -i "  + context.outputs['rRNA_blastout'] + " -d " + context.inputs['dbsequences']
            context.commands = [scan_cmd, blast_cmd]
            context.status = self.params.get('metapaths_steps','SCAN_rRNA')
            context.message = self._Message("SCANNING FOR rRNA USING DB " + db)
            contexts.append(context)

        return contexts

    def create_tRNA_scan_statistics(self, s):
        """SCAN_tRNA"""

        contexts = []

        '''inputs'''
        input_fasta = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        TPCsignal = self.configs.RESOURCES_DIR + PATHDELIM + 'TPCsignal'
        Dsignal = self.configs.RESOURCES_DIR+ PATHDELIM + 'Dsignal'

        '''outputs'''
        tRNA_stats_output = s.output_results_tRNA_dir + PATHDELIM + s.sample_name +  ".tRNA.stats.txt"
        tRNA_fasta_output = s.output_results_tRNA_dir + PATHDELIM + s.sample_name +  ".tRNA.fasta"


        context = contextmod.Context()
        context.name = 'SCAN_tRNA'
        context.inputs = { 'input_fasta':input_fasta, 'TPCsignal':TPCsignal, 'Dsignal':Dsignal }
        context.outputs = { 'tRNA_stats_output':tRNA_stats_output, 'tRNA_fasta_output': tRNA_fasta_output}


        pyScript = self.configs.SCAN_tRNA
        executable = self.configs.SCAN_tRNA_EXECUTABLE
        cmd = "%s --executable %s -o %s -F %s  -i %s -T %s  -D %s"\
             %(pyScript, executable, context.outputs['tRNA_stats_output'], context.outputs['tRNA_fasta_output'],\
             context.inputs['input_fasta'], context.inputs['TPCsignal'], context.inputs['Dsignal'])

        context.commands = [cmd]
        context.status = self.params.get('metapaths_steps','SCAN_tRNA')
        context.message = self._Message("SCANNING FOR tRNA USING tRNA-Scan")
        contexts.append(context)
        return contexts


    def create_annotate_genebank_cmd(self, s ):
        """ANNOTATE ORFS"""
        contexts = []

        '''inputs'''
        input_unannotated_gff = s.orf_prediction_dir + PATHDELIM + s.sample_name+".unannot.gff"
        mapping_txt =  s.preprocessed_dir + PATHDELIM + s.sample_name + ".mapping.txt"

        '''outputs'''
        output_annotated_gff  = s.genbank_dir + PATHDELIM + s.sample_name+".annot.gff"
        output_comparative_annotation  =  s.output_results_annotation_table_dir \
                                            + PATHDELIM + s.sample_name
        dbstring = self.get_dbstring()
        refdbs = [x.strip() for x in dbstring.split(",")  if len(x) != 0]

        rRNAdbstring = self.params.get('rRNA', 'refdbs', default = None)
        rRNAdbs = [x.strip() for x in rRNAdbstring.split(",")  if len(x) != 0]

        context = contextmod.Context()
        context.name = 'ANNOTATE_ORFS'

        context.inputs = {
            'input_unannotated_gff':input_unannotated_gff
        }

        context.inputs1 = {
            'mapping_txt':mapping_txt,
        }
        context.outputs = {
           'output_annotated_gff':output_annotated_gff,
        }
        context.outputs1 = {
            'output_comparative_annotation':output_comparative_annotation
        }

        context.status = self.params.get('metapaths_steps','ANNOTATE_ORFS')

        '''use rRNA stats if they are available'''
        options = ''
        for rRNArefdb in rRNAdbs:
            rRNA_stat_results = s.output_results_rRNA_dir + s.sample_name +\
                               '.' + rRNArefdb + '.rRNA.stats.txt'
            #print rRNA_stat_results
            if gutils.hasResults(rRNA_stat_results)  :
                context.inputs['rRNA_stat_results']  = rRNA_stat_results
                options += " --rRNA_16S " +  context.inputs['rRNA_stat_results']


        '''use rRNA stats if they are available'''
        tRNA_stat_results = s.output_results_tRNA_dir + PATHDELIM + s.sample_name + '.tRNA.stats.txt'
        if gutils.hasResults(tRNA_stat_results):
            context.inputs['tRNA_stat_results']  = tRNA_stat_results
            options += " --tRNA " +  context.inputs['tRNA_stat_results']


        pyScript = self.configs.ANNOTATE_ORFS
        cmd = "%s --input_gff  %s -o %s  %s --output-comparative-annotation %s \
                  --algorithm %s "\
              %(pyScript, context.inputs['input_unannotated_gff'],\
              context.outputs['output_annotated_gff'],  options,\
              context.outputs1['output_comparative_annotation'],s.algorithm )




        for refdb in refdbs:
            parsed_file =  s.blast_results_dir + PATHDELIM + s.sample_name\
                            + "." + refdb+ "." + s.algorithm + "out.parsed.txt"
            context.inputs[parsed_file] = parsed_file
#               cmd = cmd + " -b " + parsed_file + " -d " + refdb + " -w 1 "


        cmd = cmd + " -m " + context.inputs1['mapping_txt']
        cmd = cmd + " -D " + s.blast_results_dir + " -s " + s.sample_name

        context.message = self._Message("ANNOTATE ORFS")
        context.commands = [cmd]
        contexts.append(context)
        return contexts

    def create_genbank_file_cmd(self, s):
        """GENBANK_FILE"""

        contexts = []

        '''inputs'''
        input_annot_gff  = s.genbank_dir +PATHDELIM + s.sample_name+".annot.gff"
        input_nucleotide_fasta = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        input_amino_acid_fasta =  s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".qced.faa"

        '''outputs'''
        output_annot_gbk = s.genbank_dir + PATHDELIM + s.sample_name +  '.gbk'

        context = contextmod.Context()
        context.name = 'GENBANK_FILE'

        context.inputs = {
            'input_annot_gff':input_annot_gff,
            'input_nucleotide_fasta':input_nucleotide_fasta,
            'input_amino_acid_fasta':input_amino_acid_fasta
        }
        context.outputs = {
            'output_annot_gbk': output_annot_gbk
        }

        pyScript = self.configs.GENBANK_FILE
        cmd = "%s -g %s -n %s -p %s " %(pyScript, context.inputs['input_annot_gff'],\
             context.inputs['input_nucleotide_fasta'], context.inputs['input_amino_acid_fasta'])

        """GENBANK_FILE"""
        genbank_file_status = self.params.get('metapaths_steps','GENBANK_FILE')
        if genbank_file_status in ['redo'] or\
           (genbank_file_status in ['yes'] and not s.hasGenbankFile() ):
            cmd += ' --out-gbk ' + context.outputs['output_annot_gbk']
        context.message =  self._Message("GENBANK FILE" )

        context.status =  genbank_file_status = self.params.get('metapaths_steps','GENBANK_FILE')
        context.commands = [cmd]
        contexts.append(context)
        return contexts

    def create_ptinput_cmd(self, s):
        """PATHOLOGIC_INPUT"""

        contexts = []

        '''inputs'''
        input_annot_gff  = s.genbank_dir +PATHDELIM + s.sample_name+".annot.gff"
        input_nucleotide_fasta = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        input_amino_acid_fasta =  s.orf_prediction_dir + PATHDELIM +  s.sample_name + ".qced.faa"

        basencbi = self.configs.REFDBS + PATHDELIM + 'ncbi_tree'
        ncbi_tree = basencbi + PATHDELIM + 'ncbi_taxonomy_tree.txt'
        taxonomy_table = s.output_results_annotation_table_dir + PATHDELIM + s.sample_name + '.functional_and_taxonomic_table.txt'
        '''outputs'''

        context = contextmod.Context()
        context1 = contextmod.Context()
        context.name = 'PATHOLOGIC_INPUT'

        context.inputs = {
                           'input_annot_gff':input_annot_gff,
                           'input_amino_acid_fasta':input_amino_acid_fasta
                         }

        context.inputs1 = {
                             'ncbi_tree': ncbi_tree,
                             'taxonomy_table': taxonomy_table
                          }


        context.inputs_optional = {
                           'input_nucleotide_fasta':input_nucleotide_fasta
                         }


        #  'output_fasta_pf_dir_fasta':s.output_fasta_pf_dir + PATHDELIM +  '0.fasta',
        input_output_gbk  = s.genbank_dir +PATHDELIM + s.sample_name+".gbk"

        """ If we include it this folder then it is deleted in the execution.py during the
            <context>.removeOutput call
            #'output_fasta_pf_dir':s.output_fasta_pf_dir,
        """
        #'output_fasta_pf_dir':s.output_fasta_pf_dir,
        context.outputs = {
          'output_fasta_pf_dir_genetic':s.output_fasta_pf_dir + PATHDELIM + 'genetic-elements.dat',
          'output_fasta_pf_dir_organism':s.output_fasta_pf_dir + PATHDELIM +  'organism-params.dat',
          'dummy_ouptut_file':s.output_fasta_pf_dir + PATHDELIM +  s.sample_name + '.dummy.txt',
          'output_annot_gbk':s.genbank_dir + PATHDELIM + s.sample_name +  '.gbk'
        }

        context1.outputs = {
          'output_fasta_pf_dir_pf':s.output_fasta_pf_dir + PATHDELIM +  '0.pf',
        }



        pyScript = self.configs.GENBANK_FILE
        ''' remove the fasta files '''
        #cmd="%s -g %s -n %s -p %s " %(pyScript, context.inputs['input_annot_gff'],\
        #     context.inputs_optional['input_nucleotide_fasta'], context.inputs_optional['input_amino_acid_fasta'])

        cmd = "%s -g %s " %(pyScript, context.inputs['input_annot_gff'])

        #"""PATHOLOGIC_INPUT"""
        context.status = self.params.get('metapaths_steps','PATHOLOGIC_INPUT')

           
        if context.status in ['redo'] or (context.status in ['yes'] and not s.hasPToolsInput() ):
            cmd += ' --out-ptinput ' + s.output_fasta_pf_dir
            cmd += ' -n ' + context.inputs_optional['input_nucleotide_fasta']
            cmd += ' --ncbi-tree ' + context.inputs1['ncbi_tree']
            cmd += ' --taxonomy-table ' + context.inputs1['taxonomy_table']
            cmd += ' -p ' + context.inputs['input_amino_acid_fasta']
            cmd += ' --out-gbk ' + context.outputs['output_annot_gbk']

        context.message = self._Message("PATHOLOGIC INPUT" )

        context.commands = [cmd]
        contexts.append(context)
        return contexts


    def create_report_files_cmd(self, s):
        """CREATE_ANNOT_REPORTS"""

        contexts = []
        '''input'''
        input_annot_gff  = s.genbank_dir +PATHDELIM + s.sample_name+ ".annot.gff"

        '''output'''
        output_annot_table = s.output_results_annotation_table_dir +\
                             PATHDELIM + s.sample_name + '.functional_and_taxonomic_table.txt'

        context = contextmod.Context()
        context.name = 'CREATE_ANNOT_REPORTS'
        basefun  =  self.configs.REFDBS + PATHDELIM + 'functional_categories'
        basencbi = self.configs.REFDBS + PATHDELIM + 'ncbi_tree'
        context.inputs = {
            'input_annot_gff':input_annot_gff,
            'KO_classification':basefun + PATHDELIM +  'KO_classification.txt',
            'COG_categories':basefun + PATHDELIM +  'COG_categories.txt',
            'SEED_subsystems':basefun + PATHDELIM + 'SEED_subsystems.txt',
            'CAZY_hierarchy':basefun + PATHDELIM + 'CAZY_hierarchy.txt',
            'ncbi_taxonomy_tree': basencbi + PATHDELIM + 'ncbi_taxonomy_tree.txt',
            'ncbi_megan_map': basencbi + PATHDELIM + 'ncbi.map'
        }

        context.outputs = {
            'output_results_annotation_table_dir':s.output_results_annotation_table_dir,
            'output_annot_table':output_annot_table,
        }

        dbstring = self.get_dbstring()
        refdbs = [x.strip() for x in dbstring.split(",")  if len(x) != 0]

        #db_argument_string = ''
        for dbname in refdbs:
            parsed_file =  s.blast_results_dir + PATHDELIM + s.sample_name\
                            + "." + dbname+ "." + s.algorithm + "out.parsed.txt"
            context.inputs[parsed_file] = parsed_file

            #db_argument_string += ' -d ' + dbname
            #db_argument_string += ' -b ' + parsed_file


        pyScript = self.configs.CREATE_ANNOT_REPORTS
        #  --ncbi-taxonomy-map %s --ncbi-megan-map %s  --lca-gi-to-taxon-map %s"\
        #  --ncbi-taxonomy-map %s --ncbi-megan-map %s"\

        cmd = "%s  --input-annotated-gff %s  --input-kegg-maps %s \
               --input-cog-maps %s --input-seed-maps %s --input-cazy-maps %s --output-dir %s \
               --ncbi-taxonomy-map %s --ncbi-megan-map %s"\
             %(\
                pyScript, \
                context.inputs['input_annot_gff'], \
                context.inputs['KO_classification'], \
                context.inputs['COG_categories'],  \
                context.inputs['SEED_subsystems'], \
                context.inputs['CAZY_hierarchy'], \
                context.outputs['output_results_annotation_table_dir'], \
                context.inputs['ncbi_taxonomy_tree'], \
                context.inputs['ncbi_megan_map']
             )
        cmd = cmd + " -D " + s.blast_results_dir + " -s " + s.sample_name + " -a "  + s.algorithm

        #add the command now, remove to disable in a hackish way
        context.commands = [cmd]

        #context.status = self.params.get('metapaths_steps', 'CREATE_ANNOT_REPORTS')

        context.status = self.params.get('metapaths_steps','ANNOTATE_ORFS')
        context.message = self._Message("CREATING REPORT FILE FOR ORF ANNOTATION")

        contexts.append(context)
        return contexts


    def create_rpkm_cmd(self, s):
        """RPKM CALCULATION"""

        contexts = []

        '''input'''
        rpkm_input = s.rpkm_input_dir
        bwaFolder = s.bwa_folder
        output_gff = s.genbank_dir + s.sample_name + ".annot.gff"
        output_fas = s.preprocessed_dir + PATHDELIM + s.sample_name + ".fasta"
        rpkmExec = self.configs.RPKM_EXECUTABLE

        bwaExec = self.configs.BWA_EXECUTABLE


        '''output'''
        rpkm_output = s.output_results_rpkm_dir  + PATHDELIM + s.sample_name + ".orf_rpkm.txt"
        microbecensus_output = s.output_results_rpkm_dir  + PATHDELIM + s.sample_name + ".microbe_census.txt"
        stats_file = s.output_results_rpkm_dir  + PATHDELIM + s.sample_name + ".orf_read_counts_stats.txt"

        context = contextmod.Context()
        context1 = contextmod.Context()
        context.name = 'COMPUTE_RPKM'
        context.inputs = {
                           'rpkm_input':rpkm_input,
                           'output_gff': output_gff,
                           'output_fas':output_fas,
                           'bwaFolder': bwaFolder,
                         }

        context1.inputs = {
                           'rpkmExec': rpkmExec,
                           'bwaExec': bwaExec,
                         }
        context1.outputs = {
                           'rpkm_output': rpkm_output
                          }

        context.outputs = {
                           #'microbecensusoutput': microbecensus_output,
                           'stats_file': stats_file
                          }

        pyScript = self.configs.RPKM_CALCULATION

        #cmd = "%s -c %s  --rpkmExec %s --readsdir %s -O %s -o %s --sample_name  %s --stats %s --bwaFolder %s --bwaExec %s -m %s"\
        cmd = "%s -c %s  --rpkmExec %s --readsdir %s -O %s -o %s --sample_name  %s --stats %s --bwaFolder %s --bwaExec %s"\
              % (pyScript, context.inputs['output_fas'], context1.inputs['rpkmExec'],\
                 context.inputs['rpkm_input'], context.inputs['output_gff'],\
               context1.outputs['rpkm_output'],  s.sample_name, context.outputs['stats_file'],\
               context.inputs['bwaFolder'], context1.inputs['bwaExec'])

        context.status = self.params.get('metapaths_steps', 'COMPUTE_RPKM')

        context.commands = [cmd]
        contexts.append(context)
        context.message = self._Message("RUNNING RPKM_CALCULATION")
        return contexts

    def __init__(self, params, configs):

        self.params = gutils.Singleton(Params)(params)
        self.configs = gutils.Singleton(Configs)(configs)
        self.initFactoryList()

    def getContexts(self, s, stage):
        stageList  = {}
        for stageBlock in self.stageList[s.getType()]:
            for _stage in  stageBlock:
                stageList[_stage] = True

        if stage in stageList:
            return self.factory[stage](s)

    def getStageLists(self, type):
        return self.stageList[type]

    def initFactoryList(self):
        self.factory['GBK_TO_FNA_FAA_GFF'] = self.convert_gbk_to_fna_faa_gff_unannotated
        self.factory['GBK_TO_FNA_FAA_GFF_ANNOT'] = self.convert_gbk_to_fna_faa_gff_annotated
        self.factory['PREPROCESS_INPUT'] = self.create_quality_check_cmd
        self.factory['PREPROCESS_AMINOS'] = self.create_preprocess_input_aminos_cmd
        self.factory['ORF_PREDICTION'] = self.create_orf_prediction_cmd
        self.factory['ORF_TO_AMINO'] = self.create_aa_orf_sequences_cmd
        self.factory['FILTER_AMINOS'] = self.create_create_filtered_amino_acid_sequences_cmd
        self.factory['COMPUTE_REFSCORES'] = self.create_refscores_compute_cmd
        self.factory['FUNC_SEARCH'] = self.create_blastp_against_refdb_cmd
        self.factory['PARSE_FUNC_SEARCH'] = self.create_parse_blast_cmd
        self.factory['SCAN_rRNA'] = self.create_scan_rRNA_seqs_cmd
        self.factory['SCAN_tRNA'] = self.create_tRNA_scan_statistics
        self.factory['ANNOTATE_ORFS'] = self.create_annotate_genebank_cmd
        self.factory['PATHOLOGIC_INPUT'] = self.create_ptinput_cmd
        self.factory['GENBANK_FILE'] = self.create_genbank_file_cmd
        self.factory['CREATE_ANNOT_REPORTS'] = self.create_report_files_cmd
        self.factory['COMPUTE_RPKM'] = self.create_rpkm_cmd

        self.stageList['AMINO-FASTA'] = [
             [ 'PREPROCESS_AMINOS',
               'FILTER_AMINOS'
             ],
             ['FUNC_SEARCH'],
             ['COMPUTE_REFSCORES',
              'PARSE_FUNC_SEARCH',
              'ANNOTATE_ORFS',
              'CREATE_ANNOT_REPORTS',
              'PATHOLOGIC_INPUT'
             ]
             ]

        self.stageList['NUCL-FASTA'] = [
             ['PREPROCESS_INPUT',
              'ORF_PREDICTION',
              'ORF_TO_AMINO',
              'FILTER_AMINOS'],
             ['FUNC_SEARCH'],
             ['COMPUTE_REFSCORES' ,
              'PARSE_FUNC_SEARCH',
              'SCAN_rRNA',
              'SCAN_tRNA',
              'ANNOTATE_ORFS',
              'CREATE_ANNOT_REPORTS',
              "GENBANK_FILE",
              'PATHOLOGIC_INPUT',
              'COMPUTE_RPKM']
            ]

        self.stageList['AMINO-GENBANK-UNANNOT'] = [
            ['GBK_TO_FNA_FAA_GFF',
             'FILTER_AMINOS'],
            ['FUNC_SEARCH'],
            ['COMPUTE_REFSCORES' ,
             'PARSE_FUNC_SEARCH',
             'SCAN_rRNA',
             'SCAN_tRNA',
             'ANNOTATE_ORFS',
             'CREATE_ANNOT_REPORTS',
             'PATHOLOGIC_INPUT'
            ]
        ]

        self.stageList['AMINO-GENBANK-ANNOT'] = [
            'GBK_TO_FNA_FAA_GFF_ANNOT',
            'FILTER_AMINOS',
            'SCAN_rRNA',
            'SCAN_tRNA',
            'PATHOLOGIC_INPUT',
        ]

class JobCreator():

    """ This class looks up the steps status redo, yes, skip and decided which steps
        needs to be added into the job/context list """
    params = {}
    configs = {}
    stagesInOrder = []

    def  __init__(self, params, configs):
        self.params = params
        self.configs = configs
        #print self.configs

    def addJobs(self, s, block_mode = False):
        contextCreator = ContextCreator(self.params, self.configs)

        contextBlock = []
        for stageList in contextCreator.getStageLists(s.getType()):
            if block_mode == True:
                contextBlock = []

            for stage in stageList:
                if stage in self.params['metapaths_steps'] or\
                    stage in ['ORF_TO_AMINO', 'GBK_TO_FNA_FAA_GFF', 'GBK_TO_FNA_FAA_GFF_ANNOT',\
                           'COMPUTE_REFSCORES', 'PREPROCESS_AMINOS', 'PATHOLOGIC_INPUT',\
                           'CREATE_ANNOT_REPORTS']:
                    # if self.params['INPUT']['format'] =='gbk-unannotated':
                    #   if stage =='PREPROCESS_INPUT':
                    #     stage = 'GBK_TO_FNA_FAA_GFF'
                    #   self.stageList['gbk-unannotated'] = ['GBK_TO_FNA_FAA_GFF',
                    contexts = contextCreator.getContexts(s, stage)
                    contextBlock.append(contexts)
            if block_mode == True:
                s.addContexts(contextBlock)

        # block stages
        if block_mode == False:
            s.addContexts(contextBlock)



