Full Usage
**********

Build Reference Databases
-------------------------

Metapathways requires reference databases to perform functional/taxonomic annotation.
Below provides the commands for building currently supported database.

.. note::
   
   `uniref90` and `uniref50` are the largest databases at ~270 GB and ~30 GB respectively after set up.
   Others are less than ~10GB each.
   
   `SILVA` is the only supported taxonomic reference database at this time 
   and therefore is built automatically along with functional references.

Minimal Usage:
~~~~~~~~~~~~~~

.. code-block:: bash

   metapathways build-db \
      -t ${threads} \
      -d ${path/to/save/reference_databases} \
      --func swissprot \
      -a fast

`build-db` parameters:
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   usage: metapathways [-h] [-d PATH] [--func [CATEGORICAL ...]] [-a ALIGNER]
                       [-t INT] [--dryrun] [--snakemake [SNAKEMAKE ...]] [--test]

   automated database install

   options:
     -h, --help            show this help message and exit
     -t INT, --threads INT
                           max number of cores to use in multithreaded steps [1]
     --dryrun              dry run snakemake
     --snakemake [SNAKEMAKE ...]
                           additional snakemake cli args in the form of KEY="VALUE"
                           or KEY (no leading dashes)
     --test                use test values for all arguments

   database arguments:
     -d PATH, --refdb_dir PATH
                           path to save the reference DB, [DEFAULT "./"]
     --func [CATEGORICAL ...]
                           functional references, space-delimited list from 
                           ['metacyc', 'swissprot', 'cazy', 'eggnog', 'uniref50', 'uniref90'],
                           [DEFAULT ['metacyc', 'swissprot']]
     -a ALIGNER, --aligner ALIGNER
                           local aligner to index for, select one of ['fast', 'blast'],
                           [DEFAULT fast]

Supported Functional Databases:
===============================

+-----------+---------------------------------------------------+---------------------+
| Database  | Description                                       | Size (after setup)  |
+===========+===================================================+=====================+
| uniref90  | UniRef90 functional annotation database           | ~270 GB             |
+-----------+---------------------------------------------------+---------------------+
| uniref50  | UniRef50 functional annotation database           | ~30 GB              |
+-----------+---------------------------------------------------+---------------------+
| swissprot | SwissProt functional annotation database          | <10 GB              |
+-----------+---------------------------------------------------+---------------------+
| metacyc   | MetaCyc functional annotation database            | <10 GB              |
+-----------+---------------------------------------------------+---------------------+
| cazy      | CAZymes functional annotation database            | <10 GB              |
+-----------+---------------------------------------------------+---------------------+

.. note::

   **METACYC**: When building MetaCyc you require a username and password due the their license.
   During the build process there will be a prompt to enter your username, followed by
   a prompt to enter your password. This will allow MetaPathways to access the download
   for MetaCyc and build the database.

Annotate a Metagenome
---------------------

MetaPathways has many parameters and flags to allow for explicit control of many aspects of proccessing.
However, the minimal (default) analysis requires very few user inputs.

Minimal Usage:
~~~~~~~~~~~~~~

.. code-block:: bash

   metapathways run \
      -i $[input_metagenome.fa] \
      -d ${path/to/save/reference_databases} \
      -o ${path/to/output} \
      -t ${threads}

`run` parameters:
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   usage: Metapathways run [options]

   Minimum REQUIRED Command:
   MetaPathways run -i INPUT_FILE -o OUTPUT_DIR -d REFDB_DIR
   
   options:
     -h, --help            show this help message and exit
   
   Minimum Required Arguments:
     -i INPUT_FILE, --input_file INPUT_FILE
                           path to the input fasta file/input dir [REQUIRED]
     -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                           path to the output directory [REQUIRED]
     -d REFDB_DIR, --refdb_dir REFDB_DIR
                           path to the reference DB [REQUIRED]
   
   Quality Controls Arguments:
     --input_format {fasta,fasta-amino}
                           Input format, FASTA support only [fasta]
     --qc_min_length QC_MIN_LENGTH
                           Minimum length for quality control [180]
     --qc_delete_replicates {yes,no}
                           Delete replicates in quality control [yes]
   
   ORF Prediction Arguments:
     --orf_strand {pos,neg,both}
                           Strand for ORF prediction [both]
     --orf_algorithm ORF_ALGORITHM
                           Algorithm for ORF prediction, Prodigal support only [prodigal]
     --orf_min_length ORF_MIN_LENGTH
                           Minimum ORF length [60]
     --orf_translation_table ORF_TRANSLATION_TABLE
                           Translation table for ORF prediction, see Prodigal for translation tables [11] 
     --orf_mode {single,meta}
                           Mode for ORF prediction [meta]
   
   Functional Annotation Arguments:
     --annotation_algorithm {FAST,BLAST}
                           Algorithm for ORF annotation [FAST]
     --annotation_dbs ANNOTATION_DBS [ANNOTATION_DBS ...]
                           Database(s) for annotation, space-separated list [swissprot]
     --annotation_min_bsr ANNOTATION_MIN_BSR
                           Minimum BSR for annotation [0.4]
     --annotation_max_evalue ANNOTATION_MAX_EVALUE
                           Maximum e-value for annotation [0.000001]
     --annotation_min_score ANNOTATION_MIN_SCORE
                           Minimum score for annotation [20]
     --annotation_min_length ANNOTATION_MIN_LENGTH
                           Minimum length for annotation [45]
     --annotation_max_hits ANNOTATION_MAX_HITS
                           Maximum hits for annotation [5]
     --annotation_run_mode {default,pervol}
                           Run mode for annotation, FAST only [pervol]
   
   rRNA Annotation Arguments:
     --rRNA_refdbs RRNA_REFDBS [RRNA_REFDBS ...]
                           Reference databases for rRNA annotation, space-separated list
                           [SILVA_138.1_LSURef_NR99_tax_silva_trunc SILVA_138.1_SSURef_NR99_tax_silva_trunc]
     --rRNA_max_evalue RRNA_MAX_EVALUE
                           Maximum e-value for rRNA annotation [0.000001]
     --rRNA_min_identity RRNA_MIN_IDENTITY
                           Minimum identity for rRNA annotation [20]
     --rRNA_min_bitscore RRNA_MIN_BITSCORE
                           Minimum bitscore for rRNA annotation [50]
   
   Read Mapping Arguments (single sample support only):
     -1 FWD_FASTQ, --fastq FWD_FASTQ
                           location of the raw fastq file, either forward or interleaved
     -2 REV_FASTQ, --rev_fastq REV_FASTQ
                           location of the raw reverse fastq file, if separate paired-end
     --interleaved         if paired-end is interleaved [False]

   Pipeline Step Arguments:
     --PREPROCESS_INPUT {yes,skip,redo}
                           Step: PREPROCESS_INPUT [yes]
     --ORF_PREDICTION {yes,skip,redo}
                           Step: ORF_PREDICTION [yes]
     --FILTER_AMINOS {yes,skip,redo}
                           Step: FILTER_AMINOS [yes]
     --SCAN_rRNA {yes,skip,redo}
                           Step: SCAN_rRNA [yes]
     --SCAN_tRNA {yes,skip,redo}
                           Step: SCAN_tRNA [yes]
     --FUNC_SEARCH {yes,skip,redo}
                           Step: FUNC_SEARCH [yes]
     --PARSE_FUNC_SEARCH {yes,skip,redo}
                           Step: PARSE_FUNC_SEARCH [yes]
     --ANNOTATE_ORFS {yes,skip,redo}
                           Step: ANNOTATE_ORFS [yes]
     --GENBANK_FILE {yes,skip,redo}
                           Step: GENBANK_FILE [yes]
     --CREATE_ANNOT_REPORTS {yes,skip,redo}
                           Step: CREATE_ANNOT_REPORTS [yes]
     --PATHOLOGIC_INPUT {yes,skip,redo}
                           Step: PATHOLOGIC_INPUT [yes]
     --COMPUTE_TPM {yes,skip,redo}
                           Step: COMPUTE_TPM [yes]
     --force_redo          Redo all steps [False]

   Miscellaneous Arguments:
     -s SAMPLES [SAMPLES ...], --samples SAMPLES [SAMPLES ...]
                           process only specific samples, space-separated list
     -t THREADS, --threads THREADS
                           max number of cores to use in multithreaded steps [1]
     -v, --verbose         print more information on the stdout
     --test                use test values for all arguments

Add MAGs to MP output
---------------------

MetaPathways allows you to add MAG contig mappings to the annotated metagenome output.
Simple using the mag splitting command with a mapping file that maps the original contig
IDs (column 1) to their respective MAGs (column 2).

.. note::

   Be sure to annotated the metagenome first, then use the output dirctory as the (-o) of 
   the `mag_split` command.

Minimal Usage:
~~~~~~~~~~~~~~

.. code-block:: bash

   metapathways mag_split \
      -o ${path/to/output}/{metagenome_id} \
      -m ${contig2mag_mapping}
      
`mag_split` parameters:
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   usage: Metapathways mag_split [options]

   Minimum REQUIRED Command:
   Metapathways mag_split -o output_dir -m contig_mag_map

   options:
     -h, --help            show this help message and exit
     -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                           path where MP output was saved [REQUIRED]
     -m MAG_MAP, --contig_mag_map MAG_MAP
                           TSV file that contains contig-to-MAG mapping [REQUIRED]

Contig mapping format:
======================

+---------+------+
| contig1 | MAG1 |
+---------+------+
| config2 | MAG2 |
+---------+------+
| contig3 | MAG2 |
+---------+------+
| contig4 | MAG3 |
+---------+------+
| contig5 | MAG1 |
+---------+------+

Environmental Pathway Genome Databases (ePGDBs)
-----------------------------------------------

MetaPathways supports the creation of ePGDBs by interfacing with the 
`Pathway Tools software <http://brg.ai.sri.com/ptools/>`_.
In order to use this feature, you must first install Pathway Tools which requires a license.
Pathway Tools is freely available for research purposes to academic, non-profit, and government institutions, 
and is available for a fee to commercial institutions. For academic use, a license and the download can be
found here: `Download <https://biocyc.org/download.shtml>`_.

.. note::

   Although MetaPathways is installed in a conda/mamba environment, Pathway Tools is installed in the 
   $HOME directory by default and is globally available in the commandline. MetaPathways expects that
   the Pathway Tools software has been installed using defaults and may not work as expected if installed
   in a custom way.

Minimal Usage:
~~~~~~~~~~~~~~

.. code-block:: bash

   metapathways ptools \
      -o ${path/to/output}/{metagenome_id}	
      
`ptools` parameters:
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   usage: Metapathways ptools [options]

   Minimum REQUIRED Command:
   Metapathways ptools -o output_dir

   options:
     -h, --help            show this help message and exit
     -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                           path where MP output was saved [REQUIRED]
     --tag TAG             Custom name for ePGDB [optional]
     --container           Flag only used in containerized env [special flag]
     --taxprune            Set taxonomic pruning in pathway tools to True

.. note::
   The container option only functions if both MetaPathways and Pathway Tools are installed
   in the same container. But due to licensing constraints this cannot be distributed and therefore
   is only experimental at this point.
   
   *There will be more work on this in the near future*
   
   
