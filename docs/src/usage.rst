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

Supported Functional Databases
==============================

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

Minimal Example:
~~~~~~~~~~~~~~~~

.. code-block:: bash

   metapathways build-db \
      -t ${threads} \
      -d ${path/to/save/reference_databases} \
      --func swissprot \
      -a fast

Annotate a Metagenome
---------------------

MetaPathways has many parameters and flags to allow for explicit control of many aspects of proccessing.
However, the minimal (default) analysis requires very few user inputs.

Minimal Input
=============

MetaPathways inputs are fasta files provided in an input folder. The file names must end with 
a `.fasta` or `.fas`. These fasta files contains the contigs or DNA sequences from assembling.

..
   Parameter File 
   ==============

   The parameter file must indicate the setting for any MetaPathways run. An example paramter file 
   can be downloaded as
   ::

    $ wget  https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/text/template_param.txt

   Below we describe the settings in the parameter file. 

   Run
   ===

    As an illustration we donwload a small input file `testsample1.fasta`
    in a folder named `mp_input` and we want the output in a folder names `mp_output`

   ::

    $ mkdir mp_input
    $ cd mp_input
    $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/testdata/testsample1.fasta
    $ cd ..

   Now we kick off a run  as
   ::

     $ MetaPathways --input mp_input --output mp_output -p template_param.txt -d ~/MetaPathways_DBs/






