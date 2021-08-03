Installation
************

Setup a virtual environment
===========================

**Create a python virtual environment** 
python Virtual enviroments `venv` (for Python 3) allow you to manage separate 
package installations for different projects. They essentially allow you to create 
a “virtual” isolated Python installation and install packages into that virtual 
installation. When you switch projects, you can simply create a new virtual 
environment and not have to worry about breaking the packages installed in 
the other environments. It is always recommended to use a virtual environment 
while trying out new Python applications.

The following command creates a new virtual environment with a name *mynewenv* with Python 3
::

 $ virtualenv -p /usr/bin/python3 mynewenv

Activate the new virtual environment by running 
::

 $ source mynewenv/bin/activate

Deactivate If you want to switch projects or otherwise leave your virtual environment, simply run:
::

  $ deactivate

pip install MetaPathways
========================
Install MetaPathways by running:
::

 $ pip3 install metapathways

To make sure MetaPathways is installed type
::

 $ MetaPathways --version

which, if MetaPathways, is properly installed, will print a version number. For example
::

  MetaPathways: Version 3.5.0


Install Binaries
================

Next we install ``trnascan-1.4``, ``rpkm``, ``prodigal``, ``FAST`` and ``bwa``

..
   Download the source code as
   ::

    $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/c_cpp_sources.1.0.tar.gz

   untar the files, make and install, which takes a few minutes 
   ::

     $ tar -zxvf c_cpp_sources.1.0.tar.gz
     $ cd c_cpp_sources
     $ make`
     $ sudo make install

::
   make extensions-build
   make extensions-install

NOTE: If you lack the privileges to install these binaries system-wide, then you will need to pass the ``DESTDIR`` environment variable when calling ``make``:
::

   $ make DESTDIR=/home/username extensions-install

In this example, ``make`` will install the binaries into ``/home/username/bin``, so you must have the requisite permissions to modify the ``/home/username`` directory. It will create the ``bin`` directory if it does not already exist, and copy the binaries there.


NOTE: if you would like to unstall then type
::
   
  $ sudo make uninstall


Install ``ncbi-blast+`` locally. Visit the `download page
<https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`_.

For Ubuntu/Debian
::

  $ sudo apt-get install ncbi-blast+


Reference Sequences
===================

We use ``Snakemake`` to automate the staging of reference databases
needed by MetaPathways. We have installed ``Snakemake`` via Conda. If
you are using the Docker container, then Conda is already
initialized.

If you are using the container via Singularity, you must
first initialize Conda as follows (note the space between the period
character, and the first slash character):
::
  
   . /opt/conda/etc/profile.d/conda.sh

Now, for both Docker and Singularity, we can activate the ``Snakemake`` environment:
::
   conda activate snakemake

Next, we need to switch to the MetaPathways install directory,  with the ``Snakefile`` file:
::
   cd /opt/mp_repo

   
Now, we have ``Snakemake`` automate the installation of the required files:
::
   snakemake --cores 1 --config ref_db_dir=/path/to/db/dir -- stage_blast_full

Using the `--config` command, we can specify the desired root
directory for installing the MetaPathways reference databases
(replacing `/path/to/db/dir` with a real directory path on your
system). By default, the ``Snakemake`` configuration file sets
`/tmp/mp_ref_dbs` as the root directory for the MetaPathways reference
database, if you leave off the `--config` option.

Instead of using a single core, you can use the `--cores` option with
a greater number of specified cores to parallelize the staging of the
requested datasets.




      



Create the following reference folder structure under a folder. Here we use the 
example name ``MetaPathways_DBs``
::

 $ mkdir -p MetaPathways_DBs/taxonomic/formatted
 $ mkdir -p MetaPathways_DBs/functional/formatted 
 $ mkdir -p MetaPathways_DBs/ncbi_tree 
 $ mkdir -p MetaPathways_DBs/functional_categories

::

   MetaPathways_DBs/
   ├── functional
   │   ├── formatted
   ├── functional_categories
   ├── ncbi_tree
   └── taxonomic
       └── formatted

Download and unzip the NCBI taxonomy file to the ``MetaPathways_DBs/ncbi_tree`` folder
::

 $ cd MetaPathways_DBs/ncbi_tree
 $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/ncbi_taxonomy_tree.txt.gz
 $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/ncbi.map.gz

Download and unzip functional classification files to ``MetaPathways_DBs/functional_hierarchy`` folder
::

$ cd MetaPathways_DBs/functional_hierarchy
$ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/CAZY_hierarchy.txt.gz
$ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/COG_categories.txt.gz
$ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/KO_classification.txt.gz
$ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/SEED_subsystems.txt.gz

and we should see the following structure 
::

   MetaPathways_DBs/
   ├── functional
   │   ├── formatted
   ├── functional_categories
   │   ├── CAZY_hierarchy.txt.gz
   │   ├── COG_categories.txt.gz
   │   ├── KO_classification.txt.gz
   │   ├── SEED_subsystems.txt.gz
   ├── ncbi_tree
   │   ├── ncbi_taxonomy_tree.txt.gz
   │   ├── ncbi.map.gz
   └── taxonomic
       └── formatted

Functional Reference 
++++++++++++++++++++

The functional references are protein reference sequences used for functional and taxonomic
annotation. Any set of protein references in the FASTA format can be used, e.g., we show 
a few lines
::

  >WP_096046812.1 hypothetical protein [Sulfurospirillum sp. JPD-1]
  MSKKAFLFLILLVMSLQSLLVACGGSCLECHSKLRPYINDQNHAILNECITCHNQPSKNGQCGRDCFDCHSQEKVYAQKDVNAHQELKT
  CGTCHKEKVDFTTPKQSIISNQQNLIHLFK
  >WP_096046815.1 hypothetical protein [Sulfurospirillum sp. JPD-1]
  MKKLLIILALISRLIAEDSSDLDEIKEEDIPKILSIIKDGTKEHLPMMLDDYTTLVDIVSVNNAIEYRNRINSANEHVKTILKADKGTLI
  KTTFDNNKSYLCSDYETRSLLKKGAVFIYVFYDMNNAELFKFSIQEKDCQ
  >WP_016244176.1 hypothetical protein [Escherichia coli]
  MTDITDRHTLRRMSWSELFTAAQEAEFQRDYERARIVWSFALHVATTTINKNLSIAHIRRCDTLLHKSKTVPGNNTGGRSVCLRPQHPRR 
  ...........


Formatting Reference Sequences
++++++++++++++++++++++++++++++

 For the purpose of demonstration we walk you through the process of preparing a
 small set of protein reference sequences from the NCBI Refseq protein databases.
 Download the example protein reference sequence file `refseq-mini.fasta.gz`
 to the functional folder as follows

::

 $ cd MetaPathways_DBs/functional
 $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/refdata/refseq-mini.fasta.gz
 $ gunzip refseq-mini.fasta.gz

rename to remove the `fasta` suffix
::

 $ mv refseq-mini.fasta  refseq-mini
 $ cat refseq-mini | grep ">" > formatted/refseq-mini-names.txt

FAST
----

BLAST
-----
Format the database for `blastp` as follows:
::

  $ cd MetaPathways_DBs/functional
  $ makeblastdb -dbtype prot -in refseq-mini -out formatted/refseq-mini

Taxonomic Reference 
+++++++++++++++++++
