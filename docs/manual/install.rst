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

You'll want to install these large reference databases not within the
container, though. You should have a directory on a disk with plenty
of capacity, and use Docker's and Singularity's bind options to mount
that external directory within the container. Here's an example using
Singularity:

::
   singularity shell --bind /mnt/sandbox/user:/data docker://quay.io/hallamlab/metapathways:dev

The above example binds the host operating system's
`/mnt/sandbox/user` directory within the running container as
`/data`.

Warning: Circa 2021-10, using a beefy computer with many cores and
plenty of RAM, performing the staging of the full Blast databases may
take an hour, and staging the full set of FAST databases will
take around *24 hours*. The Blast `refseq_protein` databases take up ~90 GB of disk
capacity, while the FAST `refseq_protein` database takes up ~375
GB. The combination of other staged databases (including both Blast
and FAST versions) consumes an additional ~20 GB. Please make sure you
have adequate disk capacity before starting the database staging.

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

Above we issued the `stage_blast_full` command to Snakemake. There are
actually four options for staging the data:

* All databases, indexed for use with Blast: `stage_blast_full`
* All databases except RefSeq Proteome, indexed for use with Blast:
  `stage_blast_lite`
* All databases, indexed for use with FAST: `stage_fast_full`
* All databases except RefSeq Proteome, indexed for use with FAST:
  `stage_fast_lite`

So, first decide whether you want to use Blast or FAST, and then
decide whether you have the disk space and the install time to install
the NCBI RefSeq Proteome reference database. If you have plenty of
both, you can actually install both full sets for Blast and FAST by
first issuing the `stage_blast_full` and then `stage_fast_full` to
Snakemake, one command at a time.
      




After everything is staged, we should see the following structure:
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


