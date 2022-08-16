Installation
************

MetaPathways supports installing the software using Conda and Pip in a 
64-bit Linux environment, or from a
container image that can be used with Docker or Singularity. If you do
not have administrator (i.e., "root") access to your computer, we
recommend that users install MiniConda if they do not already have it
set up. For users wanting to use MetaPathways in an academic grid
computing environment, we recommend using the container image *via*
Singularity. Below please find a description of how to install MetaPathways 
using the two supported options:


Container Install
=================

Our container images are hosted at `Quay.io <https://quay.io/repository/hallamlab/metapathways?tab=info>`_. 
The following commands assume that you are already familiar with installing and running Docker containers 
via the ``docker`` or ``singularity`` executables:

Using `Docker <https://sylabs.io/>`_::
     sudo docker pull quay.io/hallamlab/metapathways

Using `Singularity <https://sylabs.io/>`_::
   singularity build metapathways.sif docker://quay.io/hallamlab/metapathways:latest

More advanced container-related commands are available as Make targets in the ``Makefile``.

Installing with Pip and Conda
=============

We currently offer a way to use Pip to install the MetaPathways Python package,
along with using `Conda <https://conda.io`_ to install all dependencies. We do not yet have a Conda 
package for MetaPathways. It is in the works for a future release.

For this to work, we assume that you have the following already set up in your
command line environment:

* You have basic ``C`` development tools in your path, including: ``make``, ``ar``, ``gcc``, ``g++``, etc.
* You have a development verson of the ``zlib`` library
* You have a working version of ``git``
* You have Python 3 (``python3``) and ``pip3`` installed
* You have already `installed Conda <https://docs.conda.io/en/latest/miniconda.html>`_, and it is activated
* You have ``wget`` installed

If you are using a version of Linux that uses ``apt``, and you have root access, then you can execute the 
following to get all of the dependencies except Conda:
::
   sudo apt-get update -y
   sudo apt-get install -y make \
                  binutils \
                  python3 \
                  python3-pip \
                  zlib1g-dev \
                  wget

First, clone the MetaPathways repo into the current directory:
::
   git clone https://taltman1@bitbucket.org/BCB2/metapathways.git
   
Then, change into the top-level directory of the MetaPathways repository:
::
   cd metapathways
   
If you have root or administrator access, execute the following Make target in the ``Makefile`` in the top
level of the MetaPathways repository to use Conda to install some of the dependencies. Aside from the 
Conda-based dependencies, wo binary executables will be installed in the ``/usr/local/bin`` directory:
::
  make conda-install
  
If you do not have root or administrator access, use this version, where ``/home/user/example`` 
is a directory where you want the binaries to be installed. In that directory, a ``bin``
directory will be created, if it does not already exist. It is assumed that, following
this example, that ``/home/user/example/bin`` is in your ``$PATH`` environment variable:
::
   make DESTDIR=/home/user/example conda-install

If you have root/administrator access, install the MetaPathways Python package using the following command:
::
   pip3 install git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways

Else, use this form to install the package to the user's home directory:
::
   pip3 install --user git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways
   
If installing with the ``--user`` option, make sure to add ``$HOME/.local/bin`` to your ``$PATH`` environment variable. This will allow you 
to use the programs without having to type the full path each time.


Reference Sequences
===================

MetaPathways relies on reference databases of sequences to assign 
functional and taxonomic annotations to the user's sequences. The 
reference databases, and the index files for each database, take up
a significant amount of disk storage. See below for an anecdotal 
example.


You cannot install these large reference databases within the
container, though. You should have a directory on a disk with plenty
of capacity, and use Docker's and Singularity's bind options to mount
that external directory within the container. Here's an example using
Singularity:
::
   singularity shell --bind /mnt/sandbox/user:/data docker://quay.io/hallamlab/metapathways:latest

The above example binds the host operating system's
``/mnt/sandbox/user`` directory within the running container as
``/data``.

*Warning*: Circa 2021-10, using a beefy computer with many cores and
plenty of RAM, performing the staging of the full Blast databases may
take an hour, and staging the full set of FAST databases will
take around *24 hours*. The Blast ``refseq_protein`` databases take up ~90 GB of disk
capacity, while the FAST ``refseq_protein`` database takes up ~375
GB. The combination of other staged databases (including both Blast
and FAST versions) consumes an additional ~20 GB. Please make sure you
have adequate disk capacity before starting the database staging.

We use ``Snakemake`` to automate the staging of reference databases
needed by MetaPathways. We have installed ``Snakemake`` via Conda. If
you are using the Docker container, then Conda is already
initialized. If you are using the container via Singularity, you must
first initialize Conda as follows (note the space between the period
character, and the first slash character):
::
  
   . /opt/conda/etc/profile.d/conda.sh
   
Now we can activate the ``Snakemake`` environment:
::
   conda activate snakemake

Next, we need to switch to the MetaPathways install directory,  with the ``Snakefile`` file:
::
   cd /opt/mp_repo

Now, we have ``Snakemake`` automate the installation of the required files:
::
   snakemake --cores 1 --config ref_db_dir=/path/to/db/dir -- stage_blast_lite

Using the ```--config`` command, we can specify the desired root
directory for installing the MetaPathways reference databases
(replacing ``/path/to/db/dir`` with a real directory path on your
system). By default, the ``Snakemake`` configuration file sets
``/tmp/mp_ref_dbs`` as the root directory for the MetaPathways reference
database, if you leave off the ``--config`` option. 

Instead of using a single core, you can use the ``--cores`` option with
a greater number of specified cores to parallelize the staging of the
requested datasets.

Above we issued the ``stage_fast_lite`` command to Snakemake, as an example
that runs quickly. There are actually four options for staging the data:

* All databases, indexed for use with Blast: ``stage_blast_full``
* All databases except RefSeq Proteome, indexed for use with Blast:
  ``stage_blast_lite``
* All databases, indexed for use with FAST: ``stage_fast_full``
* All databases except RefSeq Proteome, indexed for use with FAST:
  ``stage_fast_lite``

So, first decide whether you want to use Blast or FAST, and then
decide whether you have the disk space and the install time to install
the NCBI RefSeq Proteome reference database. FAST runs faster than Blast,
with comparable sensitivity. And the RefSeq Proteome is currently required
for MetaPathways to accurately annotate contigs taxonomically. Thus, we 
recommend running ``stage_fast_full``, if you have the disk storage and 
the time to let it run.
      
