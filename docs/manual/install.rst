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

Using `Docker <https://sylabs.io/>`_:
::
     sudo docker pull quay.io/hallamlab/metapathways

Using `Singularity <https://sylabs.io/>`_:
::
   singularity build metapathways.sif docker://quay.io/hallamlab/metapathways:latest

More advanced container-related commands are available as Make targets in the ``Makefile``.

Installing with Pip and Conda
=============================

Summary
-------

Assuming that you have all prerequisites satisfied, installng can be as simple as:
::
   conda create --name metapathways
   conda activate metapathways
   pip3 install git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways
   metapathways-install-deps.sh
   
Read on to learn the details.

Detailed Install
----------------

We currently offer a way to use Pip to install the MetaPathways Python package,
along with using `Conda <https://conda.io>`_ to install all dependencies. We do not yet have a Conda 
package for MetaPathways. It is in the works for a future release.

For this to work, we assume that you have the following already set up in your
command line environment:

* You have Python 3 (``python3``) and ``pip3`` installed
* You have already `installed Conda <https://docs.conda.io/en/latest/miniconda.html>`_, and it is activated
* Development files for ``zlib``, ``liblzma` and ``libbz2`` (required to install PySAM via pip)
* You have ``wget`` installed

If you are using a version of Linux that uses ``apt``, and you have root access, then you can execute the 
following to get all of the dependencies except Conda:
::
   sudo apt-get update -y
   sudo apt-get install -y \
                  python3 \
                  python3-pip \
		  zlib1g-dev \
		  liblzma-dev \
		  libbz2-dev \
                  wget

Installing Python Package as Root
................................

If you have root/administrator access, install the MetaPathways Python package using the following command:
::
   pip3 install git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways

Installing Python Package as an Unpriviledged User
..................................................

Use this form to install the package to the user's home directory:
::
   pip3 install --user git+https://bitbucket.org/BCB2/metapathways.git@dev#egg=MetaPathways
   
Make sure to add ``$HOME/.local/bin`` to your ``$PATH`` environment variable. This will allow you 
to use the programs without having to type the full path each time.

Conda-Based Setup
.................

Once you have installed the Python package, you will have the following executables either in the system Python install path, or in ``~/.local/bin``, so be sure to add those paths to your ``$PATH`` environment variable.
::
   MetaPathways
   metapathways-install-deps.sh
   metapathways-data-install.sh
   metacount
   fastal
   fastdb

Execute ``metapathways-install-deps.sh`` to install pipeline dependencies using Conda.


Reference Sequences
===================

Summary
-------
Assuming that you have MetaPathways installed, installing the reference DB can be as simple as:
::
   metapathways-data-install.sh /media/ref-db-dir stage_fast_full
   
Read on for detailed instructons.

Details
-------

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

Now we can run the ``metapathways-data-install.sh`` script
::
   metapathways-data-install.sh /media/ref-db-dir stage_fast_lite

... where ``/media/ref-db-dir`` is the reference database installation directory (make sure this directory has adequate capacity for the data to be installed).

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
      
