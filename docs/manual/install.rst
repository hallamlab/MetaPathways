Installation
************

MetaPathways supports installing the software *via* Conda or from a
container image that can be used with Docker or Singularity. If you do
not have administrator (i.e., "root") access to your computer, we
recommend that users install MiniConda if they do not already have it
set up. For users wanting to use MetaPathways in an academic grid
computing environment, we recommend using the container image *via*
Singularity.


Container Install
=================


Using Docker::
     sudo docker pull quay.io/hallamlab/metapathways

Using Singularity::
   singularity build metapathways.sif docker://quay.io/hallamlab/metapathways:dev

More container-related commands are available as Make targets in the ``Makefile``.

Conda Install
=============

Execute the following Make targets in the ``Makefile`` in the top
level of the MetaPathways repository. Assuming that you do not have
root access on your computer, you can use the ``DESTDIR`` environment
variable to specify where to put the executables::
  make DESTDIR=/home/user/bin conda-install

In this example, we assume that the user ``user`` has a ``bin``
directory in their home directory, and that this path in in their
``$PATH`` environment variable.

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
      
