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




      

