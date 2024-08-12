Quick Start (TLDR;)
*******************
Are you a reviewer or just need to test this out quickly?
The commands below will allow you to install and test MP-3.5 quickly and efficiently.
Sorry you still need Conda/Mamba.

.. note::

   If you run the last two commands, make sure to replace the ${variables} with real values.

Try it with Mamba (recommended)
===============================

.. code-block:: bash

   # install metapathways with Mamba
   mamba create -n metapathways_env -c hallamlab -c bioconda -c conda-forge metapathways
   mamba activate metapathways_env
   
   # test the install
   metapathways build_db --test
   metapathways run --test # results in `${CWD}/test`
      
   ### Try your own data! ###
   
   # build the minimal reference DB
   metapathways build_db -d ${path/to/MPDB} --func swissprot -a fast
   
   # run minimal DB on your personal data
   metapathways run -i ${your_metagenome.fa} -o ${path/to/output_dir} -d ${path/to/MPDB}

.. note::

   The build_db test does not create a DB that is to be used for downstream DB builds.
   For the full reference build the use is to provide a location with acceptable storage
   capacity and accessibility.
   
Try it with Docker
=====================

.. code-block:: bash

   # pull down with docker
   docker pull quay.io/hallamlab/metapathways
   docker run -it quay.io/hallamlab/metapathways /bin/bash
   
   # test the install (within docker)
   metapathways build_db --test
   metapathways run --test # results in `${CWD}/test`
   
   ### Try your own data! ###
   
   # build the minimal reference DB
   metapathways build_db -d ${path/to/MPDB} --func swissprot -a fast
   
   # run minimal DB on your personal data
   metapathways run -i ${your_metagenome.fa} -o ${path/to/output_dir} -d ${path/to/MPDB}

