Installation
************

Metapathways can be installed with conda or obtained as a container. In either case, reference
databases must also be downloaded as indicated below. 

`build-db` parameters:
   - `--func` databases to set up, at least one must be selected
   - `-a` aligner for which to set up the databases for, can be either `fast` or `blast`
   - `-d` path ot directory to save reference databases to
   - `-t` threads to Use

.. note::
   
   `uniref90` and `uniref50` are the largest databases at ~270 GB and ~30 GB respectively after set up.
   Others are less than ~10GB each.

Conda
-----

.. code-block:: bash

   conda install -c hallamlab -c bioconda -c conda-forge metapathways
   metapathways build-db \
      -t ${threads} \
      -d ${path_to_save_reference_databases_to} \
      --func metacyc swissprot cazy eggnog uniref50 uniref90 \
      -a fast # or blast

Apptainer
---------

.. code-block:: bash
   
   apptainer build metapathways.sif docker://quay.io/hallamlab/metapathways:latest
   apptainer run --bind ${path_to_save_reference_databases_to}:/ref \
      metapathways build-db \
         -t ${threads} \
         -d /ref \
         --func metacyc swissprot cazy eggnog uniref50 uniref90 \
         -a fast # or blast
