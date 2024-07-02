Build Reference Databases
************

Metapathways requires reference databases to perform functional/taxonomic annotation.
Below provides the commands for building currently supported database.
.. note::
   
   `uniref90` and `uniref50` are the largest databases at ~270 GB and ~30 GB respectively after set up.
   Others are less than ~10GB each.
   
   `SILVA` is the only supported taxonomic reference database at this time 
   and therefore is built automatically along with functional references.

`build-db` parameters:
   - `--func` databases to set up, at least one must be selected
   - `-a` aligner for which to set up the databases for, can be either `fast` or `blast`
   - `-d` path ot directory to save reference databases to
   - `-t` threads to Use


Conda
-----

.. code-block:: bash

   metapathways build-db \
      -t ${threads} \
      -d ${path_to_save_reference_databases_to} \
      --func swissprot # minimal example\
      -a fast # or blast


Apptainer
---------

.. code-block:: bash
   
   apptainer build metapathways.sif docker://quay.io/hallamlab/metapathways:latest
   apptainer run --bind ${path_to_save_reference_databases_to}:/ref \
      metapathways build-db \
         -t ${threads} \
         -d /ref \
         --func swissprot # minimal example \
         -a fast # or blast


Supported Functional Databases
------------------------------
| Database   | Description                                       | Size (after setup)  |
|------------|---------------------------------------------------|---------------------|
| uniref90   | UniRef90 functional annotation database           | ~270 GB             |
| uniref50   | UniRef50 functional annotation database           | ~30 GB              |
| swissprot  | SwissProt functional annotation database          | <10 GB              |
| metacyc    | MetaCyc functional annotation database            | <10 GB              |
| cazy       | CAZymes functional annotation database            | <10 GB              |




