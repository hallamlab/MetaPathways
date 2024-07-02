Installation
************

Metapathways can be installed with conda or obtained as a container.
The fastest install is via Mamba (recommended) or by pulling down the container.

.. note::
   
   Metapathways has been built and tested on Linux OS.
   Lightly tested on Mac OS.
   And not currently supported on native Windows.
   

Mamba (recommended)
-------------------

.. code-block:: bash

   mamba create -n metapathways_env -c hallamlab -c bioconda -c conda-forge metapathways

Conda
-----

.. code-block:: bash

   conda create -n metapathways_env -c hallamlab -c bioconda -c conda-forge metapathways

Docker
---------

.. code-block:: bash
   
   docker pull quay.io/hallamlab/metapathways
   
