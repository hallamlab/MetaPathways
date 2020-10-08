.. MetaPathways documentation master file, created by
   sphinx-quickstart on Wed Sep 23 13:09:29 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MetaPathways's documentation!
========================================

.. toctree::
   :maxdepth: 2
   :caption: Documentation

   manual/overview
   manual/install
   manual/processing

.. toctree::
   :maxdepth: 2
   :caption: Data Exploration

   manual/phandi_overview


.. toctree::
   :maxdepth: 2
   :caption: API References


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. MetaPathways -i tests/data/lagoon-sample/input/  -o mp_output/ -s lagoon-sample -p template_param.txt  -d ~/MetaPathways_DBs/ -v
   pytest --import-mode importlib -v
