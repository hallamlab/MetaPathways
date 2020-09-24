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

    ``$ virtualenv -p /usr/bin/python3 mynewenv``

Activate the new virtual environment by running 

   ``$ source mynewenv/bin/activate``

Deactivate If you want to switch projects or otherwise leave your virtual environment, simply run:

   ``$ deactivate``

pip install MetaPathways
========================
Install MetaPathways by running:

  ``$ pip3 install metapathways``

To make sure MetaPathways is installed type

  ``$ MetaPathways --version``

which, if MetaPathways, is properly installed, will print a version number. For example

::

  MetaPathways: Version 3.5.0


Install Binaries
================

Next we install ``trnascan-1.4``, ``rpkm``, ``prodigal``, ``FAST`` and ``bwa``

Download the source code as
  
  ``$ wget --version``

untar the files, make and install

::
   tar -zxvf 

   cd 
   make 
   sudo make install
   sudo make uninstall


ncbi-blast


Reference Sequences
===================

Create the reference folders as

::

   MetaPathways_DBs/
   ├── functional
   │   ├── formatted
   ├── functional_categories
   └── taxonomic
       └── formatted

`CreateDBFolderStructure  -d  <name of a folder path>`
::

   MetaPathways_DBs/
   ├── functional
   │   ├── formatted
   ├── functional_categories
   │   ├── ncbi.map
   │   ├── ncbi_taxonomy_tree.txt
   └── taxonomic
       └── formatted

Functional Reference 
++++++++++++++++++++


Taxonomic Reference 
+++++++++++++++++++


Formatting Reference Sequences
++++++++++++++++++++++++++++++


FAST
----


BLAST
-----


