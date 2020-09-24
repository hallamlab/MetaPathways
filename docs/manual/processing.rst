Running MetaPathways
********************

Input
=====
MetaPathways inputs are fasta files provided in an input folder. The file names must end with 
a `.fasta` or `.fas`. These fasta files contains the contigs or DNA sequences from assembling.

Parameter File 
==============

The parameter file must indicate the setting for any MetaPathways run. An example paramter file 
can be downloaded as
::

 $ wget  https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/text/template_param.txt

Below we describe the settings in the parameter file. 

Run
===

 As an illustration we donwload a small input file `testsample1.fasta`
 in a folder named `mp_input` and we want the output in a folder names `mp_output`

::

 $ mkdir mp_input
 $ cd mp_input
 $ wget https://github.com/kishori82/MetaPathways_Python.3.0/raw/kmk-develop/data/testdata/testsample1.fasta
 $ cd ..

Now we kick off a run  as
::

  $ MetaPathways --input mp_input --output mp_output -p template_param.txt -d ~/MetaPathways_DBs/
 





