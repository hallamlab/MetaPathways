Overview 
********
.. figure:: http://i.imgur.com/knn8bBb.png
    :width: 250px
    :align: center
    :height: 200px
    :alt: alternate text
    :figclass: align-center

MetaPathways [CIT2002]_ is a meta'omic analysis pipeline for the annotation and analysis for environmental sequence information.
MetaPathways include metagenomic or metatranscriptomic sequence data in one of several file formats 
(.fasta, .gff, or .gbk). The pipeline consists of five operational stages including 

Pipeline Overview
~~~~~~~~~~~~~~~~~

MetaPathways is composed of five general stages, encompassing a number of analytical or data handling steps **(Figure 1)**:

#. **QC and ORF Prediction**: Here MetaPathways performs basic quality control (QC) including removing 
duplicate sequences and sequence trimming. Open Reading Frame (ORF) prediction is then performed on the 
QC'ed sequences using Prodigal [[2](#references-1)] or GeneMark [[3](#references-1)]. The final translated 
ORFs are now also trimmed according to a user-defined setting. 

    * MetaPathways steps: `PREPROCESS INPUT`, `ORF PREDICTION`,  and `FILTER AMINOS`


#. **Functional and Taxonomic Annotation**: Using seed-and-extend homology search algorithms (B)LAST 
[[4,5](#references-1)], MetaPathways can be used to conduct searches against functional and taxonomic databases. 

    * MetaPathways steps: `FUNC SEARCH`, `PARSE FUNC SEARCH`, `SCAN rRNA`, and `ANNOTATE ORFS`

#. **Analyses**: After sequence annotation, MetaPathways performs further taxonomic analyses including the
    [Lowest Common Ancestor (LCA) algorithm](http://ab.inf.uni-tuebingen.de/software/megan/) 
    [[6](#references-1)] and [tRNA Scan](http://lowelab.ucsc.edu/tRNAscan-SE/) [[7](#references-1)], and prepares detected annotations for environmental Pathway/Genome database (ePGDB) creation via Pathway Tools.

    * MetaPathways Steps: `PATHOLOGIC INPUT`, `CREATE ANNOT REPORTS`, and `COMPUTE RPKM`.

#. **ePGDB Creation**: MetaPathways then predicts [MetaCyc pathways](http://www.metacyc.com)
using the Pathway Tools software (http://brg.ai.sri.com/ptools/)
and its pathway prediction algorithm PathoLogic [[10](#references-1)], resulting in 
the creation of an environmental Pathway/Genome Database (ePGDB), an integrative data structure of 
sequences, genes, pathways, and literature annotations for integrative interpretation.

    * MetaPathways Steps: `BUILD ePGDB`

#. **Pathway Export**: Here MetaCyc pathways or reactions are exported in a tabular format for downstream 
    analysis. *As of the v2.5 release, MetaPathways will perform this step automatically.*

    * MetaPathways Steps: `BUILD ePGDB`

.. figure:: http://i.imgur.com/HOacG2l.png
    :align: center
    :figclass: align-center

Ouptut Format
~~~~~~~~~~~~~


Visualizing Output
~~~~~~~~~~~~~~~~~~
.. [CIT2002] K. M. Konwar, N. W. Hanson, A. P. Pagé, S. J. Hallam, MetaPathways: a modular pipeline for constructing pathway/genome databases from environmental sequence information. BMC Bioinformatics 14, 202 (2013)  http://www.biomedcentral.com/1471-2105/14/202
