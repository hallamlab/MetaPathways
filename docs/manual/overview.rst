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


#. **QC and ORF Prediction**: Here MetaPathways performs basic quality control (QC) including removing duplicate 
   sequences and sequence trimming. Open Reading Frame (ORF) prediction is then performed on the QC'ed sequences 
   using Prodigal [PRODIGAL2010]_ or GeneMark [GeneMark12]_. The final translated ORFs are 
   now also trimmed according to a user-defined setting. 

   * MetaPathways steps: `PREPROCESS INPUT`, `ORF PREDICTION`,  and `FILTER AMINOS`


#. **Functional and Taxonomic Annotation**: Using seed-and-extend homology search algorithms (B)LAST 
   [BLAST90]_, [LAST11]_, MetaPathways can be used to conduct searches against functional and taxonomic databases. 

   * MetaPathways steps: `FUNC SEARCH`, `PARSE FUNC SEARCH`, `SCAN rRNA`, and `ANNOTATE ORFS`

#. **Analyses**: After sequence annotation, MetaPathways performs further taxonomic analyses including 
   the `Lowest Common Ancestor (LCA) 
   <http://ab.inf.uni-tuebingen.de/software/megan/>`_ algorithm 
   [MEGAN07]_ and `tRNA Scan <http://lowelab.ucsc.edu/tRNAscan-SE/>`_ [TRNASCAN97]_, and 
   prepares detected annotations for environmental Pathway/Genome database (ePGDB) creation via Pathway Tools.

   * MetaPathways Steps: `PATHOLOGIC INPUT`, `CREATE ANNOT REPORTS`, and `COMPUTE RPKM`.

#. **ePGDB Creation**: MetaPathways then predicts `MetaCyc pathways
   <http://www.metacyc.com>`_ using 
   the `Pathway Tools software
   <http://brg.ai.sri.com/ptools/>`_ 
   and its pathway prediction algorithm 
   PathoLogic [KARP11]_, resulting in the creation of an environmental Pathway/Genome 
   Database (ePGDB), an integrative data structure of sequences, genes, pathways, and literature 
   annotations for integrative interpretation.

   * MetaPathways Steps: `BUILD ePGDB`

#. **Pathway Export**: Here MetaCyc pathways or reactions are exported in a tabular format for downstream 
   analysis. *As of the v2.5 release, MetaPathways will perform this step automatically.*

   * MetaPathways Steps: `BUILD ePGDB`

.. figure:: http://i.imgur.com/HOacG2l.png
    :align: center
    :figclass: align-center
   
Output Format
~~~~~~~~~~~~~


Visualizing Output
~~~~~~~~~~~~~~~~~~
.. [CIT2002] K. M. Konwar, N. W. Hanson, A. P. Pagé, S. J. Hallam, MetaPathways: a modular 
   pipeline for constructing pathway/genome databases from environmental sequence information. 
   BMC Bioinformatics 14, 202 (2013)  http://www.biomedcentral.com/1471-2105/14/202

.. [PRODIGAL2010] D. Hyatt et al., Prodigal: prokaryotic gene recognition and translation 
   initiation site identification. BMC Bioinformatics 11, 119 (2010).

.. [GeneMark12] D. Hyatt, P. F. LoCascio, L. J. Hauser, E. C. Uberbacher, Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28, 2223–2230 (2012).

.. [BLAST90] S. F. Altschul, W. Gish, W. Miller, E. W. Myers, D. J. Lipman, Basic local alignment search tool. J Mol Biol 215, 403–410 (1990).
.. [LAST11]  S. M. Kiełbasa, R. Wan, K. Sato, P. Horton, M. C. Frith, Adaptive seeds tame genomic sequence comparison. Genome Res 21, 487–493 (2011).

.. [MEGAN07] D. H. Huson, A. F. Auch, J. Qi, S. C. Schuster, MEGAN analysis of metagenomic data. Genome Res 17, 377–386 (2007).
.. [TRNASCAN97] T. M. Lowe, S. R. Eddy, tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence. Nucleic Acids Research 25, 0955–0964 (1997).

..   R. Caspi et al., The MetaCyc database of metabolic pathways and enzymes and the BioCyc collection of pathway/genome databases. Nucleic Acids Research 38, D473–D479 (2009).
   P. D. Karp, S. Paley, P. Romero, The pathway tools software. Bioinformatics 18, S225–S232 (2002).

.. [KARP11] P. D. Karp, M. Latendresse, R. Caspi, The pathway tools pathway prediction algorithm. Stand Genomic Sci 5, 424–429 (2011).

..  K. D. Pruitt, T. Tatusova, D. R. Maglott, NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Research 35, D61–5 (2007).
   H. Li, R. Durbin, Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics 26, 589–595 (2010).
   R. L. Tatusov et al., The COG database: an updated version includes eukaryotes. BMC Bioinformatics 4, 41 (2003).
   M. Kanehisa, S. Goto, KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Research 28, 27–30 (2000).
   F. Meyer et al., The metagenomics RAST server - a public resource for the automatic phylogenetic and functional analysis of metagenomes. BMC Bioinformatics 9, 386 (2008).
   R. K. Aziz et al., SEED servers: high-performance access to the SEED genomes, annotations, and metabolic models. PLoS ONE 7, e48053 (2012).
   B. L. Cantarel et al., The Carbohydrate-Active EnZymes database (CAZy): an expert resource for Glycogenomics. Nucleic Acids Research 37, D233–D238 (2009).
