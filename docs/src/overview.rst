Overview 
********

MetaPathways [MP2013]_ is a meta'omic analysis pipeline for the annotation and analysis for environmental sequence information.
MetaPathways include metagenomic or metatranscriptomic sequence data in one of several file formats 
(.fasta, .gff, or .gbk). The pipeline consists of five operational stages including 

.. figure:: static/glyph.png
   :align: center
   :alt: alternate text
   :figclass: align-center

Pipeline
~~~~~~~~

MetaPathways is composed of four general stages, encompassing a number of analytical or data handling steps **(Figure 1)**:

.. |nbsp| unicode:: 0xA0 
   :trim:

|nbsp|
	
#. **Quality Control**: 
   Basic quality control (QC) is performed with includes filtering out sequences below a set
   length threshold (default 180bp). At this stage any duplicate sequences are removed (optional).
    
#. **Feature Prediction**:
   Several sequence features can be predicted on the QC'ed contigs. Open-reading frames
   (ORFs) are predicted by default and (optionally) ribosomal subunits (rRNAs) and
   transfer RNAs (tRNAs) can be predicted. To improve the runtime and efficency, Prodigal
   [PRODIGAL2010]_ is run through a parallel version (pProdigal) [PPRODIGAL2022]_ and
   tRNAscan-SE [TRNASCANSE2021]_ is run using a wrapper script that allow for more
   efficient multi-threading. BARRNAP [BARRNAP2019]_ is used for the prediction of rRNAs
   including: 16S, 23S, and 5S. MetaPathways provides an overlap-aware identification so
   that users can make informed decisions about features when they overlap eachother.
   Addtionally, users can define the minimum length of ORFs to keep for downstream analysis.
   
#. **Functional Annotation**:
   Using a seed-and-extend homology search algorithm, either BLAST [BLAST]_ or FAST [FAST]_,
   users can conduct searches against both functional and taxonomic (optional) databases. 
   Currently supported databases include: Uniprot SwissProt [SWISSPROT]_, Uniprot UniRef90
   [UNIREF90]_, MetaCyc [METACYC]_, and CAZymes [CAZy]_. However, users can create custom
   databases for any preferred databases. Optionally, reads can be used to calculate abundance
   information at both the contig and ORF-level.
      
#. **Pathway Inference**:
   MetaPathways then predicts `MetaCyc pathways <http://www.metacyc.com>`_ using the
   `Pathway Tools software <http://brg.ai.sri.com/ptools/>`_ and its pathway prediction
   algorithm PathoLogic [KARP11]_, resulting in the creation of a community-level environmental
   Pathway/Genome Database (ePGDB), an integrative data structure of sequences, genes, pathways,
   and literature annotations for integrative interpretation. Optionally, if metagenome-assembled
   genomes (MAGs) are available for the metagenome, these MAGs can be used to create
   population-level ePGDBs. MetaCyc pathways are exported in a tabular format for
   downstream analysis.



Bibliography
~~~~~~~~~~~~

Please see the following `Zotero Library <https://www.zotero.org/groups/5284390/bcb2/library>`_ for a bibliography.

..
   .. [MP2013] K. M. Konwar, N. W. Hanson, A. P. Pagé, S. J. Hallam, MetaPathways: a modular 
      pipeline for constructing pathway/genome databases from environmental sequence information. 
      BMC Bioinformatics 14, 202 (2013)  http://www.biomedcentral.com/1471-2105/14/202

   .. [PRODIGAL2010] D. Hyatt et al., Prodigal: prokaryotic gene recognition and translation 
      initiation site identification. BMC Bioinformatics 11, 119 (2010).

   .. [PPRODIGAL2022] Jaenicke, S. (2022). pprodigal (Version 1.0.1) [Software].
      Available from https://pypi.org/project/pprodigal/

   .. [TRNASCANSE2021] Chan, P. P., Lin, B. Y., Mak, A. J., & Lowe, T. M. (2021).
      tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes.
      Nucleic Acids Research, 49(16), 9077-9096. https://doi.org/10.1093/nar/gkab688

   .. [BARRNAP2019] Seemann, T. (2019). barrnap (Version 0.9) [Software]. Available from https://github.com/tseemann/barrnap




   fast
   swissprot
   uniref90
   metacyc
   cazymes


   .. [GeneMark12] D. Hyatt, P. F. LoCascio, L. J. Hauser, E. C. Uberbacher, Gene and translation initiation site prediction in metagenomic sequences. Bioinformatics 28, 2223–2230 (2012).

   .. [BLAST90] S. F. Altschul, W. Gish, W. Miller, E. W. Myers, D. J. Lipman, Basic local alignment search tool. J Mol Biol 215, 403–410 (1990).
   .. [LAST11]  S. M. Kiełbasa, R. Wan, K. Sato, P. Horton, M. C. Frith, Adaptive seeds tame genomic sequence comparison. Genome Res 21, 487–493 (2011).

   .. [MEGAN07] D. H. Huson, A. F. Auch, J. Qi, S. C. Schuster, MEGAN analysis of metagenomic data. Genome Res 17, 377–386 (2007).
   .. [TRNASCAN97] T. M. Lowe, S. R. Eddy, tRNAscan-SE: a program for improved detection of transfer RNA genes in genomic sequence. Nucleic Acids Research 25, 0955–0964 (1997).

   ..   R. Caspi et al., The MetaCyc database of metabolic pathways and enzymes and the BioCyc collection of pathway/genome databases. Nucleic Acids Research 38, D473–D479 (2009).
      P. D. Karp, S. Paley, P. Romero, The pathway tools software. Bioinformatics 18, S225–S232 (2002).

   .. [KARP11] P. D. Karp, M. Latendresse, R. Caspi, The pathway tools pathway prediction algorithm. Stand Genomic Sci 5, 424–429 (2011).

   ..  K. D. Pruitt, T. Tatusova, D. R. Maglott, NCBI reference sequences (RefSeq): a curated non-redundant sequence database of genomes, transcripts and proteins. Nucleic Acids Research 35, D61–5 (2007).
   ..  H. Li, R. Durbin, Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics 26, 589–595 (2010).
      R. L. Tatusov et al., The COG database: an updated version includes eukaryotes. BMC Bioinformatics 4, 41 (2003).
      M. Kanehisa, S. Goto, KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Research 28, 27–30 (2000).
      F. Meyer et al., The metagenomics RAST server - a public resource for the automatic phylogenetic and functional analysis of metagenomes. BMC Bioinformatics 9, 386 (2008).
      R. K. Aziz et al., SEED servers: high-performance access to the SEED genomes, annotations, and metabolic models. PLoS ONE 7, e48053 (2012).
      B. L. Cantarel et al., The Carbohydrate-Active EnZymes database (CAZy): an expert resource for Glycogenomics. Nucleic Acids Research 37, D233–D238 (2009).
