MetaPathways: A modular pipeline for constructing Pathway/Genome Databases from environmental sequence information
============

* Update: MetaPathways v2.0 is [available from its GitHub repository](https://github.com/hallamlab/metapathways2).
  * Installation instructions and tutorials can be found [on the wiki](https://github.com/hallamlab/MetaPathways/wiki)
* Update: MetaPathways v3.5 is [available from its Bitbucket repository](https://bitbucket.org/BCB2/metapathways).
  * [Documentation](https://metapathways.readthedocs.io/en/latest/quick_start.html), [install via conda](https://anaconda.org/hallamlab/metapathways), [run with docker](https://quay.io/repository/hallamlab/metapathways)

### Abstract

Background: A central challenge to understanding the ecological and biogeochemical roles of microorganisms in natural and human engineered ecosystems is the reconstruction of metabolic interaction networks from environmental sequence information. The dominant paradigm in metabolic reconstruction is to assign functional annotations using BLAST. Functional annotations are then projected onto symbolic representations of metabolism in the form of KEGG pathways or SEED subsystems. 

Results: Here we present MetaPathways, an open source pipeline for pathway inference that uses the PathoLogic algorithm to construct environmental Pathway/Genome Databases (ePGDBs) compatible with the editing and navigation features of Pathway Tools. The pipeline accepts assembled or unassembled nucleotide sequences, performs quality assessment and control, predicts and annotates noncoding genes and open reading frames, and produces inputs to PathoLogic. In addition to constructing ePGDBs, MetaPathways uses MLTreeMap to build phylogenetic trees for selected taxonomic anchor and functional gene markers, converts General Feature Format (GFF) files into concatenated GenBank files for ePGDB construction based on third-party annotations and generates useful file formats including Sequin files for direct GenBank submission and gene feature tables summarizing annotations, MLTreeMaps and ePGDB pathway coverage summaries for statistical comparisons. 

Conclusions: Metapathways provides users with a modular annotation and analysis pipeline for predicting metabolic interaction networks from environmental sequence information using an alternative to KEGG pathways and SEED subsystems mapping. It is extensible to genomic and transcriptomic datasets from a wide range of sequencing platforms, and generates useful data products for microbial community structure and function analysis. The MetaPathways software package, installation instructions, and example data can be obtained from http://hallam.microbiology.ubc.ca/MetaPathways

Keywords: Environmental Pathway/Genome Database (ePGDB), metagenome, Pathway Tools, PathoLogic, MetaCyc, microbial community, metabolism, metabolic interaction networks 

Please cite: [Konwar, Kishori M., et al. "MetaPathways: a modular pipeline for constructing pathway/genome databases from environmental sequence information." BMC bioinformatics 14.1 (2013): 202.](http://www.biomedcentral.com/1471-2105/14/202)
