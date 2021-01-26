RPKM

C++ implementation of reads per kilobase mapped statistic. Functional analysis of de novo assembled environmental sequence
information is impeded by the lack of quantitative ORF annotations. ORF counts are affected by both sequencing depth and ORF
length, longer ORFs naturally encompass more reads, making quantitative comparisons between samples difficult. To resolve this, we have implemented a bwa-based version of the RPKM. Intuitively RPKM is a simple proportion of the number of reads mapped to a sequence section, normalized for sequencing depth and ORF length. This tools is now a part of the MetaPathways pipeline.

References:

1. Kishori M. Konwar, Niels W. Hanson, Maya P. Bhatia, Dongjae Kim, Shang-Ju Wu, Aria S. Hahn, Connor Morgan-Lang, Hiu Kan Cheung, Steven J. Hallam. MetaPathways v2.5: quantitative functional, taxonomic and usability improvements. Bioinformatics 31(20), pp. 3345-3347 (2015).

2. Kishori M. Konwar, Niels W. Hanson, Antoine P. Pagé, Steven J. Hallam. MetaPathways: a modular pipeline for constructing pathway/genome databases from environmental sequence information. BMC Bioinformatics 14: 202 (2013).
