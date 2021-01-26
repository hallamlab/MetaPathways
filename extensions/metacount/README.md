MetaCount

Install

git clone the project

make

example commands
1.  ./metacount  --gff data/lagoon-sample2.unannot.gff \
  --sam lagoon-sample2_batchr1.sam \
  --sam lagoon-sample2.sam \
  --sam lagoon-sample2_a.sam \
  --stats-type ALL

example output
#ORFID	COUNT	TPM	RPKM
0_0	26.00	3952.95	4116.19
100_0	47.00	6138.23	6391.70
100_1	26.00	3343.24	3481.30
101_0	93.50	4772.32	4969.39
...........

2.  ./metacount  --gff data/lagoon-sample2.unannot.gff \
  --sam lagoon-sample2_batchr1.sam \
  --sam lagoon-sample2.sam \
  --sam lagoon-sample2_a.sam \
  --stats-type TPM

#ORFID	TPM
0_0	4116.19
100_0	6391.70
100_1	3481.30
101_0	4969.39
......

3.  ./metacount  --gff data/lagoon-sample2.unannot.gff \
  --sam lagoon-sample2_batchr1.sam \
  --sam lagoon-sample2.sam \
  --sam lagoon-sample2_a.sam \
  --stats-type RPKM

#ORFID	RPKM
0_0	3952.95
100_0	6138.23
100_1	3343.24
101_0	4772.32
101_1	5893.60

4.  ./metacount  --gff data/lagoon-sample2.unannot.gff \
  --sam lagoon-sample2_batchr1.sam \
  --sam lagoon-sample2.sam \
  --sam lagoon-sample2_a.sam \
  --stats-type COUNT


example output



C++ implementation of reads per kilobase mapped statistic. Functional analysis of de novo assembled environmental sequence
information is impeded by the lack of quantitative ORF annotations. ORF counts are affected by both sequencing depth and ORF
length, longer ORFs naturally encompass more reads, making quantitative comparisons between samples difficult. To resolve this, we have implemented a bwa-based version of the RPKM. Intuitively RPKM is a simple proportion of the number of reads mapped to a sequence section, normalized for sequencing depth and ORF length. This tools is now a part of the MetaPathways pipeline.

References:

1. Kishori M. Konwar, Niels W. Hanson, Maya P. Bhatia, Dongjae Kim, Shang-Ju Wu, Aria S. Hahn, Connor Morgan-Lang, Hiu Kan Cheung, Steven J. Hallam. MetaPathways v2.5: quantitative functional, taxonomic and usability improvements. Bioinformatics 31(20), pp. 3345-3347 (2015).

2. Kishori M. Konwar, Niels W. Hanson, Antoine P. Pagé, Steven J. Hallam. MetaPathways: a modular pipeline for constructing pathway/genome databases from environmental sequence information. BMC Bioinformatics 14: 202 (2013).
