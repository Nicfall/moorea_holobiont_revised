# 2brad_moorea
Scripts for analyzing my 2bRAD data from corals in French Polynesia (modified from Dr. Mikhail Matz's scripts at https://github.com/z0on/2bRAD_denovo)</br>

Head straight to the 'walkthrough' files for each part for instructions, except for the pop structure folder which is called 'part3_popstructure.Rmd'

Part 1: Genome</br>
Placing & indexing reference genome - if you don't have a reference genome you'll have to follow Misha's instructions at https://github.com/z0on/2bRAD_denovo, the rest of my instructions will also be missing any 'de novo' specific commands

Part 2: Reads</br>
Filtering & trimming raw reads</br>
Mapping reads to reference genome to create .bam files

Part 3: Population structure</br>
First, use angsd + ibs results to identify technical replicates & clones. Then use vcftools to find average site depth per sample. Then re-run angsd without worst samples & without technical replicates/clones. Explore population structure in PCA & K plot.

Part 4: FST outliers</br>
Using .vcf file to find FST outliers using Bayescan & outFLANK

Part 5: Heterozygosity </br>