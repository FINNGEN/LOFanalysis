# LOFanalysis

We want to compare two different approaches to LOF:

1) our "in-house" method where we  generate VCFs based on a custom merging of LOF variants. In this case, the genotype will be the probabily of carrying at least one LOF variants, based in turn on GPs.

2) the default genie "maps" that merge variants based on allele counts.



1 pheno   1 cpu 
4 pheno   4 cpu 
16 pheno  8 cpu
16 pheno 16 cpu
4 phenos  2 cpu 