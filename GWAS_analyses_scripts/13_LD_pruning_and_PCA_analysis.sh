#!/bin/bash

###############################################################################
### Perform LD pruning to remove SNPs with high LD with each other
###############################################################################
## General information

# Prerequisites:
# 1) data/Your_genotype_data in PLINK format (.bed, .bim, .fam)
###############################################################################

plink2 \
   --bfile data/Your_genotype_data \
   --indep-pairwise 50 5 0.15 \
   --out data/Pruned_SNPs
   
# Perform PCA analysis
plink2 \
  --bfile data/Your_genotype_data \
  --extract data/Pruned_SNPs.prune.in \
  --pca \
  --out data/PCA/Post_PCA