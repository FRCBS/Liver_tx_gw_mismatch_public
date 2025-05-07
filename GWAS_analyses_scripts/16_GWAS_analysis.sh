#!/bin/bash
##############################################################################
### Perform logistic regression analysis for four end points
###############################################################################
## General information

# Prerequisites:
# 1) data/Your_genotype_data
# 2) Run script 15_Covariate_phenotype_creation.R
#    data/GWAS/R_covariates.txt
#    data/GWAS/D_covariates.txt
#    data/GWAS/R_pheno.txt
#    data/GWAS/D_pheno.txt
###############################################################################

# For recipients
plink2 \
   --bfile data/Your_genotype_data \
   --glm hide-covar \
   --ci 0.95 \
   --1 \
   --covar-variance-standardize \
   --covar data/GWAS/R_covariates.txt \
   --pheno data/GWAS/R_pheno.txt \
   --out results/GWAS/Logistic_regression_recipients

# For donors
plink2 \
   --bfile data/Your_genotype_data \
   --glm hide-covar \
   --ci 0.95 \
   --1 \
   --covar-variance-standardize \
   --covar data/GWAS/D_covariates.txt \
   --pheno data/GWAS/D_pheno.txt\
   --out results/GWAS/Logistic_regression_donors