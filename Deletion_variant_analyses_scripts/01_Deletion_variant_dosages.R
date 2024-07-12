###############################################################################
## Genotype data to dosage file
###############################################################################
### General information

# Prerequisites:
# 1) folder data/Deletion_variants
# 2) data/Your_genotype_data in PLINK format (.bed, .bim, .fam)
###############################################################################

library(tidyverse)
library(data.table)

###############################################################################
# Create a list of deletion-tagging variants
# data/Deletion_variants/Deletion_tagging_variants.txt file containing one 
# column listing deletion-tagging variant IDs matching Your_genotype_data. The
# original deletion-tagging variants are presented in
# Steers NJ, Li Y, Drace Z, et al. Genomic Mismatch at LIMS1 Locus and Kidney 
# Allograft Rejection. N Engl J Med. 2019;380(20):1918-1928. 

###############################################################################

# Extract deletion-tagging variants from Your_genotype_data
system(command = paste0("plink --bfile data/Your_genotype_data --extract data/Deletion_variants/Deletion_tagging_variants.txt --make-bed --out data/Deletion_variants/Liver_dels_variants"))

# Create dosage file containing deletion-tagging variants
system(command = paste0("plink --bfile data/Deletion_variants/Liver_dels_variants --recodeA --out data/Deletion_variants/Liver_dels_dosage" ))

###############################################################################