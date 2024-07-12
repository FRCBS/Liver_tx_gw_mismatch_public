###############################################################################
## Genotype data to dosage file
###############################################################################
### General information

# Prerequisites:
# 1) data/Deletion_variants directory
# 2) data/Your_genotype_data" in PLINK format (.bed, .bim, 
#    .fam)
###############################################################################

library(tidyverse)
library(data.table)

###############################################################################
# Create a list of deletion-tagging variants
# data/Deletion_variants/Deletion_tagging_variants.txt file containing a list
# of deletion tagging variant names (from LIMS paper)

# Extract deletion-tagging variants from the genotype files
system(command = paste0("plink --bfile data/Deletion_variants/Liver_genotype_file --extract data/Deletion_tagging_variants.txt --make-bed --out data/Liver_dels_variants"))

# Create dosage file containing 40 deletion-tagging variants
system(command = paste0("plink --bfile data/Liver_dels_variants --recodeA --out data/Liver_dels_dosage" ))