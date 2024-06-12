###############################################################################
## Genotype data to dosage file
###############################################################################

library(tidyverse)
library(data.table)

###############################################################################
# Extract deletion-tagging variants from the genotype files
system(command = paste0("plink --bfile data/Liver_genotype_file --extract data/Deletion_tagging_variants.txt --make-bed --out data/Liver_dels_variants"))

# Create dosage file containing 40 deletion-tagging variants
system(command = paste0("plink --bfile data/Liver_dels_variants --recodeA --out data/Liver_dels_dosage" ))