################################################################################
### Variant list for VEP annotation
################################################################################

library(tidyverse)

################################################################################
### Variants in default VEP input format
### https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#input

### Prerequisites:
# 1) Folders data, Missense_variant_analyses_scripts, 
#    Deletion_variant_analyses_scripts, results
# 2) Foder data/Missense_variants
# 3) Your genotype file "Your_genotype_data" in PLINK format (.bed, .bim, .fam) 
# in data folder


################################################################################
# Import corresponding PLINK bim file
Liver_Tx_cohort <- read_table2("data/Your_genotype_data.bim", 
                               col_names = F)
  # X1, Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
  # X2, Variant identifier
  # X3, Position in morgans or centimorgans (safe to use dummy value of '0')
  # X4, Base-pair coordinate (1-based; limited to 231-2)
  # X5 Allele 1 (corresponding to clear bits in .bed; usually minor)
  # X6, Allele 2 (corresponding to set bits in .bed; usually major)

# Genotype in VEP format ref/alt
Liver_Tx_cohort$allele <- unite(Liver_Tx_cohort, "allele", X6, X5, 
                            sep = "/", remove = F)$allele

# Reference allele length
Liver_Tx_cohort$ref_length <- nchar(Liver_Tx_cohort$X6)
# Reduce 1 for end position calculation
Liver_Tx_cohort$ref_length <- as.numeric(Liver_Tx_cohort$ref_length -1)

# Modify and rename columns
Liver_Tx_cohort <- select(Liver_Tx_cohort, -X2, -X3, -X5, -X6) %>% 
  rename(chr = X1, start = X4)

# Compute end position column
Liver_Tx_cohort$end <- Liver_Tx_cohort$start + Liver_Tx_cohort$ref_length

# Add strand information column
Liver_Tx_cohort$strand <- "+"

# Reorder columns
Liver_Tx_cohort <- select(Liver_Tx_cohort, chr, start, end, allele, strand)

# Structure of data
str(Liver_Tx_cohort)
# tibble [8,706,949 ? 5] (S3: tbl_df/tbl/data.frame)
# $ chr   : num [1:8706949] 1 1 1 1 1 1 1 1 1 1 ...
# $ start : num [1:8706949] 702040 722408 722583 722700 727151 ...
# $ end   : num [1:8706949] 702040 722408 722583 722700 727151 ...
# $ allele: chr [1:8706949] "C/T" "C/G" "C/T" "G/A" ...
# $ strand: chr [1:8706949] "+" "+" "+" "+" ...

# Write table out wo header
write.table(Liver_Tx_cohort, "data/Missense_variants/Variables_for_VEP",
            quote = F, col.names = F, row.names = F,
            sep = " ") 