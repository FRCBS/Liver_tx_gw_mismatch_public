###############################################################################
### Creating missense variant dosage files for four different missense variant
### types (chr 1-22 wo MHC region):
### 1) For the whole data set 
### 2) For transmembrane and secreted proteins
### 3) For transmembrane proteins only
### 4) For liver-specific proteins
###############################################################################
### Prerequisites:
# 1. Run scripts 01_Variants_for_VEP.R and 02_VEP_missense_variants
# 2. PLINK 1.90 https://www.cog-genomics.org/plink/
# 3. Your genotype data "data/Your_genotype_data" in PLINK format (.bed, .bim, 
#    .fam)
# 4. data/Missense_variants/missense_variants.txt file 

###############################################################################
### 1) All missensese variants
# Match data from VEP search data/Missense_variants/missense_variants.txt
# to your genotype data (variant IDs).
# Create final file "data/Missense_variants/missense_FINAL" containing one column
# listing the matched missense variant IDs

### 2) Transmembrane and secreted proteins
# Search transcripts for both transmembrane and secreted proteins from Uniprot database
# https://www.uniprot.org/ using the following string:
# annotation:(type:transmem) OR locations:(location:"Secreted [SL-0243]") OR keyword:"Transmembrane [KW-0812]") AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"
# Match the results to missense variants and your genotype data.
# Create final file "data/Missense_variants/transm_secr_missense_variants" containing one column
# listing the matched transmembrane and secreted protein missense variant IDs

### 3) Transmembrane proteins
# Search transcripts for both transmembrane proteins from Uniprot database
# https://www.uniprot.org/ using the following string:
# annotation:(type:transmem) OR keyword:"Transmembrane [KW-0812]") AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"
# Match the results to missense variants and your genotype data.
# Create final file "data/Missense_variants/transm_missense_variants" containing one column
# listing the matched transmembrane protein missense variant IDs

### 4) Liver-related proteins
# Search liver only proteins from Uniprot database
# https://www.uniprot.org/ using the following string:
# annotation:(type:"tissue specificity" liver) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"
# Match the results to missense variants and your genotype data.
# Create final file "data/Missense_variants/liver_variants" containing one column
# listing the matched liver protein missense variant IDs

### Create "MHC_region.txt"
# The following script requires a text file "MHC_region.txt" with chromosome,
# starting position and ending position of MHC region.

# 6 28510120 33480577 R1

###############################################################################

library(tidyverse)
library(data.table)

###############################################################################
### Genotype dosage files

### 1) Dosage file for all missense variants (chr 23 and MHC region not included)

# Exclude MHC region and using only chromosomes 1-22 
system(command = paste0("plink --bfile data/Liver_cohort --chr 1-22 --exclude data/Missense_variants/MHC_region.txt --range --make-bed --out data/Missense_variants/Liver_Tx_cohort_chr1_22_without_MHC"))

# Extract the missense variants from the filtered data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/missense_final.txt --make-bed --out data/Missense_variants/LIVER_missense_variants_no_X_MHC"))
# n=28225

#Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/LIVER_missense_variants_no_X_MHC --recodeA --out data/Missense_variants/Liver_missense_dosage_wo_X_MHC"))

### 2) Dosage file for transmembrane and secreted proteins

# Extract the transmembrane and secreted missense variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_secr_transm_variants_final.txt --make-bed --out data/Missense_variants/Liver_missense_transmemb_secr"))
# n=9337

# Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_transmemb_secr --recodeA --out data/Missense_variants/Liver_missense_transmemb_secr_dosage"))

### 3) Dosage file for transmembrane-only proteins
# Extract the transmembrane variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_transmemb_variants_final.txt --make-bed --out data/Missense_variants/Liver_missense_transmemb"))
# n=6963

# Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_transmemb --recodeA --out data/Missense_variants/Liver_missense_transmembrane_dosage"))

### 4) Dosage file for liver-spesific proteins
# Extract the liver variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_liver_specific_variants_final.txt --make-bed --out data/Missense_variants/Liver_missense_liver_spesific"))
# n=3199

# Creating a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_liver_spesific --recodeA --out data/Missense_variants/Liver_missense_liver_spesific_dosage"))
