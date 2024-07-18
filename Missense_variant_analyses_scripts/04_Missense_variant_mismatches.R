###############################################################################
### Calculation of recipient-donor missense variant mismatches
###############################################################################
### General information

# Prerequisites:
# 1. Run script 03_Missense_variant_dosages.R
  # data/Missense_variants/Your_genotype_data_all_missense_dosage.raw
  # data/Missense_variants/Your_genotype_data_transm_secr_missense_dosage.raw
  # data/Missense_variants/Your_genotype_data_transm_missense_dosage.raw
  # data/Missense_variants/Your_genotype_data_liver_missense_dosage.raw
# 2.data/Your_clinical_data.txt

###############################################################################

library(data.table)
library(tidyverse)

###############################################################################
# Import the clinical data file 
Covariates <- read_table("data/Your_clinical_data")

# Structure specification example
str(Covariates)
# spec_tbl_df [666 × 21] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ Pair_number               : num [1:666] 539 540 542 543 541 544 545 546 547 548 ...
# $ R_pseudo                  : chr [1:666] "R_pseudo1" "R_pseudo2" "R_pseudo3" "R_pseudo4" ...
# $ D_pseudo                  : chr [1:666] "D_pseudo1" "D_pseudo2" "D_pseudo3" "D_pseudo4" ...
# $ R_sex                     : num [1:666] 2 1 1 2 2 2 1 1 2 1 ...
# $ R_age                     : num [1:666] 50 17 51 60 30 49 47 35 31 72 ...
# $ D_age                     : num [1:666] 53 33 74 59 24 59 53 62 47 67 ...
# $ D_sex                     : num [1:666] 2 1 1 1 1 2 1 2 2 2 ...
# $ Cold_ischemia_time_minutes: num [1:666] 320 301 715 452 268 246 350 248 266 234 ...
# $ AR_status                 : num [1:666] 0 1 0 0 1 0 1 1 0 0 ...
# $ Eplets_total_HLAI         : num [1:666] 10 19 19 19 10 11 11 29 13 26 ...
# $ Eplets_total_HLAII        : num [1:666] 4 7 8 4 12 15 9 23 10 11 ...
# $ AR_Cox_time               : num [1:666] 86.63 0.197 83.605 83.605 0.723 ...
# $ Transplantation_year      : num [1:666] 2012 2013 2014 2015 2014 ...
# $ Autoimmune_status         : num [1:666] 1 1 0 0 1 0 1 0 1 0 ...
# $ CNI_type_initial_status   : num [1:666] 1 1 1 1 1 1 1 1 1 1 ...
# $ LR_Cox_time               : num [1:666] 84.6 86 85 82.6 12.5 ...
# $ Late_rejection_status     : num [1:666] 0 0 0 0 1 0 0 0 0 0 ...
# $ Graft_loss_status         : num [1:666] 0 0 0 0 0 0 1 0 0 0 ...
# $ graft_loss_months         : num [1:666] 86.6 85 83.6 83.6 83.6 ...
# $ Death_status              : num [1:666] 0 0 0 0 0 0 0 0 0 0 ...
# $ Death_Cox_time_months     : num [1:666] 82.6 87 84.6 82.6 81.6 ...

### Create new variable for mismatch calculations 
R_D_pairs <- Covariates %>% select(Pair_number, R_pseudo, D_pseudo)

###############################################################################
### Import the missense dosage files for each protein group
dosage_all <- fread("data/Missense_variants/Your_genotype_data_all_missense_dosage.raw")

dosage_transm_secr <- fread("data/Missense_variants/Your_genotype_data_transm_secr_missense_dosage.raw")

dosage_transm <- fread("data/Missense_variants/Your_genotype_data_transm_missense_dosage.raw")

dosage_liver <- fread("data/Missense_variants/Your_genotype_data_liver_missense_dosage.raw")

###############################################################################
### Inner join imported dosage files with R_D_pairs

# All proteins
R_paired_dos_all <- inner_join(R_D_pairs, dosage_all,
                                  by = c("R_pseudo" = "IID"))

D_paired_dos_all <- inner_join(R_D_pairs, dosage_all,
                                  by = c("D_pseudo" = "IID")) 
# Transmembrane and secretory
R_paired_dos_transm_secr <- inner_join(R_D_pairs, dosage_transm_secr,
                                     by = c("R_pseudo" = "IID"))

D_paired_dos_transm_secr <- inner_join(R_D_pairs, dosage_transm_secr,
                                     by = c("D_pseudo" = "IID"))

# Transmembrane only
R_paired_dos_transm <- inner_join(R_D_pairs,
                                 dosage_transm,
                                 by = c("R_pseudo" = "IID"))

D_paired_dos_transm <- inner_join(R_D_pairs,
                                 dosage_transm,
                                 by = c("D_pseudo" = "IID"))

# Liver-related
R_paired_dos_liver <- inner_join(R_D_pairs,
                                 dosage_liver,
                                 by = c("R_pseudo" = "IID"))

D_paired_dos_liver <- inner_join(R_D_pairs,
                                 dosage_liver,
                                 by = c("D_pseudo" = "IID"))

###############################################################################
### Check that the pair order is the same in both recipient and donor data 
### frames
identical(R_paired_dos_all$Pair_number, 
          D_paired_dos_all$Pair_number)

identical(R_paired_dos_transm_secr$Pair_number, 
          D_paired_dos_transm_secr$Pair_number)

identical(R_paired_dos_transm$Pair_number, 
          D_paired_dos_transm$Pair_number)

identical(R_paired_dos_liver$Pair_number, 
          D_paired_dos_liver$Pair_number)

###############################################################################
### Calculate the mismatch sums in each protein group

## All proteins
Mismatch <- function(R,D) {
  ifelse((D>0 & R==0) | (D<2 & R==2), T, F) 
}

Mm_all_result <- sapply(2:ncol(R_paired_dos_all), 
                        function(i) {Mismatch(R_paired_dos_all[,i], 
                                              D_paired_dos_all[,i])})

Mm_all_result_df <- data.frame(Pair_number=R_paired_dos_all$Pair_number, 
                               Mm_all=rowSums(Mm_all_result))

# Match covariates and the number of mismatches
R_covariates_mm_all <- inner_join(Covariates, Mm_all_result_df,
                                  by = "Pair_number")

## Transmembrane and secreted
Mm_transm_secr_res <- sapply(2:ncol(R_paired_dos_transm_secr), 
                               function(i) {Mismatch(R_paired_dos_transm_secr[,i], 
                                                     D_paired_dos_transm_secr[,i])})

Mm_transm_secr_res_df <- data.frame(Pair_number=R_paired_dos_transm_secr$Pair_number, 
                                   Mm_trans_secr=rowSums(Mm_transm_secr_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_transm_secr <- inner_join(R_covariates_mm_all, 
                                         Mm_transm_secr_res_df,
                                         by = "Pair_number")

## Transmembrane only
Mm_transm_res <- sapply(2:ncol(R_paired_dos_transm), 
                                  function(i) {Mismatch(R_paired_dos_transm[,i], 
                                                        D_paired_dos_transm[,i])})

Mm_transm_res_df <- data.frame(Pair_number=R_paired_dos_transm$Pair_number,
                                         Mm_transmemb=rowSums(Mm_transm_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_transm <- inner_join(R_covariates_mm_transm_secr,
                                            Mm_transm_res_df,
                                            by = "Pair_number")

## Liver-related
Mm_liver_res <- sapply(2:ncol(R_paired_dos_liver), 
                          function(i) {Mismatch(R_paired_dos_liver[,i], 
                                                D_paired_dos_liver[,i])})

Mm_liver_res_df <- data.frame(Pair_number=R_paired_dos_liver$Pair_number,
                                 Mm_liver=rowSums(Mm_liver_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_liver <- inner_join(R_covariates_mm_transm,
                                    Mm_liver_res_df,
                                    by = "Pair_number")

###############################################################################
# Write out the final data file
write.table(R_covariates_mm_liver,
            file = "data/Missense_variants/R_covariates_mm_liver.txt",
            row.names = F, col.names = T, quote = F)

###############################################################################
