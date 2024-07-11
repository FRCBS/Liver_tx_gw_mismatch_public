###############################################################################
### Recipient and donor genome-wide mismatch
###############################################################################
### General information

# Prerequisites:
# 1) Run script 03_Missense_variant_dosages.R
# 2) data/Missense_variants/Liver_missense_dosage_wo_X_MHC.raw
# 3) data/Missense_variants/Liver_missense_transmemb_secr_dosage.raw
# 4) data/Missense_variants/Liver_missense_transmembrane_dosage.raw
# 5) data/Missense_variants/Liver_missense_liver_spesific_dosage.raw
# 6) data/Clinical_data_pairs file

###############################################################################

library(data.table)
library(tidyverse)

###############################################################################
### Import the missense dosage files for each protein group
dosage_all <- fread("data/Missense_variants/Liver_missense_dosage_wo_X_MHC.raw")

dosage_trans_secr <-read_table("data/Missense_variants/Liver_missense_transmemb_secr_dosage.raw")

dosage_only_transmemb <- read_table("data/Missense_variants/Liver_missense_transmembrane_dosage.raw")

dosage_liver <- read_table("data/Missense_variants/Liver_missense_liver_spesific_dosage.raw")

### Import the '23_11_13_Clinical_data_pairs' file 

Clin_data_pairs <- read_csv("data/23_11_13_Clinical_data_pairs")
# n = 666

# Check how many patients had acute rejections
sum(Clin_data_pairs$AR_status == 1)
# acute rejection cases: n = 277

# Create a phenotype file from the 'Clin_data_pairs' file:
Covariates <- select(Clin_data_pairs,
                     Pair_number, R_pseudo, D_pseudo, 
                     R_sex, R_age, D_age, D_sex, 
                     Cold_ischemia_time_minutes, AR_status,
                     Eplets_total_HLAI, Eplets_total_HLAII,
                     AR_Cox_time, Transplantation_year,
                     Autoimmune_status, CNI_type_initial_status,
                     LR_Cox_time, Late_rejection_status,
                     Graft_loss_status, graft_loss_months,
                     Death_status, Death_Cox_time_months)

write.table(Covariates,
            file = "./data/Covariates",
            row.names = F, col.names = T, quote = F)

###############################################################################
### Inner join imported dosage files with 'Covariate' file and remove
### unnecessary 

# All 
R_paired_dos_all <- inner_join(Covariates, dosage_all,
                                  by = c("R_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

D_paired_dos_all <- inner_join(Covariates, dosage_all,
                                  by = c("D_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

# Transmembrane and secretory
R_paired_dos_trans_secr <- inner_join(Covariates, dosage_trans_secr,
                                     by = c("R_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

D_paired_dos_trans_secr <- inner_join(Covariates, dosage_trans_secr,
                                     by = c("D_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

# Transmembrane only
R_paired_dos_trans <- inner_join(Covariates,
                                 dosage_only_transmemb,
                                 by = c("R_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

D_paired_dos_trans <- inner_join(Covariates,
                                 dosage_only_transmemb,
                                 by = c("D_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

# Liver specific
R_paired_dos_liver <- inner_join(Covariates,
                                 dosage_liver,
                                 by = c("R_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)


D_paired_dos_liver <- inner_join(Covariates,
                                 dosage_liver,
                                 by = c("D_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -R_pseudo, -D_pseudo, 
         -R_sex, -R_age, -D_age, -D_sex, 
         -Cold_ischemia_time_minutes, -AR_status,
         -Eplets_total_HLAI, -Eplets_total_HLAII,
         -AR_Cox_time, -Transplantation_year,
         -Autoimmune_status, -CNI_type_initial_status,
         -LR_Cox_time, -Late_rejection_status,
         -Graft_loss_status, -graft_loss_months,
         -Death_status, -Death_Cox_time_months)

###############################################################################
### Check that the pair order is the same in both recipient and donor data 
### frames
identical(R_paired_dos_all$Pair_number, 
          D_paired_dos_all$Pair_number)
# TRUE

identical(R_paired_dos_trans_secr$Pair_number, 
          D_paired_dos_trans_secr$Pair_number)
# TRUE

identical(R_paired_dos_trans$Pair_number, 
          D_paired_dos_trans$Pair_number)
# TRUE

identical(R_paired_dos_liver$Pair_number, 
          D_paired_dos_liver$Pair_number)
# TRUE

# Identify recipient-donor pairs and select pair number column

# Transmembrane and secreted
R_D_pairs_trans_secr <- inner_join(R_paired_dos_trans_secr, D_paired_dos_trans_secr,
                                   by = "Pair_number") %>% select(Pair_number)

# Select paired data
R_pairs_dos_trans_secr <- inner_join(R_D_pairs_trans_secr, R_paired_dos_trans_secr, 
                                         by = "Pair_number")


D_pairs_dos_trans_secr <- inner_join(R_D_pairs_trans_secr, D_paired_dos_trans_secr, 
                                         by = "Pair_number")

# Transmembrane only
R_D_pairs_trans <- inner_join(R_paired_dos_trans, D_paired_dos_trans,
                                   by = "Pair_number") %>% select(Pair_number)

# Select paired data
R_pairs_dos_trans <- inner_join(R_D_pairs_trans, R_paired_dos_trans, 
                                     by = "Pair_number")


D_pairs_dos_trans <- inner_join(R_D_pairs_trans, D_paired_dos_trans, 
                                     by = "Pair_number")

# Liver related
R_D_pairs_liver <- inner_join(R_paired_dos_liver, D_paired_dos_liver,
                              by = "Pair_number") %>% select(Pair_number)

# Select paired data
R_pairs_dos_liver <- inner_join(R_D_pairs_liver, R_paired_dos_liver, 
                                by = "Pair_number")


D_pairs_dos_liver <- inner_join(R_D_pairs_liver, D_paired_dos_liver, 
                                by = "Pair_number")

###############################################################################
### Calculate the mismatch sums in each protein group

# All
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

# Transmembrane and secreted
Mm_trans_secr_res <- sapply(2:ncol(R_pairs_dos_trans_secr), 
                               function(i) {Mismatch(R_pairs_dos_trans_secr[,i], 
                                                     D_pairs_dos_trans_secr[,i])})

Mm_trans_secr_res_df <- data.frame(Pair_number=R_pairs_dos_trans_secr$Pair_number, 
                                   Mm_trans_secr=rowSums(Mm_trans_secr_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_trans_secr <- inner_join(R_covariates_mm_all, 
                                         Mm_trans_secr_res_df,
                                         by = "Pair_number")

# Transmembrane only
Mm_transmemb_res <- sapply(2:ncol(R_pairs_dos_trans), 
                                  function(i) {Mismatch(R_pairs_dos_trans[,i], 
                                                        D_pairs_dos_trans[,i])})

Mm_transmemb_res_df <- data.frame(Pair_number=R_pairs_dos_trans$Pair_number,
                                         Mm_transmemb=rowSums(Mm_transmemb_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_transmemb <- inner_join(R_covariates_mm_trans_secr,
                                            Mm_transmemb_res_df,
                                            by = "Pair_number")

# Liver related
Mm_liver_res <- sapply(2:ncol(R_pairs_dos_liver), 
                          function(i) {Mismatch(R_pairs_dos_liver[,i], 
                                                D_pairs_dos_liver[,i])})

Mm_liver_res_df <- data.frame(Pair_number=R_pairs_dos_liver$Pair_number,
                                 Mm_liver=rowSums(Mm_liver_res))

# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_liver <- inner_join(R_covariates_mm_transmemb,
                                    Mm_liver_res_df,
                                    by = "Pair_number")

# Write out the final data file
write.table(R_covariates_mm_liver,
            file = "./data/Missense_variants/R_covariates_mm_liver.txt",
            row.names = F, col.names = T, quote = F)
