###############################################################################
### Create phenotype and covariate files for GWAS analysis
###############################################################################
## General information

# Prerequisites:
# 1) data/Covariates
# 2) Run script 13_LD_pruning_and_PCA_analysis.sh
#    data/PCA/Post_PCA.eigenvec
###############################################################################
library(data.table)
library(tidyverse)

###############################################################################
## Import "Covariates" files that includes both covariate and phenotype info, 
## and eigenvec values from PCA analysis

Covariates <- fread("data/Covariates")
eigenvec <- fread("data/PCA/Post_PCA.eigenvec")

## Created phenotype and phenotype files separately for recipients and donors

# Recipients
R_covariates <- select(Covariates, R_pseudo, R_sex, R_age, D_sex, D_age,
                       Cold_ischemia_time_minutes, Eplets_total_HLAI,
                       Eplets_total_HLAII, Transplantation_year,
                       Autoimmune_status, CNI_type_initial_status)

# Donors
D_covariates <- select(Covariates, D_pseudo, R_sex, R_age, D_sex, D_age,
                       Cold_ischemia_time_minutes, Eplets_total_HLAI,
                       Eplets_total_HLAII, Transplantation_year,
                       Autoimmune_status, CNI_type_initial_status)

# Recipients
R_pheno <- select(Covariates, R_pseudo, AR_status,
                  Late_rejection_status, Graft_loss_status, Death_status)

# Donors
D_pheno <- select(Covariates, D_pseudo, AR_status,
                  Late_rejection_status, Graft_loss_status, Death_status)

## Add eigenvec values and modify covariate and pheno files
R_covariates <- inner_join(eigenvec, R_covariates,
                           by = c("IID" = "R_pseudo"))

R_covariates <- rename(R_covariates, "FID" = "#FID") %>% select(FID, IID, PC1,
                                                                PC2, PC3, PC4,
                                                                R_sex, R_age,
                                                                D_sex, D_age,
                                                                Cold_ischemia_time_minutes,
                                                                Eplets_total_HLAI,
                                                                Eplets_total_HLAII,
                                                                Transplantation_year,
                                                                Autoimmune_status,
                                                                CNI_type_initial_status)

D_covariates <- inner_join(eigenvec, D_covariates,
                           by = c("IID" = "D_pseudo"))

D_covariates <- rename(D_covariates, "FID" = "#FID") %>% select(FID, IID, PC1,
                                                                PC2, PC3, PC4,
                                                                R_sex, R_age,
                                                                D_sex, D_age,
                                                                Cold_ischemia_time_minutes,
                                                                Eplets_total_HLAI,
                                                                Eplets_total_HLAII,
                                                                Transplantation_year,
                                                                Autoimmune_status,
                                                                CNI_type_initial_status)

R_pheno$IID <- R_pheno$R_pseudo
R_pheno <- rename(R_pheno, "FID" = "R_pseudo") %>% select(FID, IID, AR_status,
                                                          Late_rejection_status,
                                                          Graft_loss_status,
                                                          Death_status)

D_pheno$IID <- D_pheno$D_pseudo
D_pheno <- rename(D_pheno, "FID" = "D_pseudo") %>% select(FID, IID, AR_status,
                                                          Late_rejection_status,
                                                          Graft_loss_status,
                                                          Death_status)

## Write out phenotype and covariate files
write.table(R_covariates, "data/GWAS/R_covariates.txt", col.names = T,
            row.names = F, quote = F)

write.table(D_covariates, "data/GWAS/D_covariates.txt", col.names = T,
            row.names = F, quote = F)

write.table(R_pheno, "data/GWAS/R_pheno.txt", col.names = T,
            row.names = F, quote = F)

write.table(D_pheno, "data/GWAS/D_pheno.txt", col.names = T,
            row.names = F, quote = F)

