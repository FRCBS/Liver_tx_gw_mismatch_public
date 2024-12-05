###############################################################################
### Perform deletion variant sum analyses using adjusted Cox proportional 
### hazards model for late rejection (LR) endpoint
###############################################################################
### General information

# Prerequisites:
# 1. Run script 07_AR_Deletion_variant_sum_analyses.R
#    data/Deletion_variants/R_dos_pheno_dels_collision.txt
#    data/Deletion_variants/D_dos_pheno_dels_collision.txt

################################################################################

library(readr)
library(tidyverse)
library(survival)
library(survminer)
library(dplyr)

################################################################################
## Import matched covariate genomic collision files
R_dos_pheno_dels_collision <- read_table("data/Deletion_variants/R_dos_pheno_dels_collision.txt")
D_dos_pheno_dels_collision <- read_table("data/Deletion_variants/D_dos_pheno_dels_collision.txt")

##############################################################################
### Variant sum analyses
### Variants: "rs11985201", "rs2342606", "rs2174926", "rs1944862"

### LR ENDPOINT VARIANT SUMS
### The survival analysis: the association of variant sum mismatch to late
### rejection event

# Variant_sum
LR_cox_VARSUM <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~ Variant_sum +
                         R_age + D_age + R_sex + 
                         D_sex + Cold_ischemia_time_minutes + 
                         Eplets_total_HLAI + Eplets_total_HLAII +
                         Transplantation_year +
                         Autoimmune_status +
                         CNI_type_initial_status,
                       data = R_dos_pheno_dels_collision) 
summary(LR_cox_VARSUM)
LR_cox_adj_sum <- summary(LR_cox_VARSUM) %>% coef()
LR_cox_adj_CI <- as.data.frame(exp(confint(LR_cox_VARSUM)))
LR_COEF_CI_adj <- bind_cols(LR_cox_adj_sum, LR_cox_adj_CI)
LR_COEF_CI_adj <- round(LR_COEF_CI_adj, digits = 3)
LR_COEF_CI_adj$`HR(95%_CI)` <- paste0(LR_COEF_CI_adj$`exp(coef)`, "(",
                                      LR_COEF_CI_adj$`2.5 %`, "-", 
                                      LR_COEF_CI_adj$`97.5 %`,
                                      ")")
LR_COEF_CI_adj <- rename(LR_COEF_CI_adj, "HR" = `exp(coef)`,
                         "p_value" = `Pr(>|z|)`,
                         `2.5%` = `2.5 %`,
                         `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates")

# AR_variant_sum
LR_cox_VARSUM_AR <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~ AR_variant_sum +
                         R_age + D_age + R_sex + 
                         D_sex + Cold_ischemia_time_minutes + 
                         Eplets_total_HLAI + Eplets_total_HLAII +
                         Transplantation_year +
                         Autoimmune_status +
                         CNI_type_initial_status,
                       data = R_dos_pheno_dels_collision) 
summary(LR_cox_VARSUM_AR)

# LR_AR_variant_sum
LR_cox_VARSUM_AR <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~ AR_variant_sum +
                            R_age + D_age + R_sex + 
                            D_sex + Cold_ischemia_time_minutes + 
                            Eplets_total_HLAI + Eplets_total_HLAII +
                            Transplantation_year +
                            Autoimmune_status +
                            CNI_type_initial_status,
                          data = R_dos_pheno_dels_collision) 
summary(LR_cox_VARSUM_AR)