###############################################################################
### Perform deletion variant sum analyses using adjusted Cox proportional 
### hazards model for overall survival (OS) endpoint
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
## Import matched covariate genomic collision file
R_dos_pheno_dels_collision <- read_table("data/Deletion_variants/R_dos_pheno_dels_collision.txt")
D_dos_pheno_dels_collision <- read_table("data/Deletion_variants/D_dos_pheno_dels_collision.txt")

##############################################################################
### Variant sum analyses
### Variants: "rs11985201", "rs2342606", "rs2174926", "rs1944862"

### OS ENDPOINT VARIANT SUMS
### The survival analysis: the association of variant sum mismatch to death
### event

# Variant_sum
OS_cox_VARSUM <- coxph(Surv(Death_Cox_time_months, 
                            Death_status) ~ Variant_sum +
                         R_age + D_age + R_sex + 
                         D_sex + Cold_ischemia_time_minutes + 
                         Eplets_total_HLAI + Eplets_total_HLAII +
                         Transplantation_year +
                         Autoimmune_status +
                         CNI_type_initial_status + 
                         AR_status,
                       data = R_dos_pheno_dels_collision) 
summary(OS_cox_VARSUM)
OS_cox_adj_sum <- summary(OS_cox_VARSUM) %>% coef()
OS_cox_adj_CI <- as.data.frame(exp(confint(OS_cox_VARSUM)))
OS_COEF_CI_adj <- bind_cols(OS_cox_adj_sum, OS_cox_adj_CI)
OS_COEF_CI_adj <- round(OS_COEF_CI_adj, digits = 3)
OS_COEF_CI_adj$`HR(95%_CI)` <- paste0(OS_COEF_CI_adj$`exp(coef)`, "(",
                                      OS_COEF_CI_adj$`2.5 %`, "-", 
                                      OS_COEF_CI_adj$`97.5 %`,
                                      ")")
OS_COEF_CI_adj <- rename(OS_COEF_CI_adj, "HR" = `exp(coef)`,
                         "p_value" = `Pr(>|z|)`,
                         `2.5%` = `2.5 %`,
                         `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates")

# AR_variant_sum
OS_cox_VARSUM_AR <- coxph(Surv(Death_Cox_time_months, 
                            Death_status) ~ AR_variant_sum +
                         R_age + D_age + R_sex + 
                         D_sex + Cold_ischemia_time_minutes + 
                         Eplets_total_HLAI + Eplets_total_HLAII +
                         Transplantation_year +
                         Autoimmune_status +
                         CNI_type_initial_status + 
                         AR_status,
                       data = R_dos_pheno_dels_collision) 
summary(OS_cox_VARSUM_AR)

# LR_AR_variant_sum
OS_cox_VARSUM_LR_AR <- coxph(Surv(Death_Cox_time_months, 
                               Death_status) ~ LR_AR_variant_sum +
                            R_age + D_age + R_sex + 
                            D_sex + Cold_ischemia_time_minutes + 
                            Eplets_total_HLAI + Eplets_total_HLAII +
                            Transplantation_year +
                            Autoimmune_status +
                            CNI_type_initial_status + 
                            AR_status,
                          data = R_dos_pheno_dels_collision) 
summary(OS_cox_VARSUM_LR_AR)