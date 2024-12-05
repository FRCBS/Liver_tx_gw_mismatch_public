###############################################################################
### Calculate the variant sums and perform deletion variant sum analyses using
### adjusted Cox proportional hazards model for acute rejection (AR) endpoint
###############################################################################
### General information

# Prerequisites:
# 1. Run script 02_Deletion_variant_mismatches.R
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
##############################################################################
### Variants: "rs11985201", "rs2342606", "rs2174926", "rs1944862"

### Calculate the variant sums
R_dos_pheno_dels_collision$Variant_sum <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col", "rs2174926_col", "rs1944862_col")]))
D_dos_pheno_dels_collision$Variant_sum <- c(Del_sum_collision = rowSums(D_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col", "rs2174926_col", "rs1944862_col")]))

R_dos_pheno_dels_collision$AR_variant_sum <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col")]))
D_dos_pheno_dels_collision$AR_variant_sum <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col")]))

R_dos_pheno_dels_collision$LR_AR_variant_sum <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col","rs1944862_col")]))
D_dos_pheno_dels_collision$LR_AR_variant_sum <- c(Del_sum_collision = rowSums(D_dos_pheno_dels_collision[c("rs11985201_col", "rs2342606_col","rs1944862_col")]))

R_dos_pheno_dels_collision$LR_AR_combined <- R_dos_pheno_dels_collision$LR_AR_variant_sum

## Combine recipients with mismatch score 2 and 3 in one group
R_dos_pheno_dels_collision$LR_AR_combined <- gsub("2", "3", R_dos_pheno_dels_collision$LR_AR_combined)
R_dos_pheno_dels_collision$LR_AR_combined <- gsub("3", "2&3", R_dos_pheno_dels_collision$LR_AR_combined)

## Write out the files
write.table(R_dos_pheno_dels_collision,
            file = "data/Deletion_variants/R_dos_pheno_dels_collision.txt",
            row.names = F, col.names = T, quote = F)

write.table(D_dos_pheno_dels_collision,
            file = "data/Deletion_variants/D_dos_pheno_dels_collision.txt",
            row.names = F, col.names = T, quote = F)
#######################################
### Cox analyses
#######################################
## AR ENDPOINT VARIANT SUMS
## The survival analyses: the association of variant sum mismatch to acute
## rejection event

# Variant_sum
AR_cox_VARSUM <- coxph(Surv(AR_Cox_time, AR_status) ~ Variant_sum +
                                R_age + D_age + R_sex + 
                                D_sex + Cold_ischemia_time_minutes + 
                                Eplets_total_HLAI + Eplets_total_HLAII +
                                Transplantation_year +
                                Autoimmune_status +
                                CNI_type_initial_status,
                              data = R_dos_pheno_dels_collision) 

summary(AR_cox_VARSUM_cox)

# AR_variant_sum
AR_cox_VARSUM_AR <- coxph(Surv(AR_Cox_time, AR_status) ~ AR_variant_sum +
                             R_age + D_age + R_sex + 
                             D_sex + Cold_ischemia_time_minutes + 
                             Eplets_total_HLAI + Eplets_total_HLAII +
                             Transplantation_year +
                             Autoimmune_status +
                             CNI_type_initial_status,
                           data = R_dos_pheno_dels_collision) 

summary(AR_cox_VARSUM_AR)

# LR_AR_variant_sum
AR_cox_VARSUM_LR_AR <- coxph(Surv(AR_Cox_time, AR_status) ~ LR_AR_variant_sum +
                             R_age + D_age + R_sex + 
                             D_sex + Cold_ischemia_time_minutes + 
                             Eplets_total_HLAI + Eplets_total_HLAII +
                             Transplantation_year +
                             Autoimmune_status +
                             CNI_type_initial_status,
                           data = R_dos_pheno_dels_collision) 

summary(AR_cox_VARSUM_LR_AR)

AR_cox_adj_sum <- summary(AR_cox_VARSUM_LR_AR) %>% coef()
AR_cox_adj_CI <- as.data.frame(exp(confint(AR_cox_VARSUM_LR_AR)))
AR_COEF_CI_adj <- bind_cols(AR_cox_adj_sum, AR_cox_adj_CI)
AR_COEF_CI_adj <- round(AR_COEF_CI_adj, digits = 3)
AR_COEF_CI_adj$`HR(95%_CI)` <- paste0(AR_COEF_CI_adj$`exp(coef)`, "(",
                                      AR_COEF_CI_adj$`2.5 %`, "-", 
                                      AR_COEF_CI_adj$`97.5 %`,
                                      ")")
AR_COEF_CI_adj <- rename(AR_COEF_CI_adj, "HR" = `exp(coef)`,
                         "p_value" = `Pr(>|z|)`,
                         `2.5%` = `2.5 %`,
                         `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates")

## Kapplan-Meier plot with risk table for the manuscript
AR_LR_varsum_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ LR_AR_combined, 
                data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Acute rejection -free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE, fontsize = 5,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Mismatch score = 0", "Mismatch score = 1",
                  "Mismatch score = 2 or 3"),
  font.legend = c(15),
  legend = c(0.6,0.9))
AR_LR_varsum_KM

# The same plot with additional details written in the picture (HR and 95% CI)
AR_LR_varsum_KM$plot <- AR_LR_varsum_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.377, 95% CI 1.139-1.664 \n P-value 0.001",
                    size = 6)
# Customize risk table title using ggplot verbs
AR_LR_varsum_KM$table <- AR_LR_varsum_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))

AR_LR_varsum_KM

jpeg('./results/Deletion_variants/AR_LR_combined_variant_sum_analysis_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_LR_varsum_KM)
dev.off()
