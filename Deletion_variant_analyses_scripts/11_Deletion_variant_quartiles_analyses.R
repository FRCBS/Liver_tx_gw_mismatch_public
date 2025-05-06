###############################################################################
### Calculate the sum of deletions/mismatches (40 SNPs) per individual and
### analyze the association of deletion sum quartiles on acute rejection (AR)
### endpoint.
### acute rejection (AR), graft loss (GL), overall survival (OS) and 
### late rejection (LR)
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
# The value for mismatch is 1, other cases are 0
# recipient 0 + donor 0 = 0
# recipient 0 + donor 1 = 0
# recipient 1 + donor 0 = 1
# recipient 1 + donor 1 = 0

## Import matched covariate genomic collision files
R_dos_pheno_dels_collision <- read_table("data/Deletion_variants/R_dos_pheno_dels_collision.txt")
D_dos_pheno_dels_collision <- read_table("data/Deletion_variants/D_dos_pheno_dels_collision.txt")

## Calculate the deletion sum
R_dos_pheno_dels_collision$Del_sum_collision <- c(Del_sum_collision = rowSums(R_dos_pheno_dels_collision[103:142]))
D_dos_pheno_dels_collision$Del_sum_collision <- c(Del_sum_collision = rowSums(D_dos_pheno_dels_collision[103:142]))

## Create a new column with deletion sum quartiles
quantile(R_dos_pheno_dels_collision$Del_sum_collision)
  # 0%  25%  50%  75% 100% 
  # 0    2    3    4    9 

R_dos_pheno_dels_collision$Del_sum_quartiles <- as.integer(cut(R_dos_pheno_dels_collision$Del_sum_collision, 
                                                               quantile(R_dos_pheno_dels_collision$Del_sum_collision, 
                                                                        probs = 0:4/4), 
                                                               include.lowest = TRUE))
                                                               
###############################################################################
### AR DELETION SUM QUARTILES
###############################################################################
### Cox proportional hazards model for adjusted data
AR_delsum_quart <- coxph(Surv(AR_Cox_time, AR_status) ~ Del_sum_quartiles +
                       R_age + D_age + R_sex + 
                       D_sex + Cold_ischemia_time_minutes + 
                       Eplets_total_HLAI + Eplets_total_HLAII +
                       Transplantation_year +
                       Autoimmune_status +
                       CNI_type_initial_status,
                     data = R_dos_pheno_dels_collision) 

AR_del_quart_sum <- summary(AR_delsum_quart) %>% coef()
AR_del_quart_CI <- as.data.frame(exp(confint(AR_delsum_quart)))
AR_COEF_del_quart <- bind_cols(AR_del_quart_sum, AR_del_quart_CI)
AR_COEF_del_quart <- round(AR_COEF_del_quart, digits = 3)
AR_COEF_del_quart$`HR(95%_CI)` <- paste0(AR_COEF_del_quart$`exp(coef)`, "(",
                                      AR_COEF_del_quart$`2.5 %`, "-", 
                                      AR_COEF_del_quart$`97.5 %`,
                                      ")")
AR_COEF_del_quart <- rename(AR_COEF_del_quart, "HR" = `exp(coef)`,
                         "p_value" = `Pr(>|z|)`,
                         `2.5%` = `2.5 %`,
                         `97.5%` = `97.5 %`)  %>% 
  rownames_to_column("covariates")

write.table(AR_COEF_del_quart, 
            "./results/Deletion_variants/AR_Cox_adjusted_deletions_sum_mismatch_quartiles",
            sep = "\t", quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the supplementary
AR_del_quart_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ Del_sum_quartiles,
                data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "", 
  legend.labs = c("Q1: Del sum 0-2", "Q2: Del sum 2-3",
                  "Q3: Del sum 3-4", "Q4: Del sum 4-9"),
  font.legend = c(16),
  legend = c(0.2,0.25))

AR_del_quart_KM$plot <- AR_del_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 170, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.014, 95% CI 0.912-1.127 \n P-value 0.799",
                    size = 6)
AR_del_quart_KM$table <- AR_del_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
AR_del_quart_KM

jpeg('./results/Deletion_variants/AR_Deletion_sum_mismatch_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_del_quart_KM)
dev.off()

##############################################################################
### LR DELETION SUM QUARTILES
### Cox proportional hazards model for adjusted data
LR_delsum_quart <- coxph(Surv(LR_Cox_time, 
                              Late_rejection_status) ~ Del_sum_quartiles +
                           R_age + D_age + R_sex + 
                           D_sex + Cold_ischemia_time_minutes + 
                           Eplets_total_HLAI + Eplets_total_HLAII +
                           Transplantation_year +
                           Autoimmune_status +
                           CNI_type_initial_status,
                         data = R_dos_pheno_dels_collision) 

LR_del_quart_sum <- summary(LR_delsum_quart) %>% coef()
LR_del_quart_CI <- as.data.frame(exp(confint(LR_delsum_quart)))
LR_COEF_del_quart <- bind_cols(LR_del_quart_sum, LR_del_quart_CI)
LR_COEF_del_quart <- round(LR_COEF_del_quart, digits = 3)
LR_COEF_del_quart$`HR(95%_CI)` <- paste0(LR_COEF_del_quart$`exp(coef)`, "(",
                                         LR_COEF_del_quart$`2.5 %`, "-", 
                                         LR_COEF_del_quart$`97.5 %`,
                                         ")")
LR_COEF_del_quart <- rename(LR_COEF_del_quart, "HR" = `exp(coef)`,
                            "p_value" = `Pr(>|z|)`,
                            `2.5%` = `2.5 %`,
                            `97.5%` = `97.5 %`)  %>% 
  rownames_to_column("covariates")

write.table(LR_COEF_del_quart, 
            "./results/Deletion_variants/LR_Cox_adjusted_deletions_sum_mismatch_quartiles",
            sep = "\t", quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the supplementary
LR_del_quart_KM <- ggsurvplot(
  fit = survfit(Surv(LR_Cox_time, 
                     Late_rejection_status) ~ Del_sum_quartiles,
                data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "", 
  legend.labs = c("Q1: Del sum 0-2", "Q2: Del sum 2-3",
                  "Q3: Del sum 3-4", "Q4: Del sum 4-9"),
  font.legend = c(16),
  legend = c(0.2,0.25))

LR_del_quart_KM$plot <- LR_del_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 170, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.934, 95% CI 0.778-1.122 \n P-value 0.467",
                    size = 6)
LR_del_quart_KM$table <- LR_del_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
LR_del_quart_KM

jpeg('./results/Deletion_variants/LR_Deletion_sum_mismatch_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(LR_del_quart_KM)
dev.off()

###############################################################################
### GL DELETION SUM QUARTILES
### Cox proportional hazards model for adjusted data
GL_delsum_quart <- coxph(Surv(graft_loss_months, 
                              Graft_loss_status) ~ Del_sum_quartiles +
                           R_age + D_age + R_sex + 
                           D_sex + Cold_ischemia_time_minutes + 
                           Eplets_total_HLAI + Eplets_total_HLAII +
                           Transplantation_year +
                           Autoimmune_status +
                           CNI_type_initial_status,
                         data = R_dos_pheno_dels_collision) 

GL_del_quart_sum <- summary(GL_delsum_quart) %>% coef()
GL_del_quart_CI <- as.data.frame(exp(confint(GL_delsum_quart)))
GL_COEF_del_quart <- bind_cols(GL_del_quart_sum, GL_del_quart_CI)
GL_COEF_del_quart <- round(GL_COEF_del_quart, digits = 3)
GL_COEF_del_quart$`HR(95%_CI)` <- paste0(GL_COEF_del_quart$`exp(coef)`, "(",
                                         GL_COEF_del_quart$`2.5 %`, "-", 
                                         GL_COEF_del_quart$`97.5 %`,
                                         ")")
GL_COEF_del_quart <- rename(GL_COEF_del_quart, "HR" = `exp(coef)`,
                            "p_value" = `Pr(>|z|)`,
                            `2.5%` = `2.5 %`,
                            `97.5%` = `97.5 %`)  %>% 
  rownames_to_column("covariates")

write.table(GL_COEF_del_quart, 
            "./results/Deletion_variants/GL_Cox_adjusted_deletions_sum_mismatch_quartiles",
            sep = "\t", quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the supplementary
GL_del_quart_KM <- ggsurvplot(
  fit = survfit(Surv(graft_loss_months, 
                     Graft_loss_status) ~ Del_sum_quartiles,
                data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "", 
  legend.labs = c("Q1: Del sum 0-2", "Q2: Del sum 2-3",
                  "Q3: Del sum 3-4", "Q4: Del sum 4-9"),
  font.legend = c(16),
  legend = c(0.2,0.25))

GL_del_quart_KM$plot <- GL_del_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 170, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.831, 95% CI 0.616-1.122 \n P-value 0.227",
                    size = 6)
GL_del_quart_KM$table <- GL_del_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
GL_del_quart_KM

jpeg('./results/Deletion_variants/GL_Deletion_sum_mismatch_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(GL_del_quart_KM)
dev.off()

################################################################################
### OS DELETION SUM QUARTILES
### Cox proportional hazards model for adjusted data
OS_delsum_quart <- coxph(Surv(Death_Cox_time_months, 
                              Death_status) ~ Del_sum_quartiles +
                           R_age + D_age + R_sex + 
                           D_sex + Cold_ischemia_time_minutes + 
                           Eplets_total_HLAI + Eplets_total_HLAII +
                           Transplantation_year +
                           Autoimmune_status +
                           CNI_type_initial_status,
                         data = R_dos_pheno_dels_collision) 

OS_del_quart_sum <- summary(OS_delsum_quart) %>% coef()
OS_del_quart_CI <- as.data.frame(exp(confint(OS_delsum_quart)))
OS_COEF_del_quart <- bind_cols(OS_del_quart_sum, OS_del_quart_CI)
OS_COEF_del_quart <- round(OS_COEF_del_quart, digits = 3)
OS_COEF_del_quart$`HR(95%_CI)` <- paste0(OS_COEF_del_quart$`exp(coef)`, "(",
                                         OS_COEF_del_quart$`2.5 %`, "-", 
                                         OS_COEF_del_quart$`97.5 %`,
                                         ")")
OS_COEF_del_quart <- rename(OS_COEF_del_quart, "HR" = `exp(coef)`,
                            "p_value" = `Pr(>|z|)`,
                            `2.5%` = `2.5 %`,
                            `97.5%` = `97.5 %`)  %>% 
  rownames_to_column("covariates")

write.table(OS_COEF_del_quart, 
            "./results/Deletion_variants/OS_Cox_adjusted_deletions_sum_mismatch_quartiles",
            sep = "\t", quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the supplementary
OS_del_quart_KM <- ggsurvplot(
  fit = survfit(Surv(Death_Cox_time_months, 
                     Death_status) ~ Del_sum_quartiles,
                data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Rejection-free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "", 
  legend.labs = c("Q1: Del sum 0-2", "Q2: Del sum 2-3",
                  "Q3: Del sum 3-4", "Q4: Del sum 4-9"),
  font.legend = c(16),
  legend = c(0.2,0.25))

OS_del_quart_KM$plot <- OS_del_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 170, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.991, 95% CI 0.858-1.144 \n P-value 0.900",
                    size = 6)
OS_del_quart_KM$table <- OS_del_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
OS_del_quart_KM

jpeg('./results/Deletion_variants/OS_Deletion_sum_mismatch_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(OS_del_quart_KM)
dev.off()