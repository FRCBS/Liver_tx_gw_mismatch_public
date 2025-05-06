###############################################################################
### Cox survival analyses without genetic mismatch component
###############################################################################
### General information

# Prerequisites:
# 1) data/Covariates file
###############################################################################

library(data.table)
library(tidyverse)
library(survival)
library(survminer)

###############################################################################
## Import covariate file
Covariates <- fread("data/Covariates")

###############################################################################
## Perform Cox survival analyses for four endpoints
###############################################################################
# ACUTE REJECTION

# HLAI
AR_cox <- coxph(Surv(AR_Cox_time, AR_status) ~ Eplets_total_HLAI + 
                      Eplets_total_HLAII + 
                      R_sex + D_sex + R_age + D_age + 
                      Cold_ischemia_time_minutes +
                      Transplantation_year + Autoimmune_status + 
                      CNI_type_initial_status,
                    data = Covariates)

AR_cox_sum <- summary(AR_cox) %>% coef()
AR_cox_COEF <- bind_cols(AR_cox_sum, 
                              as.data.frame(exp(confint(AR_cox))))

AR_cox_COEF <- round(AR_cox_COEF, digits = 3) 

AR_cox_COEF$`HR(95%_CI)` <- paste0(AR_cox_COEF$`exp(coef)`, "(",
                                        AR_cox_COEF$`2.5 %`, "-", 
                                        AR_cox_COEF$`97.5 %`,
                                       ")")

AR_cox_COEF <- rename(AR_cox_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5_%`= `2.5 %`,
                          `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_cox_COEF, "results/Cox_wo_mismatch/AR_cox_wo_mismatch",
            quote = F, row.names = F, col.names = T)

#######################################
# LATE REJECTION

# HLAI
LR_cox <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~ Eplets_total_HLAI + 
                       Eplets_total_HLAII + 
                       R_sex + D_sex + R_age + D_age + 
                       Cold_ischemia_time_minutes +
                       Transplantation_year + Autoimmune_status + 
                       CNI_type_initial_status,
                     data = Covariates)

LR_cox_sum <- summary(LR_cox) %>% coef()
LR_cox_COEF <- bind_cols(LR_cox_sum, 
                              as.data.frame(exp(confint(LR_cox))))

LR_cox_COEF <- round(LR_cox_COEF, digits = 3) 

LR_cox_COEF$`HR(95%_CI)` <- paste0(LR_cox_COEF$`exp(coef)`, "(",
                                        LR_cox_COEF$`2.5 %`, "-", 
                                        LR_cox_COEF$`97.5 %`,
                                        ")")

LR_cox_COEF <- rename(LR_cox_COEF, `p_value` = `Pr(>|z|)`,
                           `HR` = `exp(coef)`,
                           `2.5_%`= `2.5 %`,
                           `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_cox_COEF, "results/Cox_wo_mismatch/LR_cox_wo_mismatch",
            quote = F, row.names = F, col.names = T)

#######################################
# GRAFT LOSS

# HLAI
GL_cox <- coxph(Surv(graft_loss_months, Graft_loss_status) ~ Eplets_total_HLAI + 
                       Eplets_total_HLAII + 
                       R_sex + D_sex + R_age + D_age + 
                       Cold_ischemia_time_minutes +
                       Transplantation_year + Autoimmune_status + 
                       CNI_type_initial_status,
                     data = Covariates)

GL_cox_sum <- summary(GL_cox) %>% coef()
GL_cox_COEF <- bind_cols(GL_cox_sum, 
                              as.data.frame(exp(confint(GL_cox))))

GL_cox_COEF <- round(GL_cox_COEF, digits = 3) 

GL_cox_COEF$`HR(95%_CI)` <- paste0(GL_cox_COEF$`exp(coef)`, "(",
                                        GL_cox_COEF$`2.5 %`, "-", 
                                        GL_cox_COEF$`97.5 %`,
                                        ")")

GL_cox_COEF <- rename(GL_cox_COEF, `p_value` = `Pr(>|z|)`,
                           `HR` = `exp(coef)`,
                           `2.5_%`= `2.5 %`,
                           `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_cox_COEF, "results/Cox_wo_mismatch/GL_cox_wo_mismatch",
            quote = F, row.names = F, col.names = T)

#######################################
# OVERALL SURVIVAL

# HLAI
OS_cox <- coxph(Surv(Death_Cox_time_months, Death_status) ~ Eplets_total_HLAI + 
                       Eplets_total_HLAII + 
                       R_sex + D_sex + R_age + D_age + 
                       Cold_ischemia_time_minutes +
                       Transplantation_year + Autoimmune_status + 
                       CNI_type_initial_status,
                     data = Covariates)

OS_cox_sum <- summary(OS_cox) %>% coef()
OS_cox_COEF <- bind_cols(OS_cox_sum, 
                              as.data.frame(exp(confint(OS_cox))))

OS_cox_COEF <- round(OS_cox_COEF, digits = 3) 

OS_cox_COEF$`HR(95%_CI)` <- paste0(OS_cox_COEF$`exp(coef)`, "(",
                                        OS_cox_COEF$`2.5 %`, "-", 
                                        OS_cox_COEF$`97.5 %`,
                                        ")")

OS_cox_COEF <- rename(OS_cox_COEF, `p_value` = `Pr(>|z|)`,
                           `HR` = `exp(coef)`,
                           `2.5_%`= `2.5 %`,
                           `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_cox_COEF, "results/Cox_wo_mismatch/OS_cox_wo_mismatch",
            quote = F, row.names = F, col.names = T)
