###############################################################################
### Perform missense variant analyses using adjusted Cox proportional 
### hazards model for time to overall survival (OS)
###############################################################################
### General information

# Prerequisites:
# 1) Run script 05_AR_missense_variant_analyses.R
#    data/Missense_variants/R_cov_mm_liver.txt
# 2) results/Missense_variants

###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

###############################################################################
## Import covariate file including all required covariates and mismatch sums
## for all four protein groups 
R_cov_mm_liver <- fread("data/Missense_variants/R_cov_mm_liver.txt")

str(R_cov_mm_liver)
# Classes ‘data.table’ and 'data.frame':	666 obs. of  29 variables:
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
# $ Mm_all                    : int  4782 5091 4860 4518 4647 4818 4851 5051 4909 4923 ...
# $ Mm_transm_secr             : int  1636 1709 1627 1460 1583 1666 1629 1663 1657 1632 ...
# $ Mm_transm              : int  1202 1238 1211 1058 1208 1233 1246 1226 1247 1212 ...
# $ Mm_liver                  : int  578 559 494 485 491 497 529 505 517 538 ...
# $ quartile                  : int  2 4 3 1 1 3 3 4 4 4 ...
# $ quart_t_secr              : int  3 4 2 1 1 3 2 3 3 3 ...
# $ quart_transmemb           : int  2 3 2 1 2 3 3 3 3 2 ...
# $ quart_liver               : int  4 4 1 1 1 1 3 2 2 3 ...
#- attr(*, ".internal.selfref")=<externalptr> 

################################################################################
### The survival analysis: the mismatch sum association to time to acute 
### rejection event

## All proteins

# Cox proportional hazards model for adjusted data
OS_cox_all <- coxph(Surv(Death_Cox_time_months, Death_status) ~ Mm_all + 
                      R_sex + D_sex + R_age + D_age + 
                      Cold_ischemia_time_minutes +
                      Eplets_total_HLAI + Eplets_total_HLAII + 
                      Transplantation_year + Autoimmune_status + 
                      CNI_type_initial_status + AR_status,
                    data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
OS_cox_all_sum <- summary(OS_cox_all) %>% coef()
OS_cox_all_COEF <- bind_cols(OS_cox_all_sum, as.data.frame(exp(confint(OS_cox_all))))
OS_cox_all_COEF <- round(OS_cox_all_COEF, digits = 3)
OS_cox_all_COEF$`HR(95%_CI)`<- paste0(OS_cox_all_COEF$`exp(coef)`, "(",
                                      OS_cox_all_COEF$`2.5 %`, "-",  
                                      OS_cox_all_COEF$`97.5 %`, ")" )

OS_cox_all_COEF <- rename(OS_cox_all_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5%`= `2.5 %`,
                          `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_cox_all_COEF,
            "./results/Missense_variants/OS_Cox_all_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## All proteins

## Analyze the quartiles of the mismatch data in all missense variants

# Cox proportional hazards model for adjusted data
OS_Cox_all_quartiles <- coxph(Surv(Death_Cox_time_months, 
                                   Death_status) ~ quartile + R_sex +
                                D_sex + R_age + D_age + 
                                Cold_ischemia_time_minutes + 
                                Eplets_total_HLAI + 
                                Eplets_total_HLAII + Transplantation_year + 
                                Autoimmune_status + 
                                CNI_type_initial_status + AR_status,
                              data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
OS_Cox_all_quart_sum <- summary(OS_Cox_all_quartiles) %>% coef()
OS_Cox_all_quart_COEF <- bind_cols(OS_Cox_all_quart_sum,
                                   as.data.frame(exp(confint(OS_Cox_all_quartiles))))
OS_Cox_all_quart_COEF <- round(OS_Cox_all_quart_COEF, digits = 3)
OS_Cox_all_quart_COEF$`HR(95%_CI)`<- paste0(OS_Cox_all_quart_COEF$`exp(coef)`, 
                                            "(",
                                            OS_Cox_all_quart_COEF$`2.5 %`, "-",  
                                            OS_Cox_all_quart_COEF$`97.5 %`, ")")

OS_Cox_all_quart_COEF <- rename(OS_Cox_all_quart_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_Cox_all_quart_COEF,
            "./results/Missense_variants/OS_Mm_missense_adjusted_all_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
OS_all_quart_KM <- ggsurvplot(
  fit = survfit(Surv(Death_Cox_time_months, Death_status) ~ quartile, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE, fontsize = 5,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 4461-4721", "Q2: MM sum 4722-4798",
                  "Q3: MM sum 4799-4872", "Q4: MM sum 4873-5604"),
  font.legend = c(16),
  legend = c(0.20,0.35))
OS_all_quart_KM

# Add HR and 95% CI
OS_all_quart_KM$plot <- OS_all_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.907, 95% CI 0.784-1.049 \n P-value 0.189",
                    size = 6)
# Customize risk table title using ggplot verbs
OS_all_quart_KM$table <- OS_all_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
OS_all_quart_KM

jpeg('./results/Missense_variants/OS_Missense_variant_analyses_all_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(OS_all_quart_KM)
dev.off()

###############################################################################
## Transmembrane and secreted

# Cox proportional hazards model for adjusted data
OS_cox_trans_secr <- coxph(Surv(Death_Cox_time_months, 
                                Death_status) ~ Mm_transm_secr + 
                             R_sex + D_sex + R_age + D_age + 
                             Cold_ischemia_time_minutes +
                             Eplets_total_HLAI + Eplets_total_HLAII + 
                             Transplantation_year + Autoimmune_status + 
                             CNI_type_initial_status + AR_status,
                           data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
OS_cox_trans_secr_sum <- summary(OS_cox_trans_secr) %>% coef()
OS_cox_trans_secr_COEF <- bind_cols(OS_cox_trans_secr_sum, 
                                    as.data.frame(exp(confint(OS_cox_trans_secr))))
OS_cox_trans_secr_COEF <- round(OS_cox_trans_secr_COEF, digits = 3)
OS_cox_trans_secr_COEF$`HR(95%_CI)`<- paste0(OS_cox_trans_secr_COEF$`exp(coef)`,
                                             "(", OS_cox_trans_secr_COEF$`2.5 %`,
                                             "-", OS_cox_trans_secr_COEF$`97.5 %`,
                                             ")")
OS_cox_trans_secr_COEF <- rename(OS_cox_trans_secr_COEF, `p_value` = `Pr(>|z|)`,
                                 `HR` = `exp(coef)`,
                                 `2.5%`= `2.5 %`,
                                 `97.5%` = `97.5 %`) %>%
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_cox_trans_secr_COEF,
            "./results/Missense_variants/OS_Cox_trans_secr_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Transmembrane and secreted

## Analyze the quartiles of the mismatch data in transmembrane and secreted
## missense variants

# Cox proportional hazards model for adjusted data
OS_Cox_t_secr_quart <- coxph(Surv(Death_Cox_time_months, 
                                  Death_status) ~ quart_t_secr + R_sex +
                               D_sex + R_age + D_age + Cold_ischemia_time_minutes
                             + Eplets_total_HLAI + 
                               Eplets_total_HLAII + Transplantation_year + 
                               Autoimmune_status + 
                               CNI_type_initial_status + AR_status,
                             data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
OS_Cox_t_secr_quart_sum <- summary(OS_Cox_t_secr_quart) %>% coef()
OS_Cox_t_secr_quart_COEF <- bind_cols(OS_Cox_t_secr_quart_sum,
                                      as.data.frame(exp(confint(OS_Cox_t_secr_quart))))
OS_Cox_t_secr_quart_COEF <- round(OS_Cox_t_secr_quart_COEF, digits = 3)
OS_Cox_t_secr_quart_COEF$`HR(95%_CI)`<- paste0(OS_Cox_t_secr_quart_COEF$`exp(coef)`,
                                               "(",
                                               OS_Cox_t_secr_quart_COEF$`2.5 %`,
                                               "-",
                                               OS_Cox_t_secr_quart_COEF$`97.5 %`,
                                               ")")

OS_Cox_t_secr_quart_COEF <- rename(OS_Cox_t_secr_quart_COEF, `p_value` = `Pr(>|z|)`,
                                   `HR` = `exp(coef)`,
                                   `2.5%`= `2.5 %`,
                                   `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_Cox_t_secr_quart_COEF,
            "./results/Missense_variants/OS_Cox_trans_secr_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
OS_t_secr_quart_KM <- ggsurvplot(
  fit = survfit(Surv(Death_Cox_time_months, Death_status) ~ quart_t_secr, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE, fontsize = 5,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 1460-1593", "Q2: MM sum 1594-1629",
                  "Q3: MM sum 1630-1675", "Q4: MM sum 1676-1952"),
  font.legend = c(16),
  legend = c(0.20,0.35))
OS_t_secr_quart_KM

# Add HR and 95% CI
OS_t_secr_quart_KM$plot <- OS_t_secr_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.897, 95% CI 0.777-1.036 \n P-value 0.138",
                    size = 6)
# Customize risk table title using ggplot verbs
OS_t_secr_quart_KM$table <- OS_t_secr_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
OS_t_secr_quart_KM

jpeg('./results/Missense_variants/OS_Missense_variant_analyses_transmembrane_secr_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(OS_t_secr_quart_KM)
dev.off()

###############################################################################
## Transmembrane only

# Cox proportional hazards model for adjusted data
OS_cox_transmemb <- coxph(Surv(Death_Cox_time_months, 
                               Death_status) ~ Mm_transm + 
                            R_sex + D_sex + R_age + D_age + 
                            Cold_ischemia_time_minutes +
                            Eplets_total_HLAI + Eplets_total_HLAII + 
                            Transplantation_year + Autoimmune_status + 
                            CNI_type_initial_status + AR_status,
                          data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
OS_cox_transmemb_sum <- summary(OS_cox_transmemb) %>% coef()
OS_cox_transmemb_COEF <- bind_cols(OS_cox_transmemb_sum, 
                                   as.data.frame(exp(confint(OS_cox_transmemb))))
OS_cox_transmemb_COEF <- round(OS_cox_transmemb_COEF, digits = 3)
OS_cox_transmemb_COEF$`HR(95%_CI)`<- paste0(OS_cox_transmemb_COEF$`exp(coef)`,
                                            "(", OS_cox_transmemb_COEF$`2.5 %`,
                                            "-", OS_cox_transmemb_COEF$`97.5 %`, 
                                            ")")

OS_cox_transmemb_COEF <- rename(OS_cox_transmemb_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_cox_transmemb_COEF,
            "./results/Missense_variants/OS_Cox_transmemb_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Transmembrane only

## Analyze the quartiles of the mismatch data in transmembrane missense 
## variants

# Cox proportional hazards model for adjusted data
OS_cox_trans_quart <- coxph(Surv(Death_Cox_time_months, 
                                 Death_status) ~ quart_transmemb + R_sex +
                              D_sex + R_age + D_age + Cold_ischemia_time_minutes
                            + Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status + AR_status,
                            data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
OS_Cox_trans_quart_sum <- summary(OS_cox_trans_quart) %>% coef()
OS_Cox_trans_quart_COEF <- bind_cols(OS_Cox_trans_quart_sum,
                                     as.data.frame(exp(confint(OS_cox_trans_quart))))
OS_Cox_trans_quart_COEF <- round(OS_Cox_trans_quart_COEF, digits = 3)
OS_Cox_trans_quart_COEF$`HR(95%_CI)`<- paste0(OS_Cox_trans_quart_COEF$`exp(coef)`,
                                              "(", OS_Cox_trans_quart_COEF$`2.5 %`,
                                              "-", OS_Cox_trans_quart_COEF$`97.5 %`, 
                                              ")")

OS_Cox_trans_quart_COEF <- rename(OS_Cox_trans_quart_COEF, 
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_Cox_trans_quart_COEF,
            "./results/Missense_variants/OS_Cox_transmemb_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
OS_trans_quart_KM <- ggsurvplot(
  fit = survfit(Surv(Death_Cox_time_months, Death_status) ~ quart_transmemb, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE, fontsize = 5,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 1058-1177", "Q2: MM sum 1178-1214",
                  "Q3: MM sum 1215-1250", "Q4: MM sum 1251-1488"),
  font.legend = c(16),
  legend = c(0.20,0.35))
OS_trans_quart_KM

# Add HR and 95% CI
OS_trans_quart_KM$plot <- OS_trans_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.934, 95% CI 0.809-1.077 \n P-value 0.348",
                    size = 6)
# Customize risk table title using ggplot verbs
OS_trans_quart_KM$table <- OS_trans_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
OS_trans_quart_KM

jpeg('./results/Missense_variants/OS_Missense_variant_analyses_transmembrane_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(OS_trans_quart_KM)
dev.off()

###############################################################################
## Liver-related

# Cox proportional hazards model for adjusted data
OS_cox_liver <- coxph(Surv(Death_Cox_time_months, Death_status) ~ Mm_liver + 
                        R_sex + D_sex + R_age + D_age + 
                        Cold_ischemia_time_minutes +
                        Eplets_total_HLAI + Eplets_total_HLAII + 
                        Transplantation_year + Autoimmune_status + 
                        CNI_type_initial_status + AR_status,
                      data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
OS_cox_liver_sum <- summary(OS_cox_liver) %>% coef()
OS_cox_liver_COEF <- bind_cols(OS_cox_liver_sum, 
                               as.data.frame(exp(confint(OS_cox_liver))))
OS_cox_liver_COEF <- round(OS_cox_liver_COEF, digits = 3)
OS_cox_liver_COEF$`HR(95%_CI)`<- paste0(OS_cox_liver_COEF$`exp(coef)`,
                                        "(", OS_cox_liver_COEF$`2.5 %`,
                                        "-", OS_cox_liver_COEF$`97.5 %`, 
                                        ")")

OS_cox_liver_COEF <- rename(OS_cox_liver_COEF, `p_value` = `Pr(>|z|)`,
                            `HR` = `exp(coef)`,
                            `2.5%`= `2.5 %`,
                            `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_cox_liver_COEF,
            "./results/Missense_variants/OS_Cox_liver_related_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Liver-related

## Analyze the quartiles of the mismatch data in liver-related missense 
## variants

# Cox proportional hazards model for adjusted data
OS_cox_liver_quart <- coxph(Surv(Death_Cox_time_months, 
                                 Death_status) ~ quart_liver + R_sex +
                              D_sex + R_age + D_age + Cold_ischemia_time_minutes
                            + Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status + AR_status,
                            data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
OS_Cox_liver_quart_sum <- summary(OS_cox_liver_quart) %>% coef()
OS_Cox_liver_quart_COEF <- bind_cols(OS_Cox_liver_quart_sum,
                                     as.data.frame(exp(confint(OS_cox_liver_quart))))
OS_Cox_liver_quart_COEF <- round(OS_Cox_liver_quart_COEF, digits = 3)
OS_Cox_liver_quart_COEF$`HR(95%_CI)`<- paste0(OS_Cox_liver_quart_COEF$`exp(coef)`,
                                              "(", OS_Cox_liver_quart_COEF$`2.5 %`,
                                              "-", OS_Cox_liver_quart_COEF$`97.5 %`, 
                                              ")")

OS_Cox_liver_quart_COEF <- rename(OS_Cox_liver_quart_COEF,
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(OS_Cox_liver_quart_COEF,
            "./results/Missense_variants/OS_Cox_liver_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
OS_liver_quart_KM <- ggsurvplot(
  fit = survfit(Surv(Death_Cox_time_months, Death_status) ~ quart_liver, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 438-498", "Q2: MM sum 499-518",
                  "Q3: MM sum 519-539", "Q4: MM sum 540-629"),
  font.legend = c(16),
  legend = c(0.20,0.35))
OS_liver_quart_KM

# Add HR and 95% CI
OS_liver_quart_KM$plot <- OS_liver_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.956, 95% CI 0.828-1.105 \n P-value 0.546",
                    size = 6)
# Customize risk table title using ggplot verbs
OS_liver_quart_KM$table <- OS_liver_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
OS_liver_quart_KM

jpeg('./results/Missense_variants/OS_Missense_variant_analyses_liver_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(OS_liver_quart_KM)
dev.off()

################################################################################
### Perform multiple test corrections for OS for mismatch sum and quartiles

## Mismatch sums

# Create combined results
OS_results <- list(OS_cox_all_COEF, OS_cox_trans_secr_COEF, 
                   OS_cox_transmemb_COEF, OS_cox_liver_COEF)
names(OS_results) <- c("OS_cox_all_COEF", "OS_cox_trans_secr_COEF", 
                       "OS_cox_transmemb_COEF", "OS_cox_liver_COEF") 

# Extract mismatch covariates and p values
OS_results <- map((OS_results),
                  function(x) {
                    DATA <- x[1,] 
                    DATA <- select(DATA, covariates, p_value)
                    return(DATA)
                  })

OS_results <- bind_rows(OS_results)

# Perform multiple test correction
OS_results$Bonferroni <- p.adjust(OS_results[,2], method = "bonferroni", n = 4)
OS_results$Holm <- p.adjust(OS_results[,2], method = "holm", n = 4)
OS_results$FDR <- p.adjust(OS_results[,2], method = "fdr", n = 4)

write.table(OS_results, 
            "results/Missense_variants/OS_Cox_Mm_mismatch_multiple_testing_results",
            quote = F, row.names = F)

#######################################
## Mismatch quartiles

# Create combined results
OS_quart_results <- list(OS_Cox_all_quart_COEF, OS_Cox_t_secr_quart_COEF, 
                         OS_Cox_trans_quart_COEF, OS_Cox_liver_quart_COEF)
names(OS_quart_results) <- c("OS_Cox_all_quart_COEF", "OS_Cox_t_secr_quart_COEF", 
                             "OS_Cox_trans_quart_COEF", "OS_Cox_liver_quart_COEF") 

# Extract mismatch covariates and p values
OS_quart_results <- map((OS_quart_results),
                        function(x) {
                          DATA <- x[1,] 
                          DATA <- select(DATA, covariates, p_value)
                          return(DATA)
                        })

OS_quart_results <- bind_rows(OS_quart_results)
OS_quart_results$Bonferroni <- p.adjust(OS_quart_results[,2], 
                                        method = "bonferroni", n = 4)
OS_quart_results$Holm <- p.adjust(OS_quart_results[,2], method = "holm", n = 4)
OS_quart_results$FDR <- p.adjust(OS_quart_results[,2], method = "fdr", n = 4)

write.table(OS_quart_results, 
            "results/Missense_variants/OS_Cox_Mm_quartiles_multiple_testing_results",
            quote = F, row.names = F)
###############################################################################
