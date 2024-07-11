###############################################################################
### Perform missense variant analyses using adjusted Cox proportional 
### hazards model
###############################################################################
### General information

# Prerequisites:
# 1) Run script AR_missense_variant_analyses.R
# 2) data/Missense_variants/R_covariates_mm_liver.txt
# 3) results/Missense_variants

###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

###############################################################################
## Import covariate file including all required covariates and mismatch sums
## for all four protein groups 
R_cov_mm_liver <- fread("data/Missense_variants/R_cov_mm_liver.txt")

# Classes ?data.table? and 'data.frame':	666 obs. of  25 variables:
#$ Pair_number               : int  539 540 542 543 541 544 545 546 547 548 ...
#$ R_pseudo                  : chr  "pseudo1" "pseudo2" "pseudo3" "pseudo4" ...
#$ D_pseudo                  : chr  "pseudo_1" "pseudo_2" "pseudo_3" "pseudo_4" ...
#$ R_sex                     : int  1 1 2 1 2 1 1 1 1 2 ...
#$ R_age                     : int  47 19 49 62 31 52 48 36 30 71 ...
#$ D_age                     : int  55 31 73 57 20 57 50 64 49 64 ...
#$ D_sex                     : int  1 1 2 2 1 2 2 2 1 2 ...
#$ Cold_ischemia_time_minutes: int  317 307 723 449 273 244 353 246 270 231 ...
#$ AR_status                 : int  0 1 0 0 1 0 1 1 0 0 ...
#$ Eplets_total_HLAI         : int  10 19 19 19 10 11 11 29 13 26 ...
#$ Eplets_total_HLAII        : int  4 7 8 4 12 15 9 23 10 11 ...
#$ AR_Cox_time               : num  86.63 0.197 83.605 83.605 0.723 ...
#$ Transplantation_year      : int  2013 2013 2014 2014 2014 2014 2014 2014 2014 2014 ...
#$ Autoimmune_status         : int  1 1 0 0 1 0 1 0 1 0 ...
#$ CNI_type_initial_status   : int  1 1 1 1 1 1 1 1 1 1 ...
#$ LR_Cox_time               : num  86.6 85 83.6 83.6 13.5 ...
#$ Late_rejection_status     : int  0 0 0 0 1 0 0 0 0 0 ...
#$ Graft_loss_status         : int  0 0 0 0 0 0 1 0 0 0 ...
#$ graft_loss_months         : num  86.6 85 83.6 83.6 83.6 ...
#$ Death_status              : int  0 0 0 0 0 0 0 0 0 0 ...
#$ Death_Cox_time_months     : num  86.6 85 83.6 83.6 83.6 ...
#$ Mm_all                    : int  4782 5091 4860 4518 4647 4818 4851 5051 4909 4923 ...
#$ Mm_trans_secr             : int  1636 1709 1627 1460 1583 1666 1629 1663 1657 1632 ...
#$ Mm_transmemb              : int  1202 1238 1211 1058 1208 1233 1246 1226 1247 1212 ...
#$ Mm_liver                  : int  578 559 494 485 491 497 529 505 517 538 ...
#$ quartile                  : int  2 4 3 1 1 3 3 4 4 4 ...
#$ quart_t_secr              : int  3 4 2 1 1 3 2 3 3 3 ...
#$ quart_transmemb           : int  2 3 2 1 2 3 3 3 3 2 ...
#$ quart_liver               : int  4 4 1 1 1 1 3 2 2 3 ...
#- attr(*, ".internal.selfref")=<externalptr> 

################################################################################
### The survival analysis: the mismatch sum association to late rejection 
### event

# All
LR_cox_all <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~ Mm_all + 
                      R_sex + D_sex + R_age + D_age + 
                      Cold_ischemia_time_minutes +
                      Eplets_total_HLAI + Eplets_total_HLAII + 
                      Transplantation_year + Autoimmune_status + 
                      CNI_type_initial_status,
                    data = R_cov_mm_liver)

LR_cox_all_sum <- summary(LR_cox_all) %>% coef()
LR_cox_all_COEF <- bind_cols(LR_cox_all_sum, 
                             as.data.frame(exp(confint(LR_cox_all))))

# Edit the data frame
LR_cox_all_COEF <- round(LR_cox_all_COEF, digits = 3)
LR_cox_all_COEF$`HR(95%_CI)`<- paste0(LR_cox_all_COEF$`exp(coef)`, "(",
                                      LR_cox_all_COEF$`2.5 %`, "-",  
                                      LR_cox_all_COEF$`97.5 %`, ")" )
LR_cox_all_COEF <- rename(LR_cox_all_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5%`= `2.5 %`,
                          `97.5%` = `97.5 %`) %>%
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_cox_all_COEF,
            "./results/Missense_variants/LR_Cox_all_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyze the quartiles of the mismatch data in all missense variants

# Cox proportional hazards model for adjusted data
LR_Cox_all_quartiles <- coxph(Surv(LR_Cox_time, 
                                   Late_rejection_status) ~ quartile + R_sex +
                                D_sex + R_age + D_age + 
                                Cold_ischemia_time_minutes + 
                                Eplets_total_HLAI + 
                                Eplets_total_HLAII + Transplantation_year + 
                                Autoimmune_status + 
                                CNI_type_initial_status,
                              data = R_cov_mm_liver)

LR_Cox_all_quart_sum <- summary(LR_Cox_all_quartiles) %>% coef()
LR_Cox_all_quart_COEF <- bind_cols(LR_Cox_all_quart_sum,
                                   as.data.frame(exp(confint(LR_Cox_all_quartiles))))
# Edit the data frame 
LR_Cox_all_quart_COEF <- round(LR_Cox_all_quart_COEF, digits = 3)
LR_Cox_all_quart_COEF$`HR(95%_CI)`<- paste0(LR_Cox_all_quart_COEF$`exp(coef)`, 
                                            "(",
                                            LR_Cox_all_quart_COEF$`2.5 %`, "-",  
                                            LR_Cox_all_quart_COEF$`97.5 %`, ")")
LR_Cox_all_quart_COEF <- rename(LR_Cox_all_quart_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_Cox_all_quart_COEF,
            "./results/Missense_variants/LR_Mm_missense_adjusted_all_quartiles_COEF_CI",
            quote = F, row.names = F)

# Produce Kapplan-Meier plot with risk table for the manuscript
LR_all_quart_KM <- ggsurvplot(
  fit = survfit(Surv(LR_Cox_time, Late_rejection_status) ~ quartile, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Late rejection free survival",
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
LR_all_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
LR_all_quart_KM$plot <- LR_all_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.812, 95% CI 0.674-0.977 \n P-value 0.028",
                    size = 6)
# Customize risk table title using ggplot verbs
LR_all_quart_KM$table <- LR_all_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
LR_all_quart_KM

jpeg('./results/Missense_variants/LR_Missense_variant_analyses_all_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(LR_all_quart_KM)
dev.off()
###############################################
## Transmembrane and secretory
LR_cox_trans_secr <- coxph(Surv(LR_Cox_time, 
                                Late_rejection_status) ~ Mm_trans_secr + 
                             R_sex + D_sex + R_age + D_age + 
                             Cold_ischemia_time_minutes +
                             Eplets_total_HLAI + Eplets_total_HLAII + 
                             Transplantation_year + Autoimmune_status + 
                             CNI_type_initial_status,
                           data = R_cov_mm_liver)

LR_cox_trans_secr_sum <- summary(LR_cox_trans_secr) %>% coef()
LR_cox_trans_secr_COEF <- bind_cols(LR_cox_trans_secr_sum, 
                                    as.data.frame(exp(confint(LR_cox_trans_secr))))

# Modifying the data frame
LR_cox_trans_secr_COEF <- round(LR_cox_trans_secr_COEF, digits = 3)
LR_cox_trans_secr_COEF$`HR(95%_CI)`<- paste0(LR_cox_trans_secr_COEF$`exp(coef)`,
                                             "(", LR_cox_trans_secr_COEF$`2.5 %`,
                                             "-", LR_cox_trans_secr_COEF$`97.5 %`,
                                             ")")
LR_cox_trans_secr_COEF <- rename(LR_cox_trans_secr_COEF, `p_value` = `Pr(>|z|)`,
                                 `HR` = `exp(coef)`,
                                 `2.5%`= `2.5 %`,
                                 `97.5%` = `97.5 %`) %>%
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_cox_trans_secr_COEF,
            "./results/Missense_variants/LR_Cox_trans_secr_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyzing the quartiles of the mismatch data in transmembrane and secretory
## missense variants

# Cox proportional hazards model for adjusted data
LR_Cox_t_secr_quart <- coxph(Surv(LR_Cox_time, 
                                  Late_rejection_status) ~ quart_t_secr + R_sex +
                               D_sex + R_age + D_age + Cold_ischemia_time_minutes
                             + Eplets_total_HLAI + 
                               Eplets_total_HLAII + Transplantation_year + 
                               Autoimmune_status + 
                               CNI_type_initial_status,
                             data = R_cov_mm_liver)

LR_Cox_t_secr_quart_sum <- summary(LR_Cox_t_secr_quart) %>% coef()
LR_Cox_t_secr_quart_COEF <- bind_cols(LR_Cox_t_secr_quart_sum,
                                      as.data.frame(exp(confint(LR_Cox_t_secr_quart))))

# Modifying the data frame 
LR_Cox_t_secr_quart_COEF <- round(LR_Cox_t_secr_quart_COEF, digits = 3)
LR_Cox_t_secr_quart_COEF$`HR(95%_CI)`<- paste0(LR_Cox_t_secr_quart_COEF$`exp(coef)`,
                                               "(",
                                               LR_Cox_t_secr_quart_COEF$`2.5 %`,
                                               "-",
                                               LR_Cox_t_secr_quart_COEF$`97.5 %`,
                                               ")")
LR_Cox_t_secr_quart_COEF <- rename(LR_Cox_t_secr_quart_COEF, `p_value` = `Pr(>|z|)`,
                                   `HR` = `exp(coef)`,
                                   `2.5%`= `2.5 %`,
                                   `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_Cox_t_secr_quart_COEF,
            "./results/Missense_variants/LR_Cox_trans_secr_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the manuscript
LR_t_secr_quart_KM <- ggsurvplot(
  fit = survfit(Surv(LR_Cox_time, Late_rejection_status) ~ quart_t_secr, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Late rejection free survival",
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
LR_t_secr_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
LR_t_secr_quart_KM$plot <- LR_t_secr_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.890, 95% CI 0.740-1.071 \n P-value 0.218",
                    size = 6)
# Customize risk table title using ggplot verbs
LR_t_secr_quart_KM$table <- LR_t_secr_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
LR_t_secr_quart_KM

jpeg('./results/Missense_variants/LR_Missense_variant_analyses_transmembrane_secr_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(LR_t_secr_quart_KM)
dev.off()
#####################################################
## Transmembrane only
LR_cox_transmemb <- coxph(Surv(LR_Cox_time, 
                               Late_rejection_status) ~ Mm_transmemb + 
                            R_sex + D_sex + R_age + D_age + 
                            Cold_ischemia_time_minutes +
                            Eplets_total_HLAI + Eplets_total_HLAII + 
                            Transplantation_year + Autoimmune_status + 
                            CNI_type_initial_status,
                          data = R_cov_mm_liver)

LR_cox_transmemb_sum <- summary(LR_cox_transmemb) %>% coef()
LR_cox_transmemb_COEF <- bind_cols(LR_cox_transmemb_sum, 
                                   as.data.frame(exp(confint(LR_cox_transmemb))))

# Modifying the data frame
LR_cox_transmemb_COEF <- round(LR_cox_transmemb_COEF, digits = 3)
LR_cox_transmemb_COEF$`HR(95%_CI)`<- paste0(LR_cox_transmemb_COEF$`exp(coef)`,
                                            "(", LR_cox_transmemb_COEF$`2.5 %`,
                                            "-", LR_cox_transmemb_COEF$`97.5 %`, 
                                            ")")
LR_cox_transmemb_COEF <- rename(LR_cox_transmemb_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_cox_transmemb_COEF,
            "./results/Missense_variants/LR_Cox_transmemb_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyzing the quartiles of the mismatch data in transmembrane
## missense variants

# Cox proportional hazards model for adjusted data
LR_cox_trans_quart <- coxph(Surv(LR_Cox_time, 
                                 Late_rejection_status) ~ quart_transmemb + R_sex +
                              D_sex + R_age + D_age + Cold_ischemia_time_minutes
                            + Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status,
                            data = R_cov_mm_liver)

LR_Cox_trans_quart_sum <- summary(LR_cox_trans_quart) %>% coef()
LR_Cox_trans_quart_COEF <- bind_cols(LR_Cox_trans_quart_sum,
                                     as.data.frame(exp(confint(LR_cox_trans_quart))))

# Modifying the data frame 
LR_Cox_trans_quart_COEF <- round(LR_Cox_trans_quart_COEF, digits = 3)
LR_Cox_trans_quart_COEF$`HR(95%_CI)`<- paste0(LR_Cox_trans_quart_COEF$`exp(coef)`,
                                              "(", LR_Cox_trans_quart_COEF$`2.5 %`,
                                              "-", LR_Cox_trans_quart_COEF$`97.5 %`, 
                                              ")")
LR_Cox_trans_quart_COEF <- rename(LR_Cox_trans_quart_COEF, 
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_Cox_trans_quart_COEF,
            "./results/Missense_variants/LR_Cox_transmemb_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the manuscript
LR_trans_quart_KM <- ggsurvplot(
  fit = survfit(Surv(LR_Cox_time, Late_rejection_status) ~ quart_transmemb, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Late rejection free survival",
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
LR_trans_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
LR_trans_quart_KM$plot <- LR_trans_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.859, 95% CI 0.714-1.034 \n P-value 0.109",
                    size = 6)
# Customize risk table title using ggplot verbs
LR_trans_quart_KM$table <- LR_trans_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
LR_trans_quart_KM

jpeg('./results/Missense_variants/LR_Missense_variant_analyses_transmembrane_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(LR_trans_quart_KM)
dev.off()

######################################################
## Liver related
LR_cox_liver <- coxph(Surv(LR_Cox_time,
                           Late_rejection_status) ~ Mm_liver + 
                        R_sex + D_sex + R_age + D_age + 
                        Cold_ischemia_time_minutes +
                        Eplets_total_HLAI + Eplets_total_HLAII + 
                        Transplantation_year + Autoimmune_status + 
                        CNI_type_initial_status,
                      data = R_cov_mm_liver)

LR_cox_liver_sum <- summary(LR_cox_liver) %>% coef()
LR_cox_liver_COEF <- bind_cols(LR_cox_liver_sum, 
                               as.data.frame(exp(confint(LR_cox_liver))))

# Modifying the data frame
LR_cox_liver_COEF <- round(LR_cox_liver_COEF, digits = 3)
LR_cox_liver_COEF$`HR(95%_CI)`<- paste0(LR_cox_liver_COEF$`exp(coef)`,
                                        "(", LR_cox_liver_COEF$`2.5 %`,
                                        "-", LR_cox_liver_COEF$`97.5 %`, 
                                        ")")
LR_cox_liver_COEF <- rename(LR_cox_liver_COEF, `p_value` = `Pr(>|z|)`,
                            `HR` = `exp(coef)`,
                            `2.5%`= `2.5 %`,
                            `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_cox_liver_COEF,
            "./results/Missense_variants/LR_Cox_liver_related_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyzing the quartiles of the mismatch data in liver-related
## missense variants

# Cox proportional hazards model for adjusted data
LR_cox_liver_quart <- coxph(Surv(LR_Cox_time,
                                 Late_rejection_status) ~ quart_liver + R_sex +
                              D_sex + R_age + D_age + Cold_ischemia_time_minutes
                            + Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status,
                            data = R_cov_mm_liver)

LR_Cox_liver_quart_sum <- summary(LR_cox_liver_quart) %>% coef()
LR_Cox_liver_quart_COEF <- bind_cols(LR_Cox_liver_quart_sum,
                                     as.data.frame(exp(confint(LR_cox_liver_quart))))

# Modidying the data frame
LR_Cox_liver_quart_COEF <- round(LR_Cox_liver_quart_COEF, digits = 3)
LR_Cox_liver_quart_COEF$`HR(95%_CI)`<- paste0(LR_Cox_liver_quart_COEF$`exp(coef)`,
                                              "(", LR_Cox_liver_quart_COEF$`2.5 %`,
                                              "-", LR_Cox_liver_quart_COEF$`97.5 %`, 
                                              ")")
LR_Cox_liver_quart_COEF <- rename(LR_Cox_liver_quart_COEF,
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(LR_Cox_liver_quart_COEF,
            "./results/Missense_variants/LR_Cox_liver_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the manuscript
LR_liver_quart_KM <- ggsurvplot(
  fit = survfit(Surv(LR_Cox_time, Late_rejection_status) ~ quart_liver, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Late rejection free survival",
  font.x = c(16, face = "bold"),
  font.y = c(16, face = "bold"),
  font.tickslab = c(16),
  risk.table = TRUE, fontsize = 5, 
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 438-498", "Q2: MM sum 499-518",
                  "Q3: MM sum 519-539", "Q4: MM sum 540-629"),
  font.legend = c(16),
  legend = c(0.20,0.35))
LR_liver_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
LR_liver_quart_KM$plot <- LR_liver_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.937, 95% CI 0.779-1.127 \n P-value 0.489",
                    size = 6)
# Customize risk table title using ggplot verbs
LR_liver_quart_KM$table <- LR_liver_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
LR_liver_quart_KM

jpeg('./results/Missense_variants/LR_Missense_variant_analyses_liver_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(LR_liver_quart_KM)
dev.off()

################################################################################
### Performing multiple test corrections for LR for mismatch sum and quartiles
## Mm sums

# Creating combined results
LR_results <- list(LR_cox_all_COEF, LR_cox_trans_secr_COEF, 
                   LR_cox_transmemb_COEF, LR_cox_liver_COEF)
names(LR_results) <- c("LR_cox_all_COEF", "LR_cox_trans_secr_COEF", 
                       "LR_cox_transmemb_COEF", "LR_cox_liver_COEF") 

# Extract Mm covariates and p values
LR_results <- map((LR_results),
                  function(x) {
                    DATA <- x[1,] 
                    DATA <- select(DATA, covariates, p_value)
                    return(DATA)
                  })

LR_results <- bind_rows(LR_results)

# Performing multiple test correction
LR_results$Bonferroni <- p.adjust(LR_results[,2], method = "bonferroni", n = 4)
LR_results$Holm <- p.adjust(LR_results[,2], method = "holm", n = 4)
LR_results$FDR <- p.adjust(LR_results[,2], method = "fdr", n = 4)

write.table(LR_results, 
            "results/Missense_variants/LR_Cox_Mm_mismatch_multiple_testing_results",
            quote = F, row.names = F)

## Mm quartiles
# Creating combined results
LR_quart_results <- list(LR_Cox_all_quart_COEF, LR_Cox_t_secr_quart_COEF, 
                         LR_Cox_trans_quart_COEF, LR_Cox_liver_quart_COEF)
names(LR_quart_results) <- c("LR_Cox_all_quart_COEF", "LR_Cox_t_secr_quart_COEF", 
                             "LR_Cox_trans_quart_COEF", "LR_Cox_liver_quart_COEF") 

# Extract Mm covariates and p values
LR_quart_results <- map((LR_quart_results),
                        function(x) {
                          DATA <- x[1,] 
                          DATA <- select(DATA, covariates, p_value)
                          return(DATA)
                        })

LR_quart_results <- bind_rows(LR_quart_results)

LR_quart_results$Bonferroni <- p.adjust(LR_quart_results[,2], 
                                        method = "bonferroni", n = 4)
LR_quart_results$Holm <- p.adjust(LR_quart_results[,2], method = "holm", n = 4)
LR_quart_results$FDR <- p.adjust(LR_quart_results[,2], method = "fdr", n = 4)

write.table(LR_quart_results, 
            "results/Missense_variants/LR_Cox_Mm_quartiles_multiple_testing_results",
            quote = F, row.names = F)