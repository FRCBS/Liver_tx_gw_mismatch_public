###############################################################################
### Perform missense variant analyses using adjusted Cox proportional 
### hazards model for time to graft loss (GL)
###############################################################################
### General information

# Prerequisites:
# 1) Run script LR_missense_variant_analyses.R
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

str(R_cov_mm_liver)
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
#$ GL_Cox_time               : num  86.6 85 83.6 83.6 13.5 ...
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

###############################################################################
### The survival analysis: the mismatch sum association to time to acute 
### rejection event

## All proteins

# Cox proportional hazards model for adjusted data
GL_cox_all <- coxph(Surv(graft_loss_months, Graft_loss_status) ~ Mm_all + 
                      R_sex + D_sex + R_age + D_age + 
                      Cold_ischemia_time_minutes +
                      Eplets_total_HLAI + Eplets_total_HLAII + 
                      Transplantation_year + Autoimmune_status + 
                      CNI_type_initial_status + AR_status,
                    data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
GL_cox_all_sum <- summary(GL_cox_all) %>% coef()
GL_cox_all_COEF <- bind_cols(GL_cox_all_sum, as.data.frame(exp(confint(GL_cox_all))))
GL_cox_all_COEF <- round(GL_cox_all_COEF, digits = 3)
GL_cox_all_COEF$`HR(95%_CI)`<- paste0(GL_cox_all_COEF$`exp(coef)`, "(",
                                 GL_cox_all_COEF$`2.5 %`, "-",  
                                 GL_cox_all_COEF$`97.5 %`, ")" )

GL_cox_all_COEF <- rename(GL_cox_all_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5%`= `2.5 %`,
                          `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_cox_all_COEF,
            "./results/Missense_variants/GL_Cox_all_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## All proteins

## Analyze the quartiles of the mismatch data in all missense variants

# Cox proportional hazards model for adjusted data
GL_Cox_all_quartiles <- coxph(Surv(graft_loss_months, 
                                   Graft_loss_status) ~ quartile + R_sex +
                                D_sex + R_age + D_age + 
                                Cold_ischemia_time_minutes + 
                                Eplets_total_HLAI + 
                                Eplets_total_HLAII + Transplantation_year + 
                                Autoimmune_status + 
                                CNI_type_initial_status + AR_status,
                              data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
GL_Cox_all_quart_sum <- summary(GL_Cox_all_quartiles) %>% coef()
GL_Cox_all_quart_COEF <- bind_cols(GL_Cox_all_quart_sum,
                                   as.data.frame(exp(confint(GL_Cox_all_quartiles))))
GL_Cox_all_quart_COEF <- round(GL_Cox_all_quart_COEF, digits = 3)
GL_Cox_all_quart_COEF$`HR(95%_CI)`<- paste0(GL_Cox_all_quart_COEF$`exp(coef)`, 
                                            "(",
                                      GL_Cox_all_quart_COEF$`2.5 %`, "-",  
                                      GL_Cox_all_quart_COEF$`97.5 %`, ")")

GL_Cox_all_quart_COEF <- rename(GL_Cox_all_quart_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_Cox_all_quart_COEF,
            "./results/Missense_variants/GL_Mm_missense_adjusted_all_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
GL_all_quart_KM <- ggsurvplot(
  fit = survfit(Surv(graft_loss_months, Graft_loss_status) ~ quartile, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Graft loss free survival",
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
GL_all_quart_KM

# Add HR and 95% CI
GL_all_quart_KM$plot <- GL_all_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.866, 95% CI 0.649-1.156 \n P-value 0.329",
                    size = 6)
# Customize risk table title using ggplot verbs
GL_all_quart_KM$table <- GL_all_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
GL_all_quart_KM

jpeg('./results/Missense_variants/GL_Missense_variant_analyses_all_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(GL_all_quart_KM)
dev.off()

################################################################################
## Transmembrane and secreted

# Cox proportional hazards model for adjusted data
GL_cox_trans_secr <- coxph(Surv(graft_loss_months, 
                                Graft_loss_status) ~ Mm_trans_secr + 
                             R_sex + D_sex + R_age + D_age + 
                             Cold_ischemia_time_minutes +
                             Eplets_total_HLAI + Eplets_total_HLAII + 
                             Transplantation_year + Autoimmune_status + 
                             CNI_type_initial_status + AR_status,
                           data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
GL_cox_trans_secr_sum <- summary(GL_cox_trans_secr) %>% coef()
GL_cox_trans_secr_COEF <- bind_cols(GL_cox_trans_secr_sum, 
                                    as.data.frame(exp(confint(GL_cox_trans_secr))))
GL_cox_trans_secr_COEF <- round(GL_cox_trans_secr_COEF, digits = 3)
GL_cox_trans_secr_COEF$`HR(95%_CI)`<- paste0(GL_cox_trans_secr_COEF$`exp(coef)`,
                                             "(", GL_cox_trans_secr_COEF$`2.5 %`,
                                             "-", GL_cox_trans_secr_COEF$`97.5 %`,
                                             ")")

GL_cox_trans_secr_COEF <- rename(GL_cox_trans_secr_COEF, `p_value` = `Pr(>|z|)`,
                                 `HR` = `exp(coef)`,
                                 `2.5%`= `2.5 %`,
                                 `97.5%` = `97.5 %`) %>%
  rownames_to_column("covariates") %>%
select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value) 

write.table(GL_cox_trans_secr_COEF,
            "./results/Missense_variants/GL_Cox_trans_secr_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Transmembrane and secreted

## Analyze the quartiles of the mismatch data in transmembrane and secreted
## missense variants

# Cox proportional hazards model for adjusted data
GL_Cox_t_secr_quart <- coxph(Surv(graft_loss_months, 
                                  Graft_loss_status) ~ quart_t_secr + R_sex +
                               D_sex + R_age + D_age + 
                               Cold_ischemia_time_minutes + Eplets_total_HLAI + 
                               Eplets_total_HLAII + Transplantation_year + 
                               Autoimmune_status + 
                               CNI_type_initial_status + AR_status,
                             data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
GL_Cox_t_secr_quart_sum <- summary(GL_Cox_t_secr_quart) %>% coef()
GL_Cox_t_secr_quart_COEF <- bind_cols(GL_Cox_t_secr_quart_sum,
                                      as.data.frame(exp(confint(GL_Cox_t_secr_quart))))
GL_Cox_t_secr_quart_COEF <- round(GL_Cox_t_secr_quart_COEF, digits = 3)
GL_Cox_t_secr_quart_COEF$`HR(95%_CI)`<- paste0(GL_Cox_t_secr_quart_COEF$`exp(coef)`,
                                             "(",
                                             GL_Cox_t_secr_quart_COEF$`2.5 %`,
                                             "-",
                                             GL_Cox_t_secr_quart_COEF$`97.5 %`,
                                             ")")

GL_Cox_t_secr_quart_COEF <- rename(GL_Cox_t_secr_quart_COEF, `p_value` = `Pr(>|z|)`,
                                   `HR` = `exp(coef)`,
                                   `2.5%`= `2.5 %`,
                                   `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_Cox_t_secr_quart_COEF,
            "./results/Missense_variants/GL_Cox_trans_secr_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
GL_t_secr_quart_KM <- ggsurvplot(
  fit = survfit(Surv(graft_loss_months, Graft_loss_status) ~ quart_t_secr, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Graft loss free survival",
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
GL_t_secr_quart_KM

# Add HR and 95% CI
GL_t_secr_quart_KM$plot <- GL_t_secr_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.747, 95% CI 0.559-1.000 \n P-value 0.050",
                    size = 6)
# Customize risk table title using ggplot verbs
GL_t_secr_quart_KM$table <- GL_t_secr_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
GL_t_secr_quart_KM

jpeg('./results/Missense_variants/GL_Missense_variant_analyses_transmembrane_secr_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(GL_t_secr_quart_KM)
dev.off()

###############################################################################
## Transmembrane only

# Cox proportional hazards model for adjusted data
GL_cox_transmemb <- coxph(Surv(graft_loss_months, 
                               Graft_loss_status) ~ Mm_transmemb + 
                            R_sex + D_sex + R_age + D_age + 
                            Cold_ischemia_time_minutes +
                            Eplets_total_HLAI + Eplets_total_HLAII + 
                            Transplantation_year + Autoimmune_status + 
                            CNI_type_initial_status + AR_status,
                          data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
GL_cox_transmemb_sum <- summary(GL_cox_transmemb) %>% coef()
GL_cox_transmemb_COEF <- bind_cols(GL_cox_transmemb_sum, 
                                   as.data.frame(exp(confint(GL_cox_transmemb))))
GL_cox_transmemb_COEF <- round(GL_cox_transmemb_COEF, digits = 3)
GL_cox_transmemb_COEF$`HR(95%_CI)`<- paste0(GL_cox_transmemb_COEF$`exp(coef)`,
                                            "(", GL_cox_transmemb_COEF$`2.5 %`,
                                            "-", GL_cox_transmemb_COEF$`97.5 %`, 
                                            ")")

GL_cox_transmemb_COEF <- rename(GL_cox_transmemb_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_cox_transmemb_COEF,
            "./results/Missense_variants/GL_Cox_transmemb_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Transmembrane only

## Analyze the quartiles of the mismatch data in transmembrane missense 
## variants

# Cox proportional hazards model for adjusted data
GL_cox_trans_quart <- coxph(Surv(graft_loss_months, 
                                 Graft_loss_status) ~ quart_transmemb + R_sex +
                              D_sex + R_age + D_age + 
                              Cold_ischemia_time_minutes + 
                              Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status + AR_status,
                            data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
GL_Cox_trans_quart_sum <- summary(GL_cox_trans_quart) %>% coef()
GL_Cox_trans_quart_COEF <- bind_cols(GL_Cox_trans_quart_sum,
                                     as.data.frame(exp(confint(GL_cox_trans_quart))))
GL_Cox_trans_quart_COEF <- round(GL_Cox_trans_quart_COEF, digits = 3)
GL_Cox_trans_quart_COEF$`HR(95%_CI)`<- paste0(GL_Cox_trans_quart_COEF$`exp(coef)`,
                                            "(", GL_Cox_trans_quart_COEF$`2.5 %`,
                                            "-", GL_Cox_trans_quart_COEF$`97.5 %`, 
                                            ")")

GL_Cox_trans_quart_COEF <- rename(GL_Cox_trans_quart_COEF, 
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_Cox_trans_quart_COEF,
            "./results/Missense_variants/GL_Cox_transmemb_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
GL_trans_quart_KM <- ggsurvplot(
  fit = survfit(Surv(graft_loss_months, Graft_loss_status) ~ quart_transmemb, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Graft loss free survival",
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
GL_trans_quart_KM

# Add HR and 95% CI
GL_trans_quart_KM$plot <- GL_trans_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.740, 95% CI 0.551-0.994 \n P-value 0.046",
                    size = 6)
# Customize risk table title using ggplot verbs
GL_trans_quart_KM$table <- GL_trans_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
GL_trans_quart_KM

jpeg('./results/Missense_variants/GL_Missense_variant_analyses_transmembrane_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(GL_trans_quart_KM)
dev.off()

###############################################################################
## Liver-related

# Cox proportional hazards model for adjusted data
GL_cox_liver <- coxph(Surv(graft_loss_months, Graft_loss_status) ~ Mm_liver + 
                        R_sex + D_sex + R_age + D_age + 
                        Cold_ischemia_time_minutes +
                        Eplets_total_HLAI + Eplets_total_HLAII + 
                        Transplantation_year + Autoimmune_status + 
                        CNI_type_initial_status + AR_status,
                      data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
GL_cox_liver_sum <- summary(GL_cox_liver) %>% coef()
GL_cox_liver_COEF <- bind_cols(GL_cox_liver_sum, 
                               as.data.frame(exp(confint(GL_cox_liver))))
GL_cox_liver_COEF <- round(GL_cox_liver_COEF, digits = 3)
GL_cox_liver_COEF$`HR(95%_CI)`<- paste0(GL_cox_liver_COEF$`exp(coef)`,
                                              "(", GL_cox_liver_COEF$`2.5 %`,
                                              "-", GL_cox_liver_COEF$`97.5 %`, 
                                              ")")

GL_cox_liver_COEF <- rename(GL_cox_liver_COEF, `p_value` = `Pr(>|z|)`,
                            `HR` = `exp(coef)`,
                            `2.5%`= `2.5 %`,
                            `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_cox_liver_COEF,
            "./results/Missense_variants/GL_Cox_liver_related_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## Liver-related

## Analyze the quartiles of the mismatch data in liver-related missense 
## variants

# Cox proportional hazards model for adjusted data
GL_cox_liver_quart <- coxph(Surv(graft_loss_months, 
                                 Graft_loss_status) ~ quart_liver + R_sex +
                              D_sex + R_age + D_age + 
                              Cold_ischemia_time_minutes + Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status + AR_status,
                            data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
GL_Cox_liver_quart_sum <- summary(GL_cox_liver_quart) %>% coef()
GL_Cox_liver_quart_COEF <- bind_cols(GL_Cox_liver_quart_sum,
                                     as.data.frame(exp(confint(GL_cox_liver_quart))))
GL_Cox_liver_quart_COEF <- round(GL_Cox_liver_quart_COEF, digits = 3)
GL_Cox_liver_quart_COEF$`HR(95%_CI)`<- paste0(GL_Cox_liver_quart_COEF$`exp(coef)`,
                                        "(", GL_Cox_liver_quart_COEF$`2.5 %`,
                                        "-", GL_Cox_liver_quart_COEF$`97.5 %`, 
                                        ")")

GL_Cox_liver_quart_COEF <- rename(GL_Cox_liver_quart_COEF,
                                  `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>%
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(GL_Cox_liver_quart_COEF,
            "./results/Missense_variants/GL_Cox_liver_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for quartile data
GL_liver_quart_KM <- ggsurvplot(
  fit = survfit(Surv(graft_loss_months, Graft_loss_status) ~ quart_liver, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Graft loss free survival",
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
GL_liver_quart_KM

# HR and 95% CI
GL_liver_quart_KM$plot <- GL_liver_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.163, 95% CI 0.869-1.556 \n P-value 0.309",
                    size = 6)
# Customize risk table title using ggplot verbs
GL_liver_quart_KM$table <- GL_liver_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
GL_liver_quart_KM

jpeg('./results/Missense_variants/GL_Missense_variant_analyses_liver_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(GL_liver_quart_KM)
dev.off()

################################################################################
### Perform multiple test corrections for GL for mismatch sum and quartiles

## Mismatch sums

# Create combined results
GL_results <- list(GL_cox_all_COEF, GL_cox_trans_secr_COEF, 
                   GL_cox_transmemb_COEF, GL_cox_liver_COEF)
names(GL_results) <- c("GL_cox_all_COEF", "GL_cox_trans_secr_COEF", 
                       "GL_cox_transmemb_COEF", "GL_cox_liver_COEF") 

# Extract mismatch covariates and p values
GL_results <- map((GL_results),
                function(x) {
                  DATA <- x[1,] 
                  DATA <- select(DATA, covariates, p_value)
                  return(DATA)
                })

GL_results <- bind_rows(GL_results)

# Perform multiple test correction
GL_results$Bonferroni <- p.adjust(GL_results[,2], method = "bonferroni", n = 4)
GL_results$Holm <- p.adjust(GL_results[,2], method = "holm", n = 4)
GL_results$FDR <- p.adjust(GL_results[,2], method = "fdr", n = 4)

write.table(GL_results, 
            "results/Missense_variants/GL_Cox_Mm_mismatch_multiple_testing_results",
            quote = F, row.names = F)

#######################################
## Mismatch quartiles

# Create combined results
GL_quart_results <- list(GL_Cox_all_quart_COEF, GL_Cox_t_secr_quart_COEF, 
                   GL_Cox_trans_quart_COEF, GL_Cox_liver_quart_COEF)
names(GL_quart_results) <- c("GL_Cox_all_quart_COEF", "GL_Cox_t_secr_quart_COEF", 
                       "GL_Cox_trans_quart_COEF", "GL_Cox_liver_quart_COEF") 

# Extract mismatch covariates and p values
GL_quart_results <- map((GL_quart_results),
                  function(x) {
                    DATA <- x[1,] 
                    DATA <- select(DATA, covariates, p_value)
                    return(DATA)
                  })

GL_quart_results <- bind_rows(GL_quart_results)
GL_quart_results$Bonferroni <- p.adjust(GL_quart_results[,2], 
                                        method = "bonferroni", n = 4)
GL_quart_results$Holm <- p.adjust(GL_quart_results[,2], method = "holm", n = 4)
GL_quart_results$FDR <- p.adjust(GL_quart_results[,2], method = "fdr", n = 4)

write.table(GL_quart_results, 
            "results/Missense_variants/GL_Cox_Mm_quartiles_multiple_testing_results",
            quote = F, row.names = F)
###############################################################################