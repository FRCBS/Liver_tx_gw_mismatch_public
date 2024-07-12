###############################################################################
### Perform missense variant analyses using adjusted Cox proportional 
### hazards model for time to acute rejection (AR)
###############################################################################
### General information

# Prerequisites:
# 1) Run script 04_Missense_variant_mismatches.R
  # data/Missense_variants/R_covariates_mm_liver.txt
# 2) folder results/Missense_variants

###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

###############################################################################
## Import covariate file including all required covariates and mismatch sums
## for all four protein groups 
R_cov_mm_liver <- fread("data/Missense_variants/R_covariates_mm_liver.txt")

str(R_cov_mm_liver)
# Classes ?data.table? and 'data.frame':	666 obs. of  25 variables:
#$ Pair_number               : int  539 540 542 543 541 544 545 546 547 548 ...
#$ R_pseudo                  : chr  "R_pseudo1" "R_pseudo2" "R_pseudo3" "R_pseudo4" ...
#$ D_pseudo                  : chr  "D_pseudo1" "D_pseudo2" "D_pseudo3" "D_pseudo4" ...
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
#$ Mm_transm_secr             : int  1636 1709 1627 1460 1583 1666 1629 1663 1657 1632 ...
#$ Mm_transm              : int  1202 1238 1211 1058 1208 1233 1246 1226 1247 1212 ...
#$ Mm_liver                  : int  578 559 494 485 491 497 529 505 517 538 ...
#- attr(*, ".internal.selfref")=<externalptr> 

###############################################################################
### The survival analysis: the mismatch sum association to time to acute 
### rejection event

## All proteins

# Cox proportional hazards model for adjusted data
AR_cox_all <- coxph(Surv(AR_Cox_time, AR_status) ~ Mm_all + 
                   R_sex + D_sex + R_age + D_age + Cold_ischemia_time_minutes +
                   Eplets_total_HLAI + Eplets_total_HLAII + 
                   Transplantation_year + Autoimmune_status + 
                   CNI_type_initial_status,
                 data = R_cov_mm_liver)

# Modify Cox analysis results into data frame
AR_cox_all_sum <- summary(AR_cox_all) %>% coef()
AR_cox_all_COEF <- bind_cols(AR_cox_all_sum, as.data.frame(exp(confint(AR_cox_all))))
AR_cox_all_COEF <- round(AR_cox_all_COEF, digits = 3) 
AR_cox_all_COEF$`HR(95%_CI)` <- paste0(AR_cox_all_COEF$`exp(coef)`, "(",
                                      AR_cox_all_COEF$`2.5 %`, "-", 
                                      AR_cox_all_COEF$`97.5 %`,
                                      ")")

AR_cox_all_COEF <- rename(AR_cox_all_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5_%`= `2.5 %`,
                          `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_cox_all_COEF,
            "./results/Missense_variants/AR_Cox_all_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

#######################################
## All proteins 
## Analyze the quartiles of the mismatch data in all missense variants
## Divide the mismatch sums into quartiles
quantile(R_cov_mm_liver$Mm_all)
#      0%     25%     50%     75%    100% 
# 4461.00 4721.00 4798.00 4871.75 5604.00 

# Create a new column from quartiles (1 = the lowest quartile, 2 = the second, 
# 3 = the third, 4 = the highest quartile)
R_cov_mm_liver$quartile <- as.integer(cut(R_cov_mm_liver$Mm_all, 
                                               quantile(R_cov_mm_liver$Mm_all, 
                                                        probs = 0:4/4), 
                                               include.lowest = TRUE))

# The number of pairs per each quartile
survdiff_quartile <- survdiff(Surv(AR_Cox_time, AR_status) ~ quartile, 
                              data = R_cov_mm_liver)

# Cox proportional hazards model for adjusted data
AR_Cox_all_quartiles <- coxph(Surv(AR_Cox_time, AR_status) ~ quartile + R_sex +
                             D_sex + R_age + D_age + Cold_ischemia_time_minutes
                           + Eplets_total_HLAI + 
                             Eplets_total_HLAII + Transplantation_year + 
                             Autoimmune_status + 
                             CNI_type_initial_status,
                           data = R_cov_mm_liver)

# Modify Cox analysis quartile results into data frame
AR_Cox_all_quart_sum <- summary(AR_Cox_all_quartiles) %>% coef()
AR_Cox_all_quart_COEF <- bind_cols(AR_Cox_all_quart_sum,
                                    as.data.frame(exp(confint(AR_Cox_all_quartiles))))
AR_Cox_all_quart_COEF <- round(AR_Cox_all_quart_COEF, digits = 3)
AR_Cox_all_quart_COEF$`HR(95%_CI)` <- paste0(AR_Cox_all_quart_COEF$`exp(coef)`, "(",
                                       AR_Cox_all_quart_COEF$`2.5 %`, "-", 
                                       AR_Cox_all_quart_COEF$`97.5 %`,
                                       ")")
AR_Cox_all_quart_COEF <- rename(AR_Cox_all_quart_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5_%`= `2.5 %`,
                                `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_Cox_all_quart_COEF,
            "./results/Missense_variants/AR_Mm_missense_adjusted_all_quartiles_COEF_CI",
            quote = F, row.names = F)

## Kapplan-Meier plot with risk table for quartile data
AR_all_quart_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ quartile, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Rejection free survival",
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
AR_all_quart_KM

# Add HR and 95% CI
AR_all_quart_KM$plot <- AR_all_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.000, 95% CI 0.999-1.001 \n P-value 0.891",
                    size = 6)
# Customize risk table title using ggplot verbs
AR_all_quart_KM$table <- AR_all_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))

AR_all_quart_KM

jpeg('results/Missense_variants/AR_Missense_variant_analyses_all_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_all_quart_KM)
dev.off()

###############################################################################

## Transmembrane and secreted
AR_cox_trans_secr <- coxph(Surv(AR_Cox_time, AR_status) ~ Mm_transm_secr + 
                   R_sex + D_sex + R_age + D_age + Cold_ischemia_time_minutes +
                   Eplets_total_HLAI + Eplets_total_HLAII + 
                   Transplantation_year + Autoimmune_status + 
                   CNI_type_initial_status,
                 data = R_cov_mm_liver)

AR_cox_trans_secr_sum <- summary(AR_cox_trans_secr) %>% coef()
AR_cox_trans_secr_COEF <- bind_cols(AR_cox_trans_secr_sum, 
                                    as.data.frame(exp(confint(AR_cox_trans_secr))))

AR_cox_trans_secr_COEF <- round(AR_cox_trans_secr_COEF, digits = 3)
AR_cox_trans_secr_COEF$`HR(95%_CI)`<- paste0(AR_cox_trans_secr_COEF$`exp(coef)`,
                                         "(",  AR_cox_trans_secr_COEF$`2.5 %`, 
                                         "-", AR_cox_trans_secr_COEF$`97.5 %`,
                                         ")") 
AR_cox_trans_secr_COEF <- rename(AR_cox_trans_secr_COEF, `p_value` = `Pr(>|z|)`,
                          `HR` = `exp(coef)`,
                          `2.5_%`= `2.5 %`,
                          `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_cox_trans_secr_COEF,
            "./results/Missense_variants/AR_Cox_trans_secr_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyzing the quartiles of the mismatch data in transmembrane and secretory
## missense variants

# Dividing the mismatch sum into quartiles
quantile(R_cov_mm_liver$Mm_transm_secr)
#   0%  25%  50%  75% 100% 
# 1460 1593 1629 1675 1952 

# Create a new column from quartiles (1 = the lowest quartile, 2 = the second, 
# 3 = the third, 4 = the highest quartile)
R_cov_mm_liver$quart_t_secr <- as.integer(cut(R_cov_mm_liver$Mm_transm_secr, 
                                              quantile(R_cov_mm_liver$Mm_transm_secr, 
                                                       probs = 0:4/4), 
                                              include.lowest = TRUE))

# Cox proportional hazards model for adjusted data
AR_Cox_t_secr_quart <- coxph(Surv(AR_Cox_time, AR_status) ~ quart_t_secr + R_sex +
                               D_sex + R_age + D_age + Cold_ischemia_time_minutes
                             + Eplets_total_HLAI + 
                               Eplets_total_HLAII + Transplantation_year + 
                               Autoimmune_status + 
                               CNI_type_initial_status,
                             data = R_cov_mm_liver)

AR_Cox_t_secr_quart_sum <- summary(AR_Cox_t_secr_quart) %>% coef()
AR_Cox_t_secr_quart_COEF <- bind_cols(AR_Cox_t_secr_quart_sum,
                                   as.data.frame(exp(confint(AR_Cox_t_secr_quart))))

# Edit the data frame 
AR_Cox_t_secr_quart_COEF <- round(AR_Cox_t_secr_quart_COEF, digits = 3)

AR_Cox_t_secr_quart_COEF$`HR(95%_CI)` <- paste0(AR_Cox_t_secr_quart_COEF$`exp(coef)`,
                                                "(", 
                                            AR_Cox_t_secr_quart_COEF$`2.5 %`,
                                            "-", 
                                            AR_Cox_t_secr_quart_COEF$`97.5 %`,
                                            ")")

AR_Cox_t_secr_quart_COEF <- rename(AR_Cox_t_secr_quart_COEF, `p_value` = `Pr(>|z|)`,
                                   `HR` = `exp(coef)`,
                                   `2.5_%`= `2.5 %`,
                                   `97.5_%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_Cox_t_secr_quart_COEF,
            "./results/Missense_variants/AR_Cox_trans_secr_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Kapplan-Meier plot with risk table for the manuscript
AR_t_secr_quart_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ quart_t_secr, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Rejection free survival",
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
AR_t_secr_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
AR_t_secr_quart_KM$plot <- AR_t_secr_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.041, 95% CI 0.934-1.161 \n P-value 0.466",
                    size = 6)
# Customize risk table title using ggplot verbs
AR_t_secr_quart_KM$table <- AR_t_secr_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
AR_t_secr_quart_KM

jpeg('./results/Missense_variants/AR_Missense_variant_analyses_transmembrane_secr_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_t_secr_quart_KM)
dev.off()

##############################################
## Transmembrane only
AR_cox_transmemb <- coxph(Surv(AR_Cox_time, AR_status) ~ Mm_transm + 
                         R_sex + D_sex + R_age + D_age + Cold_ischemia_time_minutes +
                         Eplets_total_HLAI + Eplets_total_HLAII + 
                         Transplantation_year + Autoimmune_status + 
                         CNI_type_initial_status,
                       data = R_cov_mm_liver)

AR_cox_transmemb_sum <- summary(AR_cox_transmemb) %>% coef()
AR_cox_transmemb_COEF <- bind_cols(AR_cox_transmemb_sum, 
                                    as.data.frame(exp(confint(AR_cox_transmemb))))

AR_cox_transmemb_COEF <- round(AR_cox_transmemb_COEF, digits = 3)

AR_cox_transmemb_COEF$`HR(95%_CI)`<- paste0(AR_cox_transmemb_COEF$`exp(coef)`,
                                           "(", AR_cox_transmemb_COEF$`2.5 %`,
                                           "-", AR_cox_transmemb_COEF$`97.5 %`,
                                           ")") 

AR_cox_transmemb_COEF <- rename(AR_cox_transmemb_COEF, `p_value` = `Pr(>|z|)`,
                                 `HR` = `exp(coef)`,
                                 `2.5%`= `2.5 %`,
                                 `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_cox_transmemb_COEF,
            "./results/Missense_variants/AR_Cox_transmemb_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyze the quartiles of the mismatch data in transmembrane missense variants

# Divide the mismatch sums into quartiles
quantile(R_cov_mm_liver$Mm_transm)
#  0%  25%  50%  75% 100% 
#1058 1177 1214 1250 1488


# Create a new column from quartiles (1 = the lowest quartile, 2 = the second, 
# 3 = the third, 4 = the highest quartile)
R_cov_mm_liver$quart_transmemb <- as.integer(cut(R_cov_mm_liver$Mm_transm, 
                                                 quantile(R_cov_mm_liver$Mm_transm, 
                                                          probs = 0:4/4), 
                                                 include.lowest = TRUE))

# Cox proportional hazards model for adjusted data
AR_cox_trans_quart <- coxph(Surv(AR_Cox_time, AR_status) ~ quart_transmemb + 
                              R_sex +
                              D_sex + R_age + D_age + 
                              Cold_ischemia_time_minutes +
                              Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status,
                            data = R_cov_mm_liver)

AR_Cox_trans_quart_sum <- summary(AR_cox_trans_quart) %>% coef()
AR_Cox_trans_quart_COEF <- bind_cols(AR_Cox_trans_quart_sum,
                                     as.data.frame(exp(confint(AR_cox_trans_quart))))

# Edit the data frame 
AR_Cox_trans_quart_COEF<- round(AR_Cox_trans_quart_COEF, digits = 3)

AR_Cox_trans_quart_COEF$`HR(95%_CI)` <- paste0(AR_Cox_trans_quart_COEF$`exp(coef)`,
                                           "(",
                                           AR_Cox_trans_quart_COEF$`2.5 %`, "-",
                                           AR_Cox_trans_quart_COEF$`97.5 %`,
                                           ")")

AR_Cox_trans_quart_COEF <- rename(AR_Cox_trans_quart_COEF, `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates") %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_Cox_trans_quart_COEF,
            "./results/Missense_variants/AR_Cox_transmemb_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Produce Kapplan-Meier plot with risk table for the manuscript
AR_trans_quart_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ quart_transmemb, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Rejection free survival",
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
AR_trans_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
AR_trans_quart_KM$plot <- AR_trans_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 1.021, 95% CI 0.916-1.138 \n P-value 0.708",
                    size = 6)

# Customize risk table title using ggplot verbs
AR_trans_quart_KM$table <- AR_trans_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
AR_trans_quart_KM

jpeg('./results/Missense_variants/AR_Missense_variant_analyses_transmembrane_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_trans_quart_KM)
dev.off()
################################################
## Liver-related
AR_cox_liver <- coxph(Surv(AR_Cox_time, AR_status) ~ Mm_liver + 
                        R_sex + D_sex + R_age + D_age + 
                        Cold_ischemia_time_minutes +
                        Eplets_total_HLAI + Eplets_total_HLAII + 
                        Transplantation_year + Autoimmune_status + 
                        CNI_type_initial_status, data = R_cov_mm_liver)

AR_cox_liver_sum <- summary(AR_cox_liver) %>% coef()
AR_cox_liver_COEF <- bind_cols(AR_cox_liver_sum, 
                                   as.data.frame(exp(confint(AR_cox_liver))))

AR_cox_liver_COEF <- round(AR_cox_liver_COEF, digits = 3)
AR_cox_liver_COEF$`HR(95%_CI)`<- paste0(AR_cox_liver_COEF$`exp(coef)`, "("
                                        ,AR_cox_liver_COEF$`2.5 %`, "-" ,
                                       AR_cox_liver_COEF$`97.5 %`, ")") 
AR_cox_liver_COEF <- rename(AR_cox_liver_COEF, `p_value` = `Pr(>|z|)`,
                                `HR` = `exp(coef)`,
                                `2.5%`= `2.5 %`,
                                `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates")  %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_cox_liver_COEF,
            "./results/Missense_variants/AR_Cox_liver_related_Mm_missense_adjusted_COEF",
            quote = F, row.names = F)

## Analyze the quartiles of the mismatch data in liver-related
## missense variants

# Divide the mismatch sum into quartiles
quantile(R_cov_mm_liver$Mm_liver)
# 0%  25%  50%  75% 100% 
#438  498  518  539  629 

# Create a new column from quartiles (1 = the lowest quartile, 2 = the second, 
# 3 = the third, 4 = the highest quartile)
R_cov_mm_liver$quart_liver <- as.integer(cut(R_cov_mm_liver$Mm_liver, 
                                             quantile(R_cov_mm_liver$Mm_liver, 
                                                      probs = 0:4/4), 
                                             include.lowest = TRUE))

# Cox proportional hazards model for adjusted data
AR_cox_liver_quart <- coxph(Surv(AR_Cox_time, AR_status) ~ quart_liver + R_sex +
                              D_sex + R_age + D_age + 
                              Cold_ischemia_time_minutes + 
                              Eplets_total_HLAI + 
                              Eplets_total_HLAII + Transplantation_year + 
                              Autoimmune_status + 
                              CNI_type_initial_status, data = R_cov_mm_liver)

AR_Cox_liver_quart_sum <- summary(AR_cox_liver_quart) %>% coef()
AR_Cox_liver_quart_COEF <- bind_cols(AR_Cox_liver_quart_sum,
                                     as.data.frame(exp(confint(AR_cox_liver_quart))))

AR_Cox_liver_quart_COEF <- round(AR_Cox_liver_quart_COEF, digits = 3)
AR_Cox_liver_quart_COEF$`HR(95%_CI)` <- paste0(AR_Cox_liver_quart_COEF$`exp(coef)`,
                                           "(",
                                           AR_Cox_liver_quart_COEF$`2.5 %`, "-",
                                          AR_Cox_liver_quart_COEF$`97.5 %`, ")")

AR_Cox_liver_quart_COEF <- rename(AR_Cox_liver_quart_COEF, `p_value` = `Pr(>|z|)`,
                                  `HR` = `exp(coef)`,
                                  `2.5%`= `2.5 %`,
                                  `97.5%` = `97.5 %`) %>% 
  rownames_to_column("covariates")  %>% 
  select(covariates, coef, HR, `se(coef)`, z, `HR(95%_CI)`, p_value)

write.table(AR_Cox_liver_quart_COEF,
            "./results/Missense_variants/AR_Cox_liver_Mm_missense_adjusted_quartiles_COEF_CI",
            quote = F, row.names = F)

# Produce Kapplan-Meier plot with risk table for the manuscript
AR_liver_quart_KM <- ggsurvplot(
  fit = survfit(Surv(AR_Cox_time, AR_status) ~ quart_liver, 
                data = R_cov_mm_liver), 
  xlab = "Time to event (months)", 
  ylab = "Rejection free survival",
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
AR_liver_quart_KM

# The same plot with additional details written in the picture (HR and 95% CI)
AR_liver_quart_KM$plot <- AR_liver_quart_KM$plot+ 
  ggplot2::annotate("text", 
                    x = 160, y = 0.25, # x and y coordinates of the text
                    label = "Adjusted HR 0.985, 95% CI 0.883-1.099 \n P-value 0.791",
                    size = 6)

# Customize risk table title using ggplot verbs
AR_liver_quart_KM$table <- AR_liver_quart_KM$table +
  theme(plot.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16))
AR_liver_quart_KM

jpeg('./results/Missense_variants/AR_Missense_variant_analyses_liver_quartiles_KM.jpeg', 
     width=10, height=15, res=600, units='in')
print(AR_liver_quart_KM)
dev.off()

## Write out the "R_cov_mm_liver" with calculated quartiles for each protein
# group
write.table(R_cov_mm_liver,
            "./data/Missense_variants/R_cov_mm_liver.txt",
            quote = F, row.names = F)

################################################################################
### Perform multiple test corrections for AR for mismatch sum and quartiles
## Mm sums

# Create combined results
AR_results <- list(AR_cox_all_COEF, AR_cox_trans_secr_COEF, 
                   AR_cox_transmemb_COEF, AR_cox_liver_COEF)
names(AR_results) <- c("AR_cox_all_COEF", "AR_cox_trans_secr_COEF", 
                       "AR_cox_transmemb_COEF", "AR_cox_liver_COEF") 

# Extract mismatch covariates and p values
AR_results <- map((AR_results),
                  function(x) {
                    DATA <- x[1,] 
                    DATA <- select(DATA, covariates, p_value)
                    return(DATA)
                  })

AR_results <- bind_rows(AR_results)

# Perform multiple test correction
AR_results$Bonferroni <- p.adjust(AR_results[,2], method = "bonferroni", n = 4)
AR_results$Holm <- p.adjust(AR_results[,2], method = "holm", n = 4)
AR_results$FDR <- p.adjust(AR_results[,2], method = "fdr", n = 4)

write.table(AR_results, 
            "results/Missense_variants/AR_Cox_Mm_mismatch_multiple_testing_results",
            quote = F, row.names = F)

#######################################
## Mm quartiles
# Creating combined results
AR_quart_results <- list(AR_Cox_all_quart_COEF, AR_Cox_t_secr_quart_COEF, 
                         AR_Cox_trans_quart_COEF, AR_Cox_liver_quart_COEF)
names(AR_quart_results) <- c("AR_Cox_all_quart_COEF", "AR_Cox_t_secr_quart_COEF", 
                             "AR_Cox_trans_quart_COEF", "AR_Cox_liver_quart_COEF") 

# Extract Mm covariates and p values
AR_quart_results <- map((AR_quart_results),
                        function(x) {
                          DATA <- x[1,] 
                          DATA <- select(DATA, covariates, p_value)
                          return(DATA)
                        })

AR_quart_results <- bind_rows(AR_quart_results)
AR_quart_results$Bonferroni <- p.adjust(AR_quart_results[,2], 
                                        method = "bonferroni", n = 4)
AR_quart_results$Holm <- p.adjust(AR_quart_results[,2], method = "holm", n = 4)
AR_quart_results$FDR <- p.adjust(AR_quart_results[,2], method = "fdr", n = 4)

write.table(AR_quart_results, 
            "results/Missense_variants/AR_Cox_Mm_quartiles_multiple_testing_results",
            quote = F, row.names = F)
###############################################################################