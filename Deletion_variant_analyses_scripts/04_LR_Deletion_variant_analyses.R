###############################################################################
### Perform regression analyses
### Late rejection
###############################################################################
### General information

# Prerequisites:
# 1. Run script 02_Deletion_variant_mismatches.R
#    data/Deletion_variants/R_dos_pheno_dels_collision.txt
#    data/Deletion_variants/D_dos_pheno_dels_collision.txt

###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(glue)

###############################################################################
## Import matched covariate genomic collision files
R_dos_pheno_dels_collision <- read_table("data/Deletion_variants/R_dos_pheno_dels_collision.txt")
D_dos_pheno_dels_collision <- read_table("data/Deletion_variants/D_dos_pheno_dels_collision.txt")

###############################################################################
## Analyzing the association of mismatch vs non-mismatch to late rejection

LR_analysis <- map(colnames(R_dos_pheno_dels_collision)[103:142], 
                   function(x) {
                     DATA_LR <- select(R_dos_pheno_dels_collision, x, 
                                       Late_rejection_status,
                                       LR_Cox_time,
                                       R_age, D_age, R_sex, D_sex, 
                                       Cold_ischemia_time_minutes, 
                                       Eplets_total_HLAI, Eplets_total_HLAII,
                                       Transplantation_year,
                                       Autoimmune_status,
                                       CNI_type_initial_status)
                     COX_LR <- coxph(Surv(LR_Cox_time, 
                                          Late_rejection_status) ~., 
                                     data = DATA_LR)
                     COEF_LR <- summary(COX_LR) %>% coef()
                     CI_LR <- as.data.frame(exp(confint(COX_LR)))
                     COEF_CI_LR <- bind_cols(COEF_LR, CI_LR)
                     COEF_CI_LR <- COEF_CI_LR %>% rownames_to_column("covariates")
                     write.table(COEF_CI_LR, paste0("./results/Deletion_variants/LR_",x,
                                                    "_adjusted_cox_coef_CI_deletions_mismatch"))
                     return(COEF_CI_LR)
                   })
names(LR_analysis) <- colnames(R_dos_pheno_dels_collision)[103:142]

### Creating a data frame with GL summary statistics for all 40 variants
LR_stats <- map(names(LR_analysis), function(x) { 
  DATA_stat <- LR_analysis[[x]] %>% data.frame()
  DATA_stats <- DATA_stat[1,]
  return(DATA_stats)
})
names(LR_stats) <- names(LR_analysis)

LR_stats_df <- bind_rows(LR_stats) %>% rename("variants" = "covariates")

LR_stats_df[,2:8] <- round(LR_stats_df[,2:8], digits = 3)

LR_stats_df$`HR(95%_CI)` <- paste0(LR_stats_df$exp.coef., "(",
                                   LR_stats_df$X2.5.., "-", LR_stats_df$X97.5..,
                                   ")")

LR_stats_df <- rename(LR_stats_df, "HR" = exp.coef.,
                      "p_value" = Pr...z..,
                      `2.5%` = X2.5..,
                      `97.5%` = X97.5..)

### Control of family-wise errors
LR_stats_df$Bonferroni <- p.adjust(LR_stats_df[,6], method = "bonferroni", 
                                   n = 40)
LR_stats_df$Holm <- p.adjust(LR_stats_df[,6], method = "holm", n = 40)
LR_stats_df$FDR <- p.adjust(LR_stats_df[,6], method = "fdr", n = 40)

write.table(LR_stats_df, 
            "results/Deletion_variants/LR_deletions_mismatch_combined_Cox_and_multiple_testing_results",
            col.names = T, row.names = F)

### Cox plots

## Selecting variants p < 0.05
LR_cox_plot_var <- filter(LR_stats_df, p_value < 0.05)

## Creating new data frame  with two rows, one for each value of collision 
## mismatch; the other covariates are fixed to their average values (if they are 
## continuous variables) or to their lowest level 
## (if they are discrete variables)

LR_cox_plots <- map((LR_cox_plot_var$variants), function(x) {
  DATA <- select(R_dos_pheno_dels_collision, x, 
                 LR_Cox_time, Late_rejection_status, 
                 R_age, D_age, R_sex, D_sex, 
                 Cold_ischemia_time_minutes, 
                 Eplets_total_HLAI, Eplets_total_HLAII,
                 Transplantation_year,
                 Autoimmune_status,
                 CNI_type_initial_status)
  COX <- coxph(Surv(LR_Cox_time, Late_rejection_status) ~., data = DATA)
  NEW_DATA <- with(R_dos_pheno_dels_collision,
                   data.frame(x = c(0,1),
                              R_age = rep(mean(R_age, na.rm = TRUE), 2), 
                              D_age = rep(mean(D_age, na.rm = TRUE), 2), 
                              R_sex = c(1,1), 
                              D_sex = c(1,1),
                              Cold_ischemia_time_minutes = rep(mean(Cold_ischemia_time_minutes,
                                                                    na.rm = TRUE), 2), 
                              Eplets_total_HLAI = rep(mean(Eplets_total_HLAI,
                                                           na.rm = TRUE), 2), 
                              Eplets_total_HLAII = rep(mean(Eplets_total_HLAII,
                                                            na.rm = TRUE), 2),
                              Transplantation_year = rep(mean(Transplantation_year,
                                                              na.rm = TRUE), 2),
                              Autoimmune_status = c(0,0),
                              CNI_type_initial_status = c(0,0)))
  colnames(NEW_DATA)[1] <- x
  FIT <- survfit(COX, newdata = NEW_DATA)
  PLOT <-ggsurvplot(FIT, data = NEW_DATA, 
                    conf.int = F, 
                    legend.labs= c("No mismatch", "Mismatch"),
                    xlab = "Time to event (months)", 
                    ylab = "Rejection-free survival",
                    font.x = c(14, face = "bold"),
                    font.y = c(14, face = "bold"),
                    font.tickslab = c(16),
                    legend.title = "", 
                    font.legend = c(15),
                    legend = c(0.6,0.7))
  
  jpeg(paste0("results/Deletion_variants/LR_",x,"_cox_plot.jpeg"),
       width=10, height=10, res=600, units='in')
  print(PLOT)
  dev.off()   
  # Kaplan Meier plots
  FORM <- as.formula(glue::glue('Surv(LR_Cox_time, Late_rejection_status) ~ {x}'))
  STATS <- LR_cox_plot_var %>% filter(variants == x)
  KM <- survfit(FORM, data = R_dos_pheno_dels_collision)
  KM$call$formula <- FORM
  KM_plot <- ggsurvplot(KM,
                        xlab = "Time to event (months)", 
                        ylab = "Rejection-free survival",
                        font.x = c(16, face = "bold"),
                        font.y = c(16, face = "bold"),
                        font.tickslab = c(16),
                        risk.table = T, fontsize = 5,
                        risk.table.y.text = FALSE,
                        legend.title = "", 
                        font.legend = c(16),
                        legend = c(0.6,0.9))
  KM_plot$plot <- KM_plot$plot +
    ggplot2::annotate("text", 
                      x = 170, y = 0.25,
                      label = paste0("Adjusted HR ", STATS[,3] , ", 95% CI ", 
                                     STATS[,7], "-", STATS[,8],
                                     ", \n P-value " ,STATS[,6]),
                      size = 6)
  KM_plot$table <- KM_plot$table +
    theme(plot.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 16))
  
  jpeg(paste0("results/Deletion_variants/LR_",x,"_KM_plot.jpeg"),
       width=10, height=15, res=600, units='in')
  print(KM_plot)
  dev.off()
})

################################################################################
