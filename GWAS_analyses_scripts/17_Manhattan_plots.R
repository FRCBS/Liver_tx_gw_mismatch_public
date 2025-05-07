###############################################################################
### Created Manhattan plots from GWAS results
###############################################################################
## General information

# Prerequisites:
# 1) Run script 16_GWAS_analysis.sh
#    results/GWAS/Logistic_regression_recipients.AR_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_recipients.Late_rejection_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_recipients.Graft_loss_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_recipients.Death_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_donors.AR_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_donors.Late_rejection_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_donors.Graft_loss_status.glm.logistic.hybrid
#    results/GWAS/Logistic_regression_donors.Death_status.glm.logistic.hybrid
###############################################################################

library(data.table)
library(tidyverse)
library(qqman)

###############################################################################
## Import the GWAS results for the recipients
AR_recipients <- fread("results/GWAS/Logistic_regression_recipients.AR_status.glm.logistic.hybrid")
LR_recipients <- fread("results/GWAS/Logistic_regression_recipients.Late_rejection_status.glm.logistic.hybrid")
GL_recipients <- fread("results/GWAS/Logistic_regression_recipients.Graft_loss_status.glm.logistic.hybrid")
OS_recipients <- fread("results/GWAS/Logistic_regression_recipients.Death_status.glm.logistic.hybrid")

## Import the GWAS results for the donors
AR_donors <- fread("results/GWAS/Logistic_regression_donors.AR_status.glm.logistic.hybrid")
LR_donors <- fread("results/GWAS/Logistic_regression_donors.Late_rejection_status.glm.logistic.hybrid")
GL_donors <- fread("results/GWAS/Logistic_regression_donors.Graft_loss_status.glm.logistic.hybrid")
OS_donors <- fread("results/GWAS/Logistic_regression_donors.Death_status.glm.logistic.hybrid")

## Filter recipient files by P value
AR_recipients <- filter(AR_recipients, P <= 10^-01)
  # lowest P values: n = 9 of 10^-6
LR_recipients <- filter(LR_recipients, P <= 10^-01)
  # lowest P values: n = 1 of 10^-7
GL_recipients <- filter(GL_recipients, P <= 10^-01)
  # lowest P values: n = 33 of 10^-7
OS_recipients <- filter(OS_recipients, P <= 10^-01)
  # lowest P values: n = 3 of 10^-7

## Filter donor files by P value
AR_donors <- filter(AR_donors, P <= 10^-01)
# lowest P values: n = 1 of 10^-7
LR_donors <- filter(LR_donors, P <= 10^-01)
# lowest P values: n = 3 of 10^-7
GL_donors <- filter(GL_donors, P <= 10^-01)
# lowest P values: n = 10 of 10^-7
OS_donors <- filter(OS_donors, P <= 10^-01)
# lowest P values: n = 56 of 10^-6

#######################################
## Edit data frames

# Recipients
AR_recipients$`#CHROM`[AR_recipients$`#CHROM` == "X"] <-23
LR_recipients$`#CHROM`[LR_recipients$`#CHROM` == "X"] <-23
GL_recipients$`#CHROM`[GL_recipients$`#CHROM` == "X"] <-23
OS_recipients$`#CHROM`[OS_recipients$`#CHROM` == "X"] <-23

AR_recipients$`#CHROM` <- as.numeric(AR_recipients$`#CHROM`)
LR_recipients$`#CHROM` <- as.numeric(LR_recipients$`#CHROM`)
GL_recipients$`#CHROM` <- as.numeric(GL_recipients$`#CHROM`)
OS_recipients$`#CHROM` <- as.numeric(OS_recipients$`#CHROM`)

# Donors
AR_donors$`#CHROM`[AR_donors$`#CHROM` == "X"] <-23
LR_donors$`#CHROM`[LR_donors$`#CHROM` == "X"] <-23
GL_donors$`#CHROM`[GL_donors$`#CHROM` == "X"] <-23
OS_donors$`#CHROM`[OS_donors$`#CHROM` == "X"] <-23

AR_donors$`#CHROM` <- as.numeric(AR_donors$`#CHROM`)
LR_donors$`#CHROM` <- as.numeric(LR_donors$`#CHROM`)
GL_donors$`#CHROM` <- as.numeric(GL_donors$`#CHROM`)
OS_donors$`#CHROM` <- as.numeric(OS_donors$`#CHROM`)

## Create Manhattan plots of the results

# Recipients
jpeg('./results/GWAS/AR_recipients_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(AR_recipients, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey30", "orange2")))
dev.off()

jpeg('./results/GWAS/LR_recipients_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(LR_recipients, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey30", "orange2")))
dev.off()

jpeg('./results/GWAS/GL_recipients_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(GL_recipients, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey30", "orange2")))
dev.off()

jpeg('./results/GWAS/OS_recipients_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(OS_recipients, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey30", "orange2")))
dev.off()

# Donors
jpeg('./results/GWAS/AR_donors_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(AR_donors, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey33", "mediumpurple1")))
dev.off()

jpeg('./results/GWAS/LR_donors_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(LR_donors, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey33", "mediumpurple1")))
dev.off()

jpeg('./results/GWAS/GL_donors_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(GL_donors, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey33", "mediumpurple1")))
dev.off()

jpeg('./results/GWAS/OS_donors_Manhattan.jpeg', 
     width=15, height=8, res=300, units='in')
print(manhattan(OS_donors, chr="#CHROM", bp = "POS", snp = "ID", p ="P",
                ylim = c(1, 8), cex = 0.8, cex.axis = 0.8, cex.lab = 1.3,
                col = c("grey33", "mediumpurple1")))
dev.off()
