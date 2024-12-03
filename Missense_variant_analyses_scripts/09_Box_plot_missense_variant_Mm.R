###############################################################################
### Create boxplots that compare grafts that sustained no Rejection/loss events
### vs those that did
###############################################################################
### General information

# Prerequisites:

# 1) Run script 05_AR_missense_variant_analyses.R
#    data/Missense_variants/R_cov_mm_liver.txt
# 2) folder results/Missense_variants

###############################################################################
library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)
library(data.table)
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

###############################################################################
### Turn variables as factors
R_cov_mm_liver$AR_status <- as.factor(R_cov_mm_liver$AR_status)
R_cov_mm_liver$Graft_loss_status <- as.factor(R_cov_mm_liver$Graft_loss_status)
R_cov_mm_liver$Death_status <- as.factor(R_cov_mm_liver$Death_status)
R_cov_mm_liver$Late_rejection_status <- as.factor(R_cov_mm_liver$Late_rejection_status)

### Create boxplots

## AR_status

# Mm_all
AR_all <- ggplot(R_cov_mm_liver, aes(x = AR_status, y = Mm_all, 
                           colour = AR_status)) + geom_boxplot() + 
  guides(color=guide_legend(title="Acute rejection status")) +
  theme_classic() +
  xlab("Acute rejection status") +
  ylab("All protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
AR_all

# Mm trans secr
AR_transm_secr <- ggplot(R_cov_mm_liver, aes(x = AR_status, y = Mm_transm_secr, 
                           colour = AR_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
  guides(color=guide_legend(title="Acute rejection status")) +
  theme_classic() +
  xlab("Acute rejection status") +
  ylab("Transmemb secreted protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
AR_transm_secr

# Mm transm 
AR_transm  <- ggplot(R_cov_mm_liver, aes(x = AR_status, y = Mm_transm , 
                           colour = AR_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#9900FF", "#FFCC33", "#56B4E9")) +
  guides(color=guide_legend(title="Acute rejection status")) +
  theme_classic() +
  xlab("Acute rejection status") +
  ylab("Transmembrane protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
AR_transm 

# Mm liver
AR_liver <- ggplot(R_cov_mm_liver, aes(x = AR_status, y = Mm_liver, 
                           colour = AR_status)) + geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  guides(color=guide_legend(title="Acute rejection status")) +
  theme_classic() +
  xlab("Acute rejection status") +
  ylab("Liver-related protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
AR_liver

AR_box_plots <- ggarrange(AR_all,AR_transm_secr,AR_transm,AR_liver,
                          labels=c("a","b","c","d"),
                          ncol=2,nrow=2)
AR_box_plots

jpeg('./results/Missense_variants/AR_box_plots.jpeg', 
     width=10, height=10, res=600, units='in')
print(AR_box_plots)
dev.off()

############################
## Graft_loss_status

# Mm all
GL_all <- ggplot(R_cov_mm_liver, aes(x = Graft_loss_status, y = Mm_all, 
                           colour = Graft_loss_status)) + geom_boxplot() +
  guides(color=guide_legend(title="Graft loss status")) +
  theme_classic()+
  xlab("Graft loss status") +
  ylab("All protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
GL_all

# Mm transm secr
GL_transm_secr <- ggplot(R_cov_mm_liver, aes(x = Graft_loss_status, 
                                            y = Mm_transm_secr, 
                           colour = Graft_loss_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
  guides(color=guide_legend(title="Graft loss status")) +
  theme_classic() +
  xlab("Graft loss status") +
  ylab("Transmemb secreted protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
GL_transm_secr

# Mm transm 
GL_transm  <- ggplot(R_cov_mm_liver, aes(x = Graft_loss_status, 
                                           y = Mm_transm , 
                           colour = Graft_loss_status)) + geom_boxplot() + 
  guides(color=guide_legend(title="Graft loss status")) +
  scale_color_manual(values=c("#9900FF", "#FFCC33", "#56B4E9")) +
  theme_classic()+
  xlab("Graft loss status") +
  ylab("Transmembrane protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
GL_transm 

# Mm liver
GL_liver <- ggplot(R_cov_mm_liver, aes(x = Graft_loss_status, y = Mm_liver, 
                           colour = Graft_loss_status)) + geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  guides(color=guide_legend(title="Graft loss status")) +
  theme_classic()+
  xlab("Graft loss status") +
  ylab("Liver-related protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
GL_liver

GL_box_plots <- ggarrange(GL_all,GL_transm_secr,GL_transm,GL_liver,
                          labels=c("a","b","c","d"),
                          ncol=2,nrow=2)
GL_box_plots

jpeg('./results/Missense_variants/GL_box_plots.jpeg', 
     width=10, height=10, res=600, units='in')
print(GL_box_plots)
dev.off()

##########################
## Death_status

# Mm all
OS_all <- ggplot(R_cov_mm_liver, aes(x = Death_status, y = Mm_all, 
                           colour = Death_status)) + geom_boxplot() + 
  guides(color=guide_legend(title="Death status")) +
  theme_classic()+
  xlab("Death status") +
  ylab("All protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
OS_all

# Mm transm secr
OS_transm_secr <- ggplot(R_cov_mm_liver, aes(x = Death_status, y = Mm_transm_secr, 
                           colour = Death_status)) + geom_boxplot()  + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
  guides(color=guide_legend(title="Death status")) +
  theme_classic()+
  xlab("Death status") +
  ylab("Transmemb secreted protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
OS_transm_secr

# Mm transm 
OS_transm  <- ggplot(R_cov_mm_liver, aes(x = Death_status, y = Mm_transm , 
                           colour = Death_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#9900FF", "#FFCC33", "#56B4E9")) +
  guides(color=guide_legend(title="Death status")) +
  theme_classic()+
  xlab("Death status") +
  ylab("Transmembrane protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
OS_transm 

# Mm liver
OS_liver <- ggplot(R_cov_mm_liver, aes(x = Death_status, y = Mm_liver, 
                           colour = Death_status)) + geom_boxplot() + 
  scale_color_brewer(palette="Set1") + 
  guides(color=guide_legend(title="Death status")) +
  theme_classic()+
  xlab("Death status") +
  ylab("Liver-related protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
OS_liver

OS_box_plots <- ggarrange(OS_all,OS_transm_secr,OS_transm,OS_liver,
                          labels=c("a","b","c","d"),
                          ncol=2,nrow=2)
OS_box_plots

jpeg('./results/Missense_variants/OS_box_plots.jpeg', 
     width=10, height=10, res=600, units='in')
print(OS_box_plots)
dev.off()

##############################
## Late_rejection_status

# Mm all
LR_all <- ggplot(R_cov_mm_liver, aes(x = Late_rejection_status, y = Mm_all, 
                           colour = Late_rejection_status)) +
  guides(color=guide_legend(title="Late rejection status")) +
  geom_boxplot() +
  theme_classic() +
  xlab("Late rejection status") +
  ylab("All protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
LR_all

# Mm transm secr
LR_transm_secr <- ggplot(R_cov_mm_liver, aes(x = Late_rejection_status, 
                                            y = Mm_transm_secr, 
                           colour = Late_rejection_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
  guides(color=guide_legend(title="Late rejection status")) +
  theme_classic() +
  xlab("Late rejection status") +
  ylab("Transmemb secreted protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
LR_transm_secr

# Mm transm 
LR_transm  <- ggplot(R_cov_mm_liver, aes(x = Late_rejection_status, 
                                           y = Mm_transm , 
                           colour = Late_rejection_status)) + geom_boxplot() + 
  scale_color_manual(values=c("#9900FF", "#FFCC33", "#56B4E9")) +
  guides(color=guide_legend(title="Late rejection status")) +
  theme_classic() +
  xlab("Late rejection status") +
  ylab("Transmembrane protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
LR_transm 

# Mm liver
LR_liver <- ggplot(R_cov_mm_liver, aes(x = Late_rejection_status, y = Mm_liver, 
                           colour = Late_rejection_status)) + geom_boxplot() + 
  scale_color_brewer(palette="Set1") +
  guides(color=guide_legend(title="Late rejection status")) +
  theme_classic() +
  xlab("Late rejection status") +
  ylab("Liver-related protein mismatches") +
  theme(plot.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
LR_liver

LR_box_plots <- ggarrange(LR_all,LR_transm_secr,LR_transm,LR_liver,
                          labels=c("a","b","c","d"),
                          ncol=2,nrow=2)
LR_box_plots

jpeg('./results/Missense_variants/LR_box_plots.jpeg', 
     width=10, height=10, res=600, units='in')
print(LR_box_plots)
dev.off()
