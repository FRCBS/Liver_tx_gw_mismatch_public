################################################################################
### Calculation of deletion mismatches and generation of covariate files
################################################################################
### General information

# Prerequisites:
# 1. Run script 01_Deletion_variant_dosages.R
#    data/Deletion_variants/Liver_dels_variants.raw
# 2.data/Your_clinical_data.txt

###############################################################################

library(tidyverse)
library(data.table)
library(janitor)
library(survminer)
library(Amelia)

################################################################################
# Import the clinical data file 
Covariates <- read_table("data/Your_clinical_data")

# Structure specification example
str(Covariates)
# spec_tbl_df [666 × 21] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ Pair_number               : num [1:666] 539 540 542 543 541 544 545 546 547 548 ...
# $ R_pseudo                  : chr [1:666] "R_pseudo1" "R_pseudo2" "R_pseudo3" "R_pseudo4" ...
# $ D_pseudo                  : chr [1:666] "D_pseudo1" "D_pseudo2" "D_pseudo3" "D_pseudo4" ...
# $ R_sex                     : num [1:666] 1 1 2 1 2 1 1 1 1 2 ...
# $ R_age                     : num [1:666] 47 19 49 62 31 52 48 36 30 71 ...
# $ D_age                     : num [1:666] 55 31 73 57 20 57 50 64 49 64 ...
# $ D_sex                     : num [1:666] 1 1 2 2 1 2 2 2 1 2 ...
# $ Cold_ischemia_time_minutes: num [1:666] 317 307 723 449 273 244 353 246 270 231 ...
# $ AR_status                 : num [1:666] 0 1 0 0 1 0 1 1 0 0 ...
# $ Eplets_total_HLAI         : num [1:666] 10 19 19 19 10 11 11 29 13 26 ...
# $ Eplets_total_HLAII        : num [1:666] 4 7 8 4 12 15 9 23 10 11 ...
# $ AR_Cox_time               : num [1:666] 86.63 0.197 83.605 83.605 0.723 ...
# $ Transplantation_year      : num [1:666] 2013 2013 2014 2014 2014 ...
# $ Autoimmune_status         : num [1:666] 1 1 0 0 1 0 1 0 1 0 ...
# $ CNI_type_initial_status   : num [1:666] 1 1 1 1 1 1 1 1 1 1 ...
# $ LR_Cox_time               : num [1:666] 86.6 85 83.6 83.6 13.5 ...
# $ Late_rejection_status     : num [1:666] 0 0 0 0 1 0 0 0 0 0 ...
# $ Graft_loss_status         : num [1:666] 0 0 0 0 0 0 1 0 0 0 ...
# $ graft_loss_months         : num [1:666] 86.6 85 83.6 83.6 83.6 ...
# $ Death_status              : num [1:666] 0 0 0 0 0 0 0 0 0 0 ...
# $ Death_Cox_time_months     : num [1:666] 86.6 85 83.6 83.6 83.6 ...

# Import dosage file
Liver_dels_dosage <- read_table("data/Deletion_variants/Liver_dels_dosage.raw")

# Structure of Liver_dels_dosage
# spec_tbl_df [1,567 × 46] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#  $ FID                  : chr [1:1567] "R_pseudo1" "R_pseudo2" "R_pseudo3" "R_pseudo4" ...
#  $ IID                  : chr [1:1567] "R_pseudo1" "R_pseudo2" "R_pseudo3" "R_pseudo4" ...
#  $ PAT                  : num [1:1567] 0 0 0 0 0 0 0 0 0 0 ...
#  $ MAT                  : num [1:1567] 0 0 0 0 0 0 0 0 0 0 ...
#  $ SEX                  : num [1:1567] 0 0 0 0 0 0 0 0 0 0 ...
#  $ PHENOTYPE            : num [1:1567] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
#  $ chr1_15839166_C_G_G  : num [1:1567] 0 1 1 2 0 1 1 1 0 0 ...
#  $ chr1_25426560_T_C_T  : num [1:1567] 1 1 1 1 0 1 1 0 2 2 ...
#  $ chr1_72346221_G_T_G  : num [1:1567] 1 2 0 1 1 1 0 1 1 0 ...
#  $ chr1_152618187_T_C_T : num [1:1567] 0 0 0 2 0 1 1 0 0 0 ...
#  $ chr1_152786728_A_C_C : num [1:1567] 0 0 0 0 0 0 0 1 0 1 ...
#  $ chr1_196854483_A_G_G : num [1:1567] 0 2 1 0 0 1 1 0 0 1 ...
#  $ chr1_222193326_G_C_C : num [1:1567] 0 1 0 0 0 0 2 1 2 2 ...
#  $ chr2_108606280_G_A_G : num [1:1567] 0 0 1 1 0 0 0 1 0 1 ...
#  $ ...

###############################################################################
# Join the dosage files with the covariate file
R_dos_pheno_dels <- inner_join(Covariates, Liver_dels_dosage,
                               by = c("R_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE) 

D_dos_pheno_dels <- inner_join(Covariates, Liver_dels_dosage,
                               by = c("D_pseudo" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# FOR RECIPIENTS

# Recode the SNPs accordingly:
# 0 = homozygous for major allele or heterozygous (NON-RISK ALLELE)
# 1 = homozygous for minor allele (RISK ALLELE)

# NOTION! Minor alleles are the ones that are tagging the SNPs (usually) and
# it needs to be checked in Your_genotype_data!

# In the dosage file, the homozygosity/heterozygosity is coded with 0, 1 or 2
# 2 means that there 2 alleles for the specific allele for which the calculation
# was done
# Recode 2 (homozygous for minor allele) into value 1 and value 1 into 0

R_dos_pheno_dels$rs10927864_0_1 <- recode(R_dos_pheno_dels$chr1_15839166_C_G_G , '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11209948_0_1 <- recode(R_dos_pheno_dels$chr1_72346221_G_T_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11249248_0_1 <- recode(R_dos_pheno_dels$chr1_25426560_T_C_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11587012_0_1 <- recode(R_dos_pheno_dels$chr1_152786728_A_C_C, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs158736_0_1 <- recode(R_dos_pheno_dels$chr1_222193326_G_C_C, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs6693105_0_1 <- recode(R_dos_pheno_dels$chr1_152618187_T_C_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs7542235_0_1 <- recode(R_dos_pheno_dels$chr1_196854483_A_G_G, '2' = 1, '1' = 0) # CFH
R_dos_pheno_dels$rs7419565_0_1 <- recode(R_dos_pheno_dels$chr2_130195449_T_C_C, '2' = 1, '1' = 0) 
R_dos_pheno_dels$rs893403_0_1 <- recode(R_dos_pheno_dels$chr2_108606280_G_A_G, '2' = 1, '1' = 0) # LIMS1
R_dos_pheno_dels$rs10053292_0_1 <- recode(R_dos_pheno_dels$chr5_149896561_T_C_C, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs2387715_0_1 <- recode(R_dos_pheno_dels$chr5_180934266_A_T_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs7703761_0_1 <- recode(R_dos_pheno_dels$chr5_148177325_C_T_C, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs17654108_0_1 <- recode(R_dos_pheno_dels$chr6_121480387_A_T_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs2160195_0_1 <- recode(R_dos_pheno_dels$chr7_38365370_A_T_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs4621754_0_1 <- recode(R_dos_pheno_dels$chr7_158706303_A_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs4729606_0_1 <- recode(R_dos_pheno_dels$chr7_100724167_T_C_C, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11985201_0_1 <- recode(R_dos_pheno_dels$chr8_39581402_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs4543566_0_1 <- recode(R_dos_pheno_dels$chr8_6968014_C_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs1523688_0_1 <- recode(R_dos_pheno_dels$chr9_104592374_T_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs2174926_0_1 <- recode(R_dos_pheno_dels$chr9_118763272_A_G_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs10885336_0_1 <- recode(R_dos_pheno_dels$chr10_112351444_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs2342606_0_1 <- recode(R_dos_pheno_dels$chr10_80034407_T_C_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs3793917_0_1 <- recode(R_dos_pheno_dels$chr10_122459759_C_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11228868_0_1 <- recode(R_dos_pheno_dels$chr11_55252994_C_T_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs1944862_0_1 <- recode(R_dos_pheno_dels$chr11_55521460_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs1478309_0_1 <- recode(R_dos_pheno_dels$chr12_10389146_T_G_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs9318648_0_1 <- recode(R_dos_pheno_dels$chr13_24566498_A_G_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs11156875_0_1 <- recode(R_dos_pheno_dels$chr14_35150610_A_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs8007442_0_1 <- recode(R_dos_pheno_dels$chr14_21924807_T_C_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs8022070_0_1 <- recode(R_dos_pheno_dels$chr14_81411038_C_T_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs10521145_0_1 <- recode(R_dos_pheno_dels$chr16_28585563_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs2244613_0_1 <- recode(R_dos_pheno_dels$chr16_55810697_G_T_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs16966699_0_1 <- recode(R_dos_pheno_dels$chr17_41343914_C_G_G, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs8064493_0_1 <- recode(R_dos_pheno_dels$chr17_41237071_A_G_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs103294_0_1 <- recode(R_dos_pheno_dels$chr19_54293995_T_C_T, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs324121_0_1 <- recode(R_dos_pheno_dels$chr19_52391192_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs3810336_0_1 <- recode(R_dos_pheno_dels$chr19_56175525_G_A_A, '2' = 1, '1' = 0)
R_dos_pheno_dels$rs4806152_0_1 <- recode(R_dos_pheno_dels$chr19_35395758_A_C_C, '2' = 1, '1' = 0)

# Minor allele different than in the original publication
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$chr7_101357735_A_G_G, '2' = 1)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474_0_1, '0' = 2)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474_0_1, '1' = 0)
R_dos_pheno_dels$rs6943474_0_1 <- recode(R_dos_pheno_dels$rs6943474_0_1, '2' = 1)


R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$chr11_48548630_A_G_A, '2' = 1)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017_0_1, '0' = 2)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017_0_1, '1' = 0)
R_dos_pheno_dels$rs4882017_0_1 <- recode(R_dos_pheno_dels$rs4882017_0_1, '2' = 1)

# FOR DONORS

# Recode the SNPs accordingly:
# 1 = homozygous for major allele or heterozygous
# 2 = homozygous for minor allele

# Recode 0 (homozygous for major allele) into value 1
# (the other values remain the same (1=1 and 2=2))
# for every SNP except for
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017)
D_dos_pheno_dels$rs10927864_0_1 <- recode(D_dos_pheno_dels$chr1_15839166_C_G_G , '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11209948_0_1 <- recode(D_dos_pheno_dels$chr1_72346221_G_T_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11249248_0_1 <- recode(D_dos_pheno_dels$chr1_25426560_T_C_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11587012_0_1 <- recode(D_dos_pheno_dels$chr1_152786728_A_C_C, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs158736_0_1 <- recode(D_dos_pheno_dels$chr1_222193326_G_C_C, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs6693105_0_1 <- recode(D_dos_pheno_dels$chr1_152618187_T_C_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs7542235_0_1 <- recode(D_dos_pheno_dels$chr1_196854483_A_G_G, '2' = 1, '1' = 0) # CFH
D_dos_pheno_dels$rs7419565_0_1 <- recode(D_dos_pheno_dels$chr2_130195449_T_C_C, '2' = 1, '1' = 0) 
D_dos_pheno_dels$rs893403_0_1 <- recode(D_dos_pheno_dels$chr2_108606280_G_A_G, '2' = 1, '1' = 0) # LIMS1
D_dos_pheno_dels$rs10053292_0_1 <- recode(D_dos_pheno_dels$chr5_149896561_T_C_C, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs2387715_0_1 <- recode(D_dos_pheno_dels$chr5_180934266_A_T_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs7703761_0_1 <- recode(D_dos_pheno_dels$chr5_148177325_C_T_C, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs17654108_0_1 <- recode(D_dos_pheno_dels$chr6_121480387_A_T_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs2160195_0_1 <- recode(D_dos_pheno_dels$chr7_38365370_A_T_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs4621754_0_1 <- recode(D_dos_pheno_dels$chr7_158706303_A_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs4729606_0_1 <- recode(D_dos_pheno_dels$chr7_100724167_T_C_C, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11985201_0_1 <- recode(D_dos_pheno_dels$chr8_39581402_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs4543566_0_1 <- recode(D_dos_pheno_dels$chr8_6968014_C_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs1523688_0_1 <- recode(D_dos_pheno_dels$chr9_104592374_T_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs2174926_0_1 <- recode(D_dos_pheno_dels$chr9_118763272_A_G_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs10885336_0_1 <- recode(D_dos_pheno_dels$chr10_112351444_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs2342606_0_1 <- recode(D_dos_pheno_dels$chr10_80034407_T_C_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs3793917_0_1 <- recode(D_dos_pheno_dels$chr10_122459759_C_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11228868_0_1 <- recode(D_dos_pheno_dels$chr11_55252994_C_T_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs1944862_0_1 <- recode(D_dos_pheno_dels$chr11_55521460_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs1478309_0_1 <- recode(D_dos_pheno_dels$chr12_10389146_T_G_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs9318648_0_1 <- recode(D_dos_pheno_dels$chr13_24566498_A_G_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs11156875_0_1 <- recode(D_dos_pheno_dels$chr14_35150610_A_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs8007442_0_1 <- recode(D_dos_pheno_dels$chr14_21924807_T_C_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs8022070_0_1 <- recode(D_dos_pheno_dels$chr14_81411038_C_T_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs10521145_0_1 <- recode(D_dos_pheno_dels$chr16_28585563_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs2244613_0_1 <- recode(D_dos_pheno_dels$chr16_55810697_G_T_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs16966699_0_1 <- recode(D_dos_pheno_dels$chr17_41343914_C_G_G, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs8064493_0_1 <- recode(D_dos_pheno_dels$chr17_41237071_A_G_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs103294_0_1 <- recode(D_dos_pheno_dels$chr19_54293995_T_C_T, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs324121_0_1 <- recode(D_dos_pheno_dels$chr19_52391192_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs3810336_0_1 <- recode(D_dos_pheno_dels$chr19_56175525_G_A_A, '2' = 1, '1' = 0)
D_dos_pheno_dels$rs4806152_0_1 <- recode(D_dos_pheno_dels$chr19_35395758_A_C_C, '2' = 1, '1' = 0)

# Minor allele different than in the original publication
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$chr7_101357735_A_G_G, '2' = 1)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474_0_1, '0' = 2)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474_0_1, '1' = 0)
D_dos_pheno_dels$rs6943474_0_1 <- recode(D_dos_pheno_dels$rs6943474_0_1, '2' = 1)

D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$chr11_48548630_A_G_A, '2' = 1)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017_0_1, '0' = 2)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017_0_1, '1' = 0)
D_dos_pheno_dels$rs4882017_0_1 <- recode(D_dos_pheno_dels$rs4882017_0_1, '2' = 1)

# Remove extra columns from both 'D_dos_pheno_dels' and 'R_dos_pheno_dels'
# covariate files (so that we only have the 'pairs' column and the 
# deletion-columns with either value 1 or 0)

R_paired_dosage_collision  <- select(R_dos_pheno_dels,
                                     "Pair_number", 62:101)
missmap(R_paired_dosage_collision)

D_paired_dosage_collision <- select(D_dos_pheno_dels,
                                    "Pair_number", 62:101)
missmap(D_paired_dosage_collision)

# Calculate the mismatch sum
# The function to calculate the sum is:
Mismatch <- function(R,D) {
  ifelse((D==0 & R==1), 1, 0) 
}

Collision <- sapply(2:ncol(R_paired_dosage_collision), 
                    function(i) {Mismatch(R_paired_dosage_collision[,i], 
                                          D_paired_dosage_collision[,i])})

Collision_df <- data.frame(Pair_number=R_paired_dosage_collision$Pair_number, 
                           Collision)

# Rename the columns in genomic collision data frame
# NOTION! The renaming requires accuracy and it needs to be made sure that
# each column actually represents the right variant!
Collision_df <- rename(Collision_df, rs10927864_col = X1)
Collision_df <- rename(Collision_df, rs11209948_col = X2)
Collision_df <- rename(Collision_df, rs11249248_col = X3)
Collision_df <- rename(Collision_df, rs11587012_col = X4)
Collision_df <- rename(Collision_df, rs158736_col = X5)
Collision_df <- rename(Collision_df, rs6693105_col = X6)
Collision_df <- rename(Collision_df, rs7542235_col = X7)
Collision_df <- rename(Collision_df, rs7419565_col = X8)
Collision_df <- rename(Collision_df, rs893403_col = X9)
Collision_df <- rename(Collision_df, rs10053292_col = X10)
Collision_df <- rename(Collision_df, rs2387715_col = X11)
Collision_df <- rename(Collision_df, rs7703761_col = X12)
Collision_df <- rename(Collision_df, rs17654108_col = X13)
Collision_df <- rename(Collision_df, rs2160195_col = X14)
Collision_df <- rename(Collision_df, rs4621754_col = X15)
Collision_df <- rename(Collision_df, rs4729606_col = X16)
Collision_df <- rename(Collision_df, rs11985201_col = X17)
Collision_df <- rename(Collision_df, rs4543566_col = X18)
Collision_df <- rename(Collision_df, rs1523688_col = X19)
Collision_df <- rename(Collision_df, rs2174926_col = X20)
Collision_df <- rename(Collision_df, rs10885336_col = X21)
Collision_df <- rename(Collision_df, rs2342606_col = X22)
Collision_df <- rename(Collision_df, rs3793917_col = X23)
Collision_df <- rename(Collision_df, rs11228868_col = X24)
Collision_df <- rename(Collision_df, rs1944862_col = X25)
Collision_df <- rename(Collision_df, rs1478309_col = X26)
Collision_df <- rename(Collision_df, rs9318648_col = X27)
Collision_df <- rename(Collision_df, rs11156875_col = X28)
Collision_df <- rename(Collision_df, rs8007442_col = X29)
Collision_df <- rename(Collision_df, rs8022070_col = X30)
Collision_df <- rename(Collision_df, rs10521145_col = X31)
Collision_df <- rename(Collision_df, rs2244613_col = X32)
Collision_df <- rename(Collision_df, rs16966699_col = X33)
Collision_df <- rename(Collision_df, rs8064493_col = X34)
Collision_df <- rename(Collision_df, rs103294_col = X35)
Collision_df <- rename(Collision_df, rs324121_col = X36)
Collision_df <- rename(Collision_df, rs3810336_col = X37)
Collision_df <- rename(Collision_df, rs4806152_col = X38)
Collision_df <- rename(Collision_df, rs6943474_col = X39)
Collision_df <- rename(Collision_df, rs4882017_col = X40)

# Match the covariate files with genomic collision result
# (for both recipients and donors)
R_dos_pheno_dels_collision <- inner_join(R_dos_pheno_dels,
                                         Collision_df, by = "Pair_number")
D_dos_pheno_dels_collision <- inner_join(D_dos_pheno_dels,
                                         Collision_df, by = "Pair_number")

###############################################################################
# Write out the files
write.table(R_dos_pheno_dels_collision,
             file = "data/Deletion_variants/R_dos_pheno_dels_collision.txt",
             row.names = F, col.names = T, quote = F)

write.table(D_dos_pheno_dels_collision,
            file = "data/Deletion_variants/D_dos_pheno_dels_collision.txt",
            row.names = F, col.names = T, quote = F)
###############################################################################
