###############################################################################
### Genotype, covariate and dosage data preparations for missense variant
### analyses
###############################################################################

library(tidyverse)
library(data.table)

###############################################################################

### General information

# Prerequisites:
# 1. PLINK 1.90 https://www.cog-genomics.org/plink/
# 2. Folders data, results and src
# 3. Chromosomal genotype data in PLINK format (.bed .bim .fam) in data folder
# 4. Chromosome genotype PLINK format (.bed .bim .fam) with excluded sex chromosomes and MHC in 
#    data folder (EI OLE PREREQUISITE)
# 5. Paired clincal data file
# 6. R_covariates_mm_liver file
# 7. MHC_region.txt file
# 8. Modified "missense_variants"
# 9. Transmembrane and secretory Uniprot protein transcripts
# 10. Transmembrane-only Uniprot protein transcripts
# 11. Liver-spesific Uniprot protein transcripts


### Object "MHC_region.txt" 
# The following script requires a text file "MHC_region.txt" with chromosome,
# starting position and ending position of MHC region.

# 6 28510120 33480577 R1

### Object "missense_variants"
# The following script requires a text file "missense_variants" with missense variants 
# extracted from the chromosomal genotype file using ensemble variant effect
# predictor (VEP)

missense_variants <- fread('data/Missense_variants/missense_variants.txt',
                           data.table = F)
str(missense_variants)
# Structure specification
#'data.frame':	28659 obs. of  7 variables:
# $ #Uploaded_variation: chr  "1_935779_G/A" "1_944029_C/G" "1_946538_G/A" "1_953279_C/T" ...
# $ Location           : chr  "1:935779" "1:944029" "1:946538" "1:953279" ...
# $ Allele             : chr  "A" "G" "A" "T" ...
# $ Feature            : chr  "ENST00000341065" "ENST00000341065" "ENST00000327044" "ENST00000327044" ...
# $ chr                : int  1 1 1 1 1 1 1 1 1 1 ...
# $ position           : int  935779 944029 946538 953279 953858 973289 973858 973929 976215 978953 ...
# $ VEP_variant_ID     : chr  "chr1_935779_G_A" "chr1_944029_C_G" "chr1_946538_G_A" "chr1_953279_C_T" ...
#> str(missense_variants)
#'data.frame':	76838 obs. of  14 variables:
# $ #Uploaded_variation: chr  "1_935779_G/A" "1_935779_G/A" "1_935779_G/A" "1_935779_G/A" ...
# $ Location           : chr  "1:935779" "1:935779" "1:935779" "1:935779" ...
# $ Allele             : chr  "A" "A" "A" "A" ...
# $ Gene               : chr  "ENSG00000187634" "ENSG00000187634" "ENSG00000187634" "ENSG00000187634" ...
# $ Feature            : chr  "ENST00000341065" "ENST00000342066" "ENST00000437963" "ENST00000616016" ...
# $ Feature_type       : chr  "Transcript" "Transcript" "Transcript" "Transcript" ...
# $ Consequence        : chr  "missense_variant" "missense_variant" "missense_variant" "missense_variant" ...
# $ cDNA_position      : int  84 403 373 1359 313 313 262 1359 313 313 ...
# $ CDS_position       : int  85 313 313 850 313 313 262 850 313 313 ...
# $ Protein_position   : int  29 105 105 284 105 105 88 284 105 105 ...
# $ Amino_acids        : chr  "G/S" "G/S" "G/S" "G/S" ...
# $ Codons             : chr  "Ggc/Agc" "Ggc/Agc" "Ggc/Agc" "Ggc/Agc" ...
# $ Existing_variation : chr  "-" "-" "-" "-" ...
# $ Extra              : chr  "IMPACT=MODERATE;STRAND=1;FLAGS=cds_start_NF" "IMPACT=MODERATE;STRAND=1" "IMPACT=MODERATE;STRAND=1;FLAGS=cds_end_NF" #"IMPACT=MODERATE;STRAND=1" ...

### Object "Transmemb_secretory"
# The following script requires a file containing the transmembrane and
# secretory Uniprot protein transcripts

Transmemb_secretory <- read_delim("data/Uniprot_annotation/uniprot-(annotation_(type_transmem)+OR+locations_(location__Se--.tab", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
str(Transmemb_secretory)
# Structure specification
# spec_tbl_df [7,034 × 8] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
# $ Entry             : chr [1:7034] "Q9BQ51" "Q92824" "O43157" "O43508" ...
#$ Ensembl transcript: chr [1:7034] "ENST00000397747;" "ENST00000376752 [Q92824-2];ENST00000545128;" "ENST00000296440;ENST00000358536;ENST00000449094 [O43157-3];ENST00000456774 [O43157-2];" "ENST00000293825;" ...
#$ Entry name        : chr [1:7034] "PD1L2_HUMAN" "PCSK5_HUMAN" "PLXB1_HUMAN" "TNF12_HUMAN" ...
#$ Status            : chr [1:7034] "reviewed" "reviewed" "reviewed" "reviewed" ...
#$ Protein names     : chr [1:7034] "Programmed cell death 1 ligand 2 (PD-1 ligand 2) (PD-L2) (PDCD1 ligand 2) (Programmed death ligand 2) (Butyroph"| __truncated__ "Proprotein convertase subtilisin/kexin type 5 (EC 3.4.21.-) (Proprotein convertase 5) (PC5) (Proprotein convert"| __truncated__ "Plexin-B1 (Semaphorin receptor SEP)" "Tumor necrosis factor ligand superfamily member 12 (APO3 ligand) (TNF-related weak inducer of apoptosis) (TWEAK"| __truncated__ ...
#$ Gene names        : chr [1:7034] "PDCD1LG2 B7DC CD273 PDCD1L2 PDL2" "PCSK5 PC5 PC6" "PLXNB1 KIAA0407 PLXN5 SEP" "TNFSF12 APO3L DR3LG UNQ181/PRO207" ...
#$ Organism          : chr [1:7034] "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" ...
#$ Length            : num [1:7034] 273 1860 2135 249 5412 ...
#- attr(*, "spec")=
#  .. cols(
#    ..   Entry = col_character(),
#    ..   `Ensembl transcript` = col_character(),
#    ..   `Entry name` = col_character(),
#    ..   Status = col_character(),
#    ..   `Protein names` = col_character(),
#    ..   `Gene names` = col_character(),
#    ..   Organism = col_character(),
#    ..   Length = col_double()
#    .. )
#- attr(*, "problems")=<externalptr>

### Object "Transmembrane"
# The following script requires a file containing the transmembrane
# Uniprot protein transcripts

Transmembrane <- read_delim("data/Uniprot_annotation/uniprot-(annotation_(type_transmem)+OR+keyword__Transmembrane+[K--.tab",
                            "\t", escape_double = FALSE, trim_ws = TRUE)

str(Transmembrane)
# Structure specification
#tibble [5,214 ? 8] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#$ Entry             : chr [1:5214] "P22223" "Q9ULX7" "Q9HCJ2" "Q643R3" ...
#$ Ensembl transcript: chr [1:5214] "ENST00000264012;ENST00000429102 [P22223-2];" "ENST00000369111;ENST00000647854;" "ENST00000278198;ENST00000527150;ENST00000528697;ENST00000530763;ENST00000619527;" "ENST00000314891;" ...
#$ Entry name        : chr [1:5214] "CADH3_HUMAN" "CAH14_HUMAN" "LRC4C_HUMAN" "LPCT4_HUMAN" ...
#$ Status            : chr [1:5214] "reviewed" "reviewed" "reviewed" "reviewed" ...
#$ Protein names     : chr [1:5214] "Cadherin-3 (Placental cadherin) (P-cadherin)" "Carbonic anhydrase 14 (EC 4.2.1.1) (Carbonate dehydratase XIV) (Carbonic anhydrase XIV) (CA-XIV)" "Leucine-rich repeat-containing protein 4C (Netrin-G1 ligand) (NGL-1)" "Lysophospholipid acyltransferase LPCAT4 (1-acylglycerol-3-phosphate O-acyltransferase 7) (1-AGP acyltransferase"| __truncated__ ...
#$ Gene names        : chr [1:5214] "CDH3 CDHP" "CA14 UNQ690/PRO1335" "LRRC4C KIAA1580 NGL1 UNQ292/PRO331" "LPCAT4 AGPAT7 AYTL3 LPEAT2" ...
#$ Organism          : chr [1:5214] "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" ...
#$ Length            : num [1:5214] 829 337 640 524 610 114 162 749 605 432 ...

### Object "Liver_spesific"
# The following script requires a tibble "Liver_spesific" containing the 
# liver-spesific Uniprot protein transcripts
Liver_specific <- read_delim("data/Uniprot_annotation/uniprot-annotation_(type__tissue+specificity_+liver)+AND+reviewed_--.tab", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
str(Liver_spesific)
# Structure spesification
#tibble [2,593 ? 8] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
#$ Entry             : chr [1:2593] "Q9ULX7" "Q9UHB6" "O75970" "P53602" ...
#$ Ensembl transcript: chr [1:2593] "ENST00000369111;ENST00000647854;" "ENST00000341247;ENST00000394943 [Q9UHB6-4];ENST00000547825 [Q9UHB6-3];ENST00000552783 [Q9UHB6-5];ENST00000552823 [Q9UHB6-2];" "ENST00000319217;ENST00000447879 [O75970-3];ENST00000536827 [O75970-5];ENST00000541718 [O75970-2];" "ENST00000301012;" ...
#$ Entry name        : chr [1:2593] "CAH14_HUMAN" "LIMA1_HUMAN" "MPDZ_HUMAN" "MVD1_HUMAN" ...
#$ Status            : chr [1:2593] "reviewed" "reviewed" "reviewed" "reviewed" ...
#$ Protein names     : chr [1:2593] "Carbonic anhydrase 14 (EC 4.2.1.1) (Carbonate dehydratase XIV) (Carbonic anhydrase XIV) (CA-XIV)" "LIM domain and actin-binding protein 1 (Epithelial protein lost in neoplasm)" "Multiple PDZ domain protein (Multi-PDZ domain protein 1)" "Diphosphomevalonate decarboxylase (EC 4.1.1.33) (Mevalonate (diphospho)decarboxylase) (MDDase) (Mevalonate pyro"| __truncated__ ...
#$ Gene names        : chr [1:2593] "CA14 UNQ690/PRO1335" "LIMA1 EPLIN SREBP3 PP624" "MPDZ MUPP1" "MVD MPD" ...
#$ Organism          : chr [1:2593] "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" "Homo sapiens (Human)" ...
#$ Length            : num [1:2593] 337 759 2070 400 329 ...


###############################################################################
### Genotype dosage files

## 1) Dosage file for all missense variants (chr 23 and MHC region included)

# Extrac the missense variants from the final data set
system(command = paste0("plink --bfile data/Liver_cohort --extract data/Missense_variants/missense_FINAL.txt --make-bed --out data/Missense_variants/Liver_missense_variants"))
# n= 28659
# Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_variants --recodeA --out data/Missense_variants/LIVER_missense_dosage"))


### 2) Dosage file for all missense variants (chr 23 and MHC region not included)

# Exclude MHC region and using only chromosomes 1-22 
system(command = paste0("plink --bfile data/Liver_cohort --chr 1-22 --exclude data/Missense_variants/MHC_region.txt --range --make-bed --out data/Missense_variants/Liver_Tx_cohort_chr1_22_without_MHC"))

# Extract the missense variants from the filtered data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/missense_FINAL.txt --make-bed --out data/Missense_variants/LIVER_missense_variants_no_X_MHC"))
# n=28225

#Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/LIVER_missense_variants_no_X_MHC --recodeA --out data/Missense_variants/LIVER_missense_dosage_without_X_MHC"))

### 3) Dosage file for transmembrane and secretory proteins

# Extract the transmembrane and secretory missense variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_secr_transm_variants_FINAL.txt --make-bed --out data/Missense_variants/Liver_missense_transmemb_secr"))
# n=9337

# Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_transmemb_secr --recodeA --out data/Missense_variants/Liver_missense_transmemb_secr_dosage"))

# 4) Dosage file for transmembrane-only proteins
# Extract the transmembrane variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_transmemb_variants_FINAL.txt --make-bed --out data/Missense_variants/Liver_missense_transmemb"))
# n=6963

# Create a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_transmemb --recodeA --out data/Missense_variants/Liver_missense_transmembrane_dosage"))

# 5) Dosage file for liver-spesific proteins
# Extract the liver variants from the data set
system(command = paste0("plink --bfile data/Missense_variants/Liver_cohort_wo_MHC --extract data/Missense_variants/List_of_liver_specific_variants_FINAL.txt --make-bed --out data/Missense_variants/Liver_missense_liver_spesific"))
# n=3199

# Creating a dosage file
system(command = paste0("plink --bfile data/Missense_variants/Liver_missense_liver_spesific --recodeA --out data/Missense_variants/liver_missense_liver_spesific_dosage"))
