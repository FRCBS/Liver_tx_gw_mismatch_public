###############################################################################
### PCA plot creation
###############################################################################
## General information

# Prerequisites:
# 1) Run script 13_LD_pruning_and_PCA_analysis.sh
#    data/PCA/Post_PCA.eigenval
###############################################################################
library(data.table)
library(tidyverse)

###############################################################################
## Upload Post_PCA.eigenval file
eigenval <- fread("data/PCA/Post_PCA.eigenval")

eigenval$PCA_number <- row_number(eigenval) %>% order(decreasing = T)

eigenval_plot <- plot(eigenval$PCA_number, eigenval$V1, col = "red")
axis(1, at = seq(1, 10, by = 1), las=2)


