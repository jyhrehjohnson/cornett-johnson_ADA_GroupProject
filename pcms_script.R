# Phylogenetic Comparative Analysis
#
# This rscript contains the functions we created for our PCMs
# package.

require(tidyverse)

f <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_orthologs.csv"
d <- read_csv(f, col_names = TRUE)

#pcms_tree <- function()
#pcms_plot <- function()