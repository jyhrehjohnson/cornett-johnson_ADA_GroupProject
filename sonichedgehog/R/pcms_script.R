# Phylogenetic Comparative Analysis
# This rscript contains the functions we created for our PCMs
# package.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#Create a function that reads in the file information
read_file <- function(f_name){
  data <-read_csv(f_name, col_names = TRUE)
  return(data)
}

#pcms_tree <- function()
#pcms_plot <- function()