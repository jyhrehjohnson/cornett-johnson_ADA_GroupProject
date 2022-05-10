# Phylogenetic Comparative Analysis
#
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

#Create a function that reads in the file information----
read_file <- function(f_name){
  data <-readr::read_csv(f_name, col_names = TRUE)
  return(data)
}
# fxn to read in multiple AA sequences----
f_file <- function(fast_file){
  #for reading multiple AA sequences from msa package
  fast <- Biostrings::readAAStringSet(fast_file, format = "fasta", use.names = TRUE) # format: biostrings, AAString set
  return(fast)
}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file
pcms_AAtree <- function(fast_file){ # commented out until completed because otherwise prevents from knitting
   #for reading multiple AA sequences from msa package
   fas <- readAAStringSet(fas_file, format = "fasta", use.names = TRUE)
  #align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package 
   msa::msa(fast, method = c("Muscle"), type = "protein", order=c("aligned", "input"))
   #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
   ape::as.AAbin(fast, show.aa = TRUE, check.names = TRUE)
     return(fast_file)
}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file
 pcms_AlignAA <- function(file){
   #for reading multiple AA sequences from msa package
   fas <- Biostrings::readAAStringSet(file, format = "fasta", use.names = TRUE)
   # align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package 
   fas_msa <- fas %>% msa::msa(method = c("Muscle"), type = "protein", order=c("aligned", "input"))
   return(fas_msa)
   }

pcms_treeAA <- function(pat){
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  dist_fas_msa <- pat %>% ape::as.AAbin(show.aa = TRUE, check.names = TRUE) %>% 
    ape::dist.aa()
  #neighbor joining method
  tree <- ape::nj(dist_fas_msa)
  ggt <- ggtree::ggtree(tree, cex = 0.8, aes(color=branch.length)) +
    scale_color_continuous(high='green',low='blue') +
    geom_tiplab(align=FALSE, size=2) +
    geom_treescale(y = 0, color = "coral4", fontsize = 4)
  return(ggt)
}


#Plot the maximum likelihood----
plot_Tree <- function(tree){
phytools::plotTree(tree, fsize=0.8,lwd=1,offset=3) #Plot the ML, set font size, create space so nodes aren't on top of each other 
}

