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
  data <-read_csv(f_name, col_names = TRUE)
  return(data)
}
# fxn to read in multiple AA sequences----
f_file <- function(fast_file){
  #for reading multiple AA sequences from msa package
  fast <- msa::readAAStringSet(fast_file, format = "fasta", use.names = TRUE) # format: biostrings, AAString set
  return(fast)
}

# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file
pcms_AAtree <- function(fast_file){
  #for reading multiple AA sequences from msa package
  #fas <- readAAStringSet(fas_file, format = "fasta", use.names = TRUE)
  #align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package 
  msa::msa(fast, method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  ape::as.AAbin(fast, show.aa = TRUE, check.names = TRUE)
  ape::as.character.AAbin(fast) #converting AAbin data to character format, 
  ape:: as.alignment(fast) # # converting from AAbin to alignment format for phylogenetic tree
  base::as.matrix(fast) # converting alignment to matrix
  # talk to tony about fxn because piping the separate functions together doesn't work, 
  # worried this fxn won't work for the same reason
  return(fast_file)
}

AAbin_labs <- as.matrix(labels(fas_AAbin)) # extraction of the species names
fas_AAbin <- dist.aa(fas_AAbin)
tree <- nj(fas_AAbin)
ggt <-ggtree(tree, cex = 0.8, aes(color=branch.length)) +
  scale_color_continuous(high='lightskyblue1',low='coral4') +
  geom_tiplab(align=TRUE, size=2) +
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt


#Plot the maximum likelihood----
plot_Tree <- function(tree){
plotTree(tree, fsize=0.8,lwd=1,offset=1) #Plot the ML, set font size, create space so nodes aren't on top of each other 
nodelabels(tree$node.label, adj = c(1, 0), frame = "none") #Label the nodes
}

