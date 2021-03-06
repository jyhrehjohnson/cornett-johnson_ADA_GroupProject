require(tidyverse)
# require(devtools)
#require(usethis)
#require(devtools)
#require(roxygen2)
# require(withr)
require(ape)
require(phangorn)
require(phytools) 
require(BiocManager) # to get packages from bioconductor (ie package not on CRAN)
require(ggtree) # extension of ggplot to make phylo trees
require(msa) # multiple sequence alignment
require(Biostrings) # goes with the msa package
require(ggplot2)


f <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_orthologs_Carnivora.csv"
d <- read_csv(f, col_names = TRUE)

#By Hand ----
# read in sequences from .fasta file
fas <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_Carnivora.fasta"
#for reading multiple AA sequences from msa package
file <- readAAStringSet(fas, format = "fasta", use.names = TRUE) # format: biostrings, AAString set

#align the fasta file using MUSCLE algorithm
fas_msa <-  msa(file, method = c("Muscle"), type = "protein", order=c("aligned", "input")) # multiple sequence alignment from msa package 
fas_AAbin <- as.AAbin(fas_msa, show.aa = TRUE, check.names = TRUE) # #read aligned data, converting to AAbin
fas_AAch <- as.character.AAbin(fas_AAbin) #converting AAbin to character strings
fas_align <- as.alignment(fas_AAch)
fas_AAmatrix <- as.matrix(fas_align) # converting alignment to matrix
# converting from AAbin to alignment format
AAbin_labs <- as.matrix(labels(fas_AAbin)) # extraction of the species names
fas_AAbin <- dist.aa(fas_AAbin)
tree <- nj(fas_AAbin)

# plot the NJ tree ---- 
ggt <-ggtree(tree, cex = 0.8, aes(color = branch.length)) +
  scale_color_continuous(high='green',low='blue') +
  geom_tiplab(align = TRUE, size = 4) +
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt


# building a function to make a polypeptide (Amino Acid) sequence phylogenetic tree from data in fasta file
pcms_AAtree <- function(fast_file){ # commented out until completed because otherwise prevents from knitting
  #for reading multiple AA sequences from msa package
  fas <- readAAStringSet(fas_file, format = "fasta", use.names = TRUE)
  #align the fasta file using MUSCLE algorithm: multiple sequence alignment from msa package 
  msa::msa(fas, method = c("Muscle"), type = "protein", order=c("aligned", "input"))
  #read aligned data, storing in AAbin format (class will be AAbin) (ape package)
  ape::as.AAbin(fast, show.aa = TRUE, check.names = TRUE)
  return(fast_file)
}


#Plot the maximum likelihood ---
plotTree(tree, fsize=0.8,lwd=1,offset=1) #Plot the ML, set font size, create space so nodes aren't on top of each other 
nodelabels(tree$node.label, adj = c(1, 0), frame = "none") #label the nodes


# %>%  version fucken works bby! ----
fas_msa <- msa(file, method = c("Muscle"), type = "protein", order=c("aligned", "input")) # multiple sequence alignment from msa package 
dist_fas_msa <- fas_msa %>% as.AAbin(show.aa = TRUE, check.names = TRUE) %>% #read aligned data, converting to AAbin
  dist.aa() # branch length calculations

fas_msa <- fas_msa %>% as.AAbin(show.aa = TRUE, check.names = TRUE) %>% 
  as.character.AAbin() %>% #converting AAbin to character strings
  as.alignment() %>% 
  as.matrix()
 

tree <- nj(dist_fas_msa)
tree # gives tree details
plot(tree) # gives a hard to read version of the tree (nice!)


# somehow at this step, i'm getting a lot of data loss (ie only one sequence still has AA in it when i get to the matrix)
fas_align <- as.alignment(fas_AAch)
fas_AAmatrix <- as.matrix(fas_align) # converting alignment to matrix
# converting from AAbin to alignment format
AAbin_labs <- as.matrix(labels(fas_AAbin)) # extraction of the species names
fas_AAbin <- dist.aa(fas_AAbin)
tree <- nj(fas_AAbin)