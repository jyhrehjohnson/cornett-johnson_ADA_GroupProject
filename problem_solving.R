require(tidyverse)
require(rentrez)
require(devtools)
require(usethis)
require(devtools)
require(roxygen2)
require(withr)
require(ape)

f <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_orthologs_Carnivora.csv"
d <- read_csv(f, col_names = TRUE)
# read in sequences from .fasta file
fas <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_Carnivora.fasta"
#for reading multiple AA sequences from msa package
file <- readAAStringSet(fas, format = "fasta", use.names = TRUE)

#align the fasta file using MUSCLE algorithm
fas_msa <-  msa(file, method = c("Muscle"), type = "protein", order=c("aligned", "input")) # multiple sequence alignment from msa package 
fas_AAbin <- as.AAbin(fas_msa, show.aa = TRUE, check.names = TRUE) # #read aligned data, converting to AAbin
fas_AAch <- as.character.AAbin(fas_AAbin) #converting AAbin to character strings

# somehow at this step, i'm getting a lot of data loss (ie only one sequence still has AA in it when i get to the matrix)
fas_align <- as.alignment(fas_AAch)
fas_AAmatrix <- as.matrix(fas_AAch) # converting alignment to matrix
 # converting from AAbin to alignment format
AAbin_labs <- as.matrix(labels(fas_AAbin)) # extraction of the species names


## an<-as.alignment(nbin)  #converting DNAbin to alignment format
## nm<-as.matrix(an)       #converting alignment to matrix
# nbinmat<-as.matrix(labels(nbin)) #extraction of the sample names
# nbin
# class(nbin)
# dnbin<-dist.dna(nbin, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
# tree<-nj(dnbin)