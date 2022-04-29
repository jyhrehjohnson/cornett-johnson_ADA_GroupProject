require(tidyverse)
require(rentrez)
require(devtools)
require(usethis)
require(devtools)
require(roxygen2)
require(withr)

f <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_orthologs_1.csv"
d <- read_csv(f, col_names = TRUE)
SHH_aa.seq <- entrez_fetch(db = "protein", 
                           id = rowwise(d$RefSeq_Protein_accessions), 
                           rettype = "fasta", 
                          # web_history = SHH_search$web_history
                           )
SHH_aa.seq
sequences <- 
d <- d %>% mutate(SHH_Peptides = SHH_aa.seq)
