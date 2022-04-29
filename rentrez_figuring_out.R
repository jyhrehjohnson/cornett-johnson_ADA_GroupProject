require(tidyverse)
require(rentrez)
require(devtools)
require(usethis)
require(devtools)
require(roxygen2)
require(withr)

f <- "https://raw.githubusercontent.com/jyhrehjohnson/cornett-johnson_ADA_GroupProject/main/SHH_orthologs.csv"
d <- read_csv(f, col_names = TRUE)

SHH_search <- entrez_search(db="gene",
                            term ="(SHH[All Fields] AND (Gnathostomata [Organism] OR Gnathostomata [Organism] OR Gnathostomata[All Fields]) AND orthologs[All Fields]) AND alive[prop]", 
                            use_history=TRUE)
SHH_search$web_history 

SHH_aa.seq <- entrez_fetch(db = "sequences", rettype = "fasta", web_history = SHH_search$web_history)

sequences <- 
d <- d %>% mutate(entrez_fetch())