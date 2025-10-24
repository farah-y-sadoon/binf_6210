##### Software Tools - Class 14 - Challenge
##***************************
##  Software Tools - Class 14 - Challenge
##
## Farah Sadoon
##
## 2025-10-23
##
##***************************

## _ Packages used -------
library(tidyverse)
library(seqinr)
library(rentrez)
library(Biostrings)

# PART 6: CHALLENGE! ----

# 1. Decide as a group a taxon you would like to search the "nuccore" database for. Each group member can then search for sequence records for a different gene for that taxon. Share with your group information about your search terms and the number of hits you received for that gene.

mammal_search <- entrez_search(db = "nuccore", term = "Mammalia[ORGN] AND COX3[Gene]", retmax = 100)

# 2. Try to get the data for your gene into dataframe format like we did in part 4. The dataframe should contain columns for the original title of the sequences (e.g. CytB title), the sequence itself (e.g. CytB sequence), and species name.

cox3_fetch <- entrez_fetch(db = "nuccore", id = mammal_search$ids, rettype = "fasta")
write(cox3_fetch, "cox3_fetch.fasta", sep = "\n")
string_set <- readDNAStringSet("cox3_fetch.fasta")
df_cox3 <- data.frame(cox3_title = names(string_set), cox3_sequence = paste(string_set))
df_cox3$species_name <- word(df_cox3$cox3_title, 2L, 3L)
df_cox3 <- df_cox3[, c("cox3_title", "species_name", "cox3_sequence")]
view(df_cox3)

# 3. As a group, decide on a feature of your sequence data that you are interested in. For example, the sequence lengths or the "AT" proportions of your sequences. Create a new column in your dataframe that holds this information. Then, compare with your group the average value of that feature for your gene.

# calculate sequence length
df_cox3$sequence_length <- nchar(df_cox3$cox3_sequence)
avg_seq_length <- mean(df_cox3$sequence_length)




