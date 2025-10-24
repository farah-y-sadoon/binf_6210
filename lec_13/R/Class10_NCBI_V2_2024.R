##### Software Tools - Class 10 - NCBI V2 ----
##***************************
##  Software Tools - Class 10 - NCBI V2
##
## Karl Cottenie
##
## 2024-10-03
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(Biostrings)

# Startup ends here

## _ Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

# Author: Jacqueline May
# Updated: October 4th, 2023 by Jessica Castellanos Labarcena 
# Updated: October 3rd, 2024, by Karl Cottenie
# This is a tutorial for getting NCBI data into R using the package "rentrez".

# References and further reading:
# https://www.ncbi.nlm.nih.gov/books/NBK25501/
# https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
# https://docs.ropensci.org/rentrez/

# PART 1: LOAD PACKAGES ----

# Install the "rentrez package".
install.packages("rentrez")
library(rentrez)
install.packages("seqinr")
library(seqinr)

# PART 2: GETTING DATABASE INFO ----

# Use the entrez_dbs() function to retrieve the names of databases that are available through EUtils.
?entrez_dbs
entrez_dbs()

# The following functions will provide you with more information about each of these databases.

# Returns summary information about the database, including name of the database, description, last update etc.
?entrez_db_summary
entrez_db_summary(db = "snp")

# Returns a list of search fields that can be used for a particular database.
# These are the fields that are available in the Advanced Search Builder on the NCBI browser.
?entrez_db_searchable
entrez_db_searchable(db = "genome")

# Returns a list of other databases that contain cross-referenced records for a particular database.
?entrez_db_links
taxonomy_links <- entrez_db_links(db = "taxonomy")
# Convert to dataframe format.
dfTaxLinks <- as.data.frame(taxonomy_links)
View(dfTaxLinks)

rm(taxonomy_links, dfTaxLinks)

# Feel free to try some on your own! #

# PART 3: USING SEARCH TERMS ----

# Let's try searching a database using the entrez_search() function.
?entrez_search

# "Search a given NCBI database with a particular query". So, let's search PubMed for records that match the key term "Cotesia"!
# db = name of the database to search for, term = the search term.
search_result <- entrez_search(db = "pubmed", term = "Cotesia")

# Check the class.
class(search_result)

# It is an "esearch" object and a list. Let's look at what it contains.
search_result

# We can see that "All Fields" were searched in the pubmed database for the term "Cotesia". There are more 'hits' for this search than IDs in this object because of the argument "retmax". If you look at the documentation for entrez_search(), you will see that retmax refers to the maximum number of hits returned by the search. Its default value is 20, but of course we can change that.

# Let's look at the components of search_result.
search_result$ids
# What kind of class do you think it will be?
class(search_result$ids)
# How many IDs do you think there are?
length(search_result$ids)

search_result$count ## But we had 619 hits! Let's increase the argument retmax to 619

# Assign the max number of hits to a variable called "hits".
maxHits <- search_result$count
# Set our retmax argument to maxHits.
search_result_2 <- entrez_search(db = "pubmed", term = "Cotesia", retmax = maxHits)

length(search_result_2$ids) ## There are now 619 IDs.

# What about using more specific search terms? We can search specific fields and use the boolean operators AND, OR, and NOT to expand and/or restrict our searches.

# Let's search the Nucleotide database ("nuccore") for some Salmonidae (salmon, trout, etc.) sequences!

# What fields are available to search for the Nucleotide Database?
entrez_db_searchable("nuccore")

# ORGN = Scientific and common names of organism, and all higher levels of taxonomy.
salmon_search <- entrez_search(db = "nuccore", term = "Salmonidae[ORGN]", retmax = 100)

# Let's take a look at salmon_search.
salmon_search$ids  ## Unique IDs returned by the search.
salmon_search$count  ## Total number of hits for the search.
salmon_search$retmax  ## Maximum number of hits returned by the search.

# Now, let's search for sequences within two specific Salmonidae genera (Oncorhynchus and Salmo).
genera_search <- entrez_search(db = "nuccore", term = "Oncorhynchus[ORGN] OR Salmo[ORGN]", retmax = 100)
genera_search

# Notice we used the "OR" operator in our search and not the "AND" operator. What happens if we use "AND"?
genera_search_AND <- entrez_search(db = "nuccore", term = "Oncorhynchus[ORGN] AND Salmo[ORGN]", retmax = 100)

genera_search_AND ## 0 hits!! Be careful with your operators.
 
# What if we want sequences from a certain date ([PDAT] field)? 
genera_search_test <- entrez_search(db = "nuccore", term = "Oncorhynchus[ORGN] OR Salmo[ORGN] AND 2016:2018[PDAT]", retmax = 100)
genera_search_test

# NOTE: You could also copy and paste your search from the NCBI browser into the argument for "term". This could be found in your History in the Advanced Search builder or under Search Details on the right side of the screen after you've run a search on NCBI.

# Remove objects we don't need anymore.
rm(maxHits, search_result, search_result_2, salmon_search, genera_search, genera_search_AND, genera_search_test)

# Let's search for data from a certain gene, cytochrome B (another mitochondrial gene).
# Note to be careful here when specifying the name of the Gene (e.g. Cyt B and CytB give you a different # of hits!).
cytb_search <- entrez_search(db = "nuccore", term = "Salmonidae[ORGN] AND CytB[Gene]", retmax = 100)
cytb_search

cyt_b_search <- entrez_search(db = "nuccore", term = "Salmonidae[ORGN] AND Cyt B[Gene]", retmax = 100)
cyt_b_search ## Comes down to how the researchers who entered the information spelled the gene name.

# When using entrez, it is important that you be as specific as possible. Because we said we want all gene sequences that contain CytB, we also end up getting whole mitochondrial genomes!
# IMPORTANT!!!: To make sure we only get CytB sequences, we could restrict our sequence length (the length you choose will depend on the gene you are searching):
cytb_search1 <- entrez_search(db = "nuccore", term = "(Salmonidae[ORGN] AND CytB[Gene] AND 600:1000[SLEN])", retmax = 100)
cytb_search1

# See the difference in the number of hits?
cytb_search
cytb_search1

# NOTE! GenBank is particular about the annotation of some genes, such as those for non-coding rRNA. If you are searching for rRNA like 18S,  you can use the PROP field:
EighteenS_search <- entrez_search(db = "nuccore", term = "Salmonidae[ORGN] AND 18S rRNA AND biomol_rRNA[PROP]", retmax = 100)
EighteenS_search

rm(cytb_search, cyt_b_search, EighteenS_search)

# PART 4: FETCHING RECORDS ----

# OK! But how do we get the data about the nucleotide sequences (and not just their IDs)? Well, there are two ways. entrez_fetch() returns full records in different formats, and entrez_summary() returns metadata about the sequence records.

# Using the entrez_summary function:
?entrez_summary
cytb_summ <- entrez_summary(db = "nuccore", id = cytb_search1$ids)  ## Notice that we are accessing the ids from the cytb_search2 object!
cytb_summ

class(cytb_summ)

# Take a closer look at some of the metadata of this particular ID.
cytb_summ$`2786555283`
class(cytb_summ$`2786555283`) ## Another list! 

# Diving deeper into the metadata..
# Genome information.
cytb_summ$`2786555283`$genome  ## Notice the double dollar signs!
class(cytb_summ$`2786555283`$genome)

# Organism information.
cytb_summ$`2786555283`$organism
class(cytb_summ$`2786555283`$organism)

cytb_summ |> 
  map_chr(~ .x$organism) |> 
  #unique()
  table()
  

# Remember lists are flexible and can contain objects of different types.

# Helper function called extract_from_esummary takes specific elements from each of your list elements.
?extract_from_esummary
extract_from_esummary(cytb_summ, "organism")

# Now, let's use entrez_fetch() to obtain sequence data in FASTA format. As useful as the summary records are, sometimes they just don’t have the information that you need. If you want a complete representation of a record you can use entrez_fetch, using the argument rettype to specify the format you’d like the record in.

?entrez_fetch 

# Let's specify that we want the return type ("rettype") in FASTA format.
cytb_fetch <- entrez_fetch(db = "nuccore", id = cytb_search1$ids, rettype = "fasta")

# What class is it?
class(cytb_fetch)
# It's an extremely long character vector.
head(cytb_fetch)

# Write to file in your current working directory, so you can keep a copy of the data (and double check it).
# Notice the new line (\n) deliminator.
write(cytb_fetch, "cytb_fetch.fasta", sep = "\n")  ## Check this file in a text editor! Is it what you expect? Are there REALLY long sequences? If so, you might have accidentally downloaded whole genomes. Make sure all of your sequences are of the expected length.

# Read it back in as DNA StringSet using the readDNAStringSet() function.
stringSet <- readDNAStringSet("cytb_fetch.fasta")
class(stringSet)


head(names(stringSet))

# Let's put these sequences into a dataframe for fun!
dfCytB <- data.frame(CytB_Title = names(stringSet), CytB_Sequence = paste(stringSet))
View(dfCytB)

# When looking at dfCytB, those names don't look very nice. We could clean them up using the handy word() function from the "stringr" package. You could also use regular expressions here if you want (or any other function that will accomplish the same task).

# Make a new column, so you can keep the original sequence title.
dfCytB$Species_Name <- word(dfCytB$CytB_Title, 2L, 3L)
# Rearrange the columns.
dfCytB <- dfCytB[, c("CytB_Title", "Species_Name", "CytB_Sequence")]
# Let's look at the dataframe.
View(dfCytB)

# Now, let's do the same thing for the nuclear protein coding gene RAG1.
rag1_search <- entrez_search(db = "nuccore", term = "(Salmonidae[ORGN] AND RAG1[Gene] AND 1000:3000[SLEN])", retmax = 100)
# Fetch the sequences in FASTA format.
rag1_fetch <- entrez_fetch(db = "nuccore", id = rag1_search$ids, rettype = "fasta")
# Write to file.
write(rag1_fetch, "rag1_fetch.fasta", sep = "\n") 
# Read back in as DNAStringSet.
stringSet <- readDNAStringSet("rag1_fetch.fasta")
# Convert to dataframe format.
dfRAG1 <- data.frame(RAG1_Title = names(stringSet), RAG1_Sequence = paste(stringSet))
# Clean up species names.
dfRAG1$Species_Name <- word(dfRAG1$RAG1_Title, 2L, 3L)
# Rearrange the columns.
dfRAG1 <- dfRAG1[, c("RAG1_Title", "Species_Name", "RAG1_Sequence")]
# Let's look at the dataframe.
View(dfRAG1)

# Note: in some cases your sequence titles may be a lot messier, and you would have to use regex to extract the species names.

# For now, let's randomly sample one sequence per species (using tidyverse)!!!
length(unique(dfCytB$Species_Name)) 

dfCytB_Subset <- dfCytB %>% 
  group_by(Species_Name) %>%  ## Group by Name.
  sample_n(1)  ## Handy function in dplyr for randomly selecting rows from tbl object.

# There should be 4 rows (one per species) in this new subset for Cytb.
# We can confirm this:
all.equal(length(unique(dfCytB$Species_Name)), nrow(dfCytB_Subset))

# Let's do the same for the RAG1 gene. There should be 32 rows. 
length(unique(dfRAG1$Species_Name))

dfRAG1_Subset <- dfRAG1 %>% 
  group_by(Species_Name) %>% 
  sample_n(1)

# Confirming it worked:
all.equal(length(unique(dfRAG1$Species_Name)), nrow(dfRAG1_Subset))

# We can merge two dataframes together using the merge() function. 
?base::merge

# We are matching rows in dfCytB that have the same names in the dfRAG1 subset (basically finding the overlapping records). BUT let's specify all = TRUE so we still keep those species that don't have sequences available for both genes.
dfAllSeqs <- merge(dfCytB_Subset, dfRAG1_Subset, by = "Species_Name", all = T)  ## Can just specify T instead of writing out TRUE.
View(dfAllSeqs)

# alternative way, using tidyverse verbs
dfCytB_Subset |> 
  full_join(dfRAG1_Subset)

# If we specified all = F (the default), we would only be keeping those rows with names in BOTH dataframes (i.e. complete cases).
dfOverlap <- merge(dfCytB_Subset, dfRAG1_Subset, by = "Species_Name", all = F)
View(dfOverlap)

# tidyverse way
dfCytB_Subset |> 
  inner_join(dfRAG1_Subset)

# PART 5: WEB HISTORY OBJECTS ----

# NOTE: If you want to extract large amounts of sequences from NCBI (e.g. all of the hits from your search), you can use a web_history object by setting the "use_history" argument to TRUE. When you are dealing with very large queries it can be time consuming to pass long vectors of unique IDs to and from the NCBI. To avoid this problem, the NCBI provides a feature called “web history” which allows users to store IDs on the NCBI servers then refer to them in future calls.For example:
squamata_search <- entrez_search(db = "nuccore", term = "Squamata[ORGN] AND RAG1[Gene] AND 1000:3000[SLEN]", use_history = T)

squamata_search
squamata_search$web_history ## This stores the ID on the NCBI server.

# If you are interested, I have provided some functions I wrote for extracting FASTA files using web_history objects on CourseLink. If you want to load them into your session you would put the R file into your current working directory and type: 
source("Entrez_Functions.R")

# Example if you want to try it yourself later (commented out for now because it may take a while):

# Fetching the sequences in FASTA format and writing to file.
#FetchFastaFiles(searchTerm = "Squamata[ORGN] AND RAG1[Gene] AND 1000:3000[SLEN]", seqsPerFile = 100, fastaFileName = "Squamata_RAG1")

# Merging the sequences into a dataframe.
# <- MergeFastaFiles(filePattern = "Squamata_RAG1*")


# PART 6: CHALLENGE! ----

# 1. Decide as a group a taxon you would like to search the "nuccore" database for. Each group member can then search for sequence records for a different gene for that taxon. Share with your group information about your search terms and the number of hits you received for that gene.

# 2. Try to get the data for your gene into dataframe format like we did in part 4. The dataframe should contain columns for the original title of the sequences (e.g. CytB title), the sequence itself (e.g. CytB sequence), and species name.

# 3. As a group, decide on a feature of your sequence data that you are interested in. For example, the sequence lengths or the "AT" proportions of your sequences. Create a new column in your dataframe that holds this information. Then, compare with your group the average value of that feature for your gene.


