##***************************
## BINF 6410 Software Tools - Assignment 4
##
## Farah Sadoon 1302190
##
## 2025-12-05
##
##***************************

##_Load Libraries ----
library(tidyverse)
library(janitor)
library(picante)
library(ape)
library(purrr)
library(httr)

##_Loading, Exploring and Cleaning Data ----

# Great lakes fish species data frame
# Downloaded from https://www.glerl.noaa.gov/data/waterlife/glwlexplorer.html 
df_fish <- read_tsv("../data/fish_great_lakes_glerl_noaa.tsv")

# Explore great lakes fish species data frame
class(df_fish)
dim(df_fish)
head(df_fish)
tail(df_fish)
str(df_fish)
names(df_fish)

# Check if any species names are missing
table(df_fish$Species)

# Create a subset of the fish data set that includes only the columns of interest
df_fish2 <- df_fish %>%
  dplyr::select("Species", "Superior", "Michigan", "Huron", "Erie", "Ontario") %>%
  mutate(across(c(Superior, Michigan, Huron, Erie, Ontario),
                ~ replace_na(.x, "0"))) %>% # convert NA values to 0 for all columns listed
  mutate(across(c(Superior, Michigan, Huron, Erie, Ontario),
                ~ ifelse(.x == "X", 1, 0))) %>%  # convert X to numeric 1 to indicate presence, and any other value, i.e., "Extirpated" to 0 for absence
  filter(!if_all(c(Superior, Michigan, Huron, Erie, Ontario),
         ~ .x == 0)) %>% # removes any rows where the species is not present in any of the lakes
  clean_names() # converts all column names to snake_case

# Inspect cleaned data frame
head(df_fish2)
print(df_fish2$species) # noticed few species names were not correctly formatted

# Adjust species names in data frame manually
df_fish2$species <- df_fish2$species %>%
  recode(
    "Chrosomuseos" = "Chrosomus eos",
    "Chrosomuserythrogaster" = "Chrosomus erythrogaster",
    "Chrosomusneogaeus" = "Chrosomus neogaeus", 
    "Sandercanadensis" = "Sander canadensis", 
    "Sandervitreus" = "Sander vitreus", 
    "Lethenteronappendix" = "Lethenteron appendix"
  )

# Inspect the new data frame after names have been adjusted
head(df_fish2)
print(df_fish2$species)

# Use BOLD API to get nucleotide sequences for all fish - Phylum Chordata
# df_chordata <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Chordata&format=tsv")
# write_tsv(df_chordata, "../data/df_chordata.tsv")

df_chordata <- read_tsv("../data/df_chordata.tsv")

# Inspect chordate data frame
class(df_chordata)
dim(df_chordata)
head(df_chordata)
tail(df_chordata)
str(df_chordata)
names(df_chordata)

# Create a subset of the chordata data that includes the relevant information for this analysis
df_chordata2 <- df_chordata %>% 
  dplyr::select(species_name, markercode, nucleotides) %>% 
  filter(!is.na(species_name)) %>% 
  filter(markercode == "COI-5P") %>% 
  rename(species = "species_name") %>% # this needs to match the df_fish2 data frame for the left join
  distinct() %>% 
  clean_names()

# Inspect chordate data frame to see how it has changed after filtering
head(df_chordata2)
tail(df_chordata2)
names(df_chordata2)

table(df_chordata2$species_name)
table(df_chordata2$markercode) # all records should be COI-5P

# Join df_chordata2 with df_fish2
df_fish3 <- inner_join(df_fish2, df_chordata2, by = "species")

# Evaluate if  any rows/species were lost in the left join
length(unique(df_fish3$species)) == length(unique(df_fish2$species)) # this evaluates to false, so species were lost
length(unique(df_fish3$species))
length(unique(df_fish2$species)) # there are three more species here, which means that 3 species did not have sequence data

table(df_fish3$species) # since the minimum number of records for a unique species is 1, only keep one record per species to ensure no bias when creating phylogeny 

# Since there are multiple records for each species, create a data frame that only includes the records with the longest nucleotide sequence for each unique species
df_fish4 <- df_fish3 %>% 
  group_by(species) %>% 
  slice_max(order_by = nchar(nucleotides), with_ties = FALSE, n = 1) %>% # order each species by the max nucleotide sequence length and keep only the longest one, keep the first entry if there are ties
  ungroup()

# evaluate the number of entries for each species - should be one for each
table(df_fish4$species)