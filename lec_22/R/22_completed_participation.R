##***************************
##  COMPLETED Software Tools - Functions, debugging, and iteration
##
## Karl Cottenie
##
## 2025-11-20
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())


# Startup ends here

## _ FUNCTIONS ----

dfMice = read_tsv("../data/mice.tsv")
dfMice

dfWhales = read_tsv("../data/whales.tsv")
dfWhales

# The following function should calculate the number of species per bin_uri

fn_uniqueSpecies <- function(df) {
  df.unique <- df %>%
    group_by(species_name) %>%
    summarize(num_species = length(unique(species_name))) %>%
    arrange(desc(num_species))
  return(df.unique)
}

dfMice %>% fn_uniqueSpecies()
dfWhales %>% fn_uniqueSpecies()

##### _ FIX THIS FUNCTION ----

fn_uniqueSpecies <- function(df) {
  df.unique <- df %>%
    group_by(bin_uri) %>%
    summarize(num_species = length(unique(species_name))) %>%
    arrange(desc(num_species))
  return(df.unique)
}

##### _ SMART CHECK ----

# create a tibble with values I can sanity check manually
check_df <- tibble(
  bin_uri = c("bin1", "bin2", "bin3", "bin1", "bin1"),
  species_name = c("a", "b", "c", "b", "c")
)

fn_uniqueSpecies(check_df) # apply function to new, smaller data frame

## _ FOR LOOP ----

# First create a list with your data sets
ls_mammals = list(Mice = dfMice, Whales = dfWhales)
ls_mammals
ls_mammals$Mice %>% class()
ls_mammals[1] %>% class()
ls_mammals[[1]] %>% class()

for( i in 1:2) {
  fn_uniqueSpecies(ls_mammals[i])
}

##### _ FIX THIS CODE ----
for( i in 1:2) {
  fn_uniqueSpecies(ls_mammals[[i]]) %>% 
  print()
}


