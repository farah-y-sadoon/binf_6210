_Setting up Libraries --------
library(tidyverse)
library(tidyr)
library(dplyr)
library(janitor)
library(stringr)

#_Import Data --------

# Use BOLD API to extract data
  df_anura <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Anura&format=tsv")

#_Inspect Data --------
  
  # Define objects and ensure they are as expected
  class(df_anura)
  head(df_anura)
  tail(df_anura)
  summary(df_anura)
  names(df_anura)
  
  # Create a subset of the data for columns of interest
  df_anura2 <- df_anura[, c(1, 8, 14, 16, 20, 22, 47, 48, 55, 56, 57, 72)]
  
  # Display the number of records that exists for qualitative columnS
  View(df_anura2 %>% count(bin_uri, sort = TRUE))
  View(df_anura2 %>% count(order_name, sort = TRUE)) # make sure only the order of interest is reported here, otherwise remove in cleaning steps
  View(df_anura2 %>% count(family_name, sort = TRUE))
  View(df_anura2 %>% count(genus_name, sort = TRUE))
  View(df_anura2 %>% count(species_name, sort = TRUE)) # interesting that many records do not report species names
  View(df_anura2 %>% count(country, sort = TRUE))
  View(df_anura2 %>% count(nucleotides, sort = TRUE))

  # Visualize distribution for numeric columns relevant to analysis
  hist(df_anura2$lat)
  hist(df_anura2$lon)
  
  