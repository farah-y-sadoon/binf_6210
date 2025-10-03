#_Setting up Libraries --------
library(tidyverse)
library(tidyr)
library(dplyr)
library(janitor)
library(stringr)
library(sf)
library(ggplot2)

#_Import Data --------

# Use BOLD API to extract data
  #df_anura <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Anura&format=tsv")
  #write_tsv(df_anura, "../data/df_anura.tsv")
  
  getwd()

#_Inspect Data --------
  df_anura <- read_tsv("../data/df_anura.tsv")
  
  # Define objects and ensure they are as expected
  class(df_anura)
  dim(df_anura)
  head(df_anura)
  tail(df_anura)
  summary(df_anura)
  names(df_anura)
  
  # Create a subset of the data for columns of interest
  df_anura2 <- df_anura %>% 
    select(processid, bin_uri, order_name, family_name, genus_name, species_name, lat, lon, country, province_state)
  
  # read_tsv("../data/df_anura2.tsv")
  dim(df_anura2)
  names(df_anura2)
  #write_tsv(df_anura2, "../data/df_anura2.tsv")
  
  # Display the number of records that exists for qualitative columns
  View(df_anura2 %>% 
         count(bin_uri, sort = TRUE))
  View(df_anura2 %>% 
         count(order_name, sort = TRUE)) # make sure only the order of interest is reported here, otherwise remove in cleaning steps
  View(df_anura2 %>% 
         count(family_name, sort = TRUE))
  View(df_anura2 %>% 
         count(genus_name, sort = TRUE))
  View(df_anura2 %>% 
         count(species_name, sort = TRUE)) # interesting that many records do not report species names
  View(df_anura2 %>% 
         count(country, sort = TRUE))

  # Analyse numeric columns of interest by checking distribution and ranges
  hist(df_anura2$lat)
  hist(df_anura2$lon)
  summary(df_anura2$lat, na.rm = TRUE) 
  summary(df_anura2$lon, na.rm = TRUE) # interesting that longitude range seems to cover the entire globe ()

  # Look at the number of countries that have the most bin data
  df_anura2 %>%
    count(country, sort = TRUE)
  
  # Make sure columns are in the datatypes we need them
  str(df_anura2)
  
#_Cleaning Data --------
  
  # Filter fields create a new column for region
   df_anura_cleaned <- df_anura2 %>%
    filter(!is.na(lat)) %>%
    filter(!is.na(lon)) %>%
    filter(!is.na(bin_uri)) %>%
    filter(lat > -66.50 & lat < 66.50) %>%
    mutate(climate_region = ifelse(lat < 23.51 & lat > -23.51, "tropical", "temperate"))
  
  # Clean names to make sure everything is in snake_case and filter duplicates
  clean_names(df_anura_cleaned)
  
  dim(df_anura_cleaned)
  df_anura_cleaned <- df_anura_cleaned %>%
    distinct() # interesting - before this there were approx. 3000 duplicate records
  dim(df_anura_cleaned)
  
  # Check that the new column categorizes region as required
  df_temp_check_region <- df_anura_cleaned %>%
    group_by(climate_region) %>%             # Replace 'region' with your column name
    summarise(min_lat = min(lat), max_lat = max(lat), med_lat = median(lat), n_records = n()) %>%
    arrange(min_lat)
    
    df_temp_check_region
    
  # Remove temporary data frame
  rm(df_temp_check_region)

  # Look at new cleaned data frame
  dim(df_anura_cleaned)
  
  # Display the number of records that exists for qualitative columns
  View(df_anura_cleaned %>% 
         count(bin_uri, sort = TRUE))
  df_anura_cleaned %>% 
         count(order_name, sort = TRUE) # make sure only the order of interest is reported here, otherwise remove in cleaning steps
  View(df_anura_cleaned %>% 
         count(family_name, sort = TRUE))
  View(df_anura_cleaned %>% 
         count(genus_name, sort = TRUE))
  View(df_anura_cleaned %>% 
         count(species_name, sort = TRUE)) # interesting that many records do not report species names
  View(df_anura_cleaned %>% 
         count(country, sort = TRUE))
  
  # Analyse numeric columns of interest by checking distribution and ranges
  hist(df_anura_cleaned$lat)
  hist(df_anura_cleaned$lon)
  summary(df_anura_cleaned$lat)
  summary(df_anura2$lon) # longitude range still seems to cover the entire globe
  
  # Look at plotted coordinates to visualize the spread of data
  ggplot(df_anura_cleaned, aes(x = lon, y = lat, colour = climate_region)) +
    geom_point(alpha = 0.1) +
    coord_fixed()

#_Data Analysis --------
##_Question 1: How does biodiversity of anurans differ between temperate and tropical regions?
  