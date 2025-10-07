#_Setting up Libraries --------
library(tidyverse)
library(tidyr)
library(dplyr)
library(janitor)
library(stringr)
library(sf)
library(ggplot2)
library(vegan)
library(iNEXT)
library(countrycode)

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
    select(processid, bin_uri, order_name, lat, lon, country)
  
  # read_tsv("../data/df_anura2.tsv")
  dim(df_anura2)
  names(df_anura2)
  write_tsv(df_anura2, "../data/df_anura2.tsv")
  
  # Display the number of records that exists for qualitative columns
  View(df_anura2 %>% 
         count(bin_uri, sort = TRUE))
  View(df_anura2 %>% 
         count(order_name, sort = TRUE)) # make sure only the order of interest is reported here, otherwise remove in cleaning steps
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
  
  # Filter fields create a new column for latitude bands
  # Check if there are any records at the equator exactly
  df_check_equator <- df_anura2 %>% 
    filter(lat == 0) # no observations returned, no specimen were recorded on the equator exactly
  rm(df_check_equator)

  df_anura_cleaned <- df_anura2 %>%
    filter(!is.na(lat)) %>%
    filter(!is.na(lon)) %>%
    filter(!is.na(bin_uri)) %>%
    filter(!is.na(country) & country != "Unrecoverable") %>%
    mutate(hemisphere = ifelse(lat >= 0, "Northern", "Southern")) %>% # add hemisphere for later analysis of biodiversity differences
    mutate(lat_band = case_when(
      lat > 0 & lat < 10   ~ "0-10°N",
      lat >= 10 & lat < 20  ~ "10-20°N",
      lat >= 20 & lat < 30  ~ "20-30°N",
      lat >= 30 & lat < 40  ~ "30-40°N",
      lat >= 40 & lat < 50  ~ "40-50°N",
      lat >= 50 & lat < 60  ~ "50-60°N",
      lat >= 60 & lat < 70  ~ "60-70°N",
      lat >= 70 & lat < 80  ~ "70-80°N",
      lat >= 80 & lat <= 90 ~ "80-90°N",
      lat < 0 & lat > -10   ~ "0-10°S",
      lat <= -10 & lat > -20 ~ "10-20°S",
      lat <= -20 & lat > -30 ~ "20-30°S",
      lat <= -30 & lat > -40 ~ "30-40°S",
      lat <= -40 & lat > -50 ~ "40-50°S",
      lat <= -50 & lat > -60 ~ "50-60°S",
      lat <= -60 & lat > -70 ~ "60-70°S",
      lat <= -70 & lat > -80 ~ "70-80°S",
      lat <= -80 & lat >= -90 ~ "80-90°S",
      TRUE ~ "INVALID"))
  
  # Add continent for later analysis of biodiversity - NA and unrecoverable country codes are filtered out earlier for countrycodes()
  df_anura_cleaned <- df_anura_cleaned %>% 
    mutate(continent = countrycode(country, origin = "country.name", destination = "continent"))
  
  # Check to see if any latitude bands were not captured by the case_when clauses
  df_check_lat_band <- df_anura_cleaned %>% 
    filter(lat_band == "INVALID") # no observations returned, no specimen in the cleaned data frame have invalid latitude values
  rm(df_check_lat_band)
  
  # Check that the new column categorizes lat_band as required, and that each band has enough records for further analysis
  df_check_lat_band_range <- df_anura_cleaned %>%
    group_by(lat_band) %>%             # Replace 'region' with your column name
    summarise(min_lat = min(lat), max_lat = max(lat), med_lat = median(lat), n_records = n()) %>%
    arrange(min_lat)
  rm(df_check_lat_band_range)
  
  # Clean names to make sure everything is in snake_case and filter duplicates
  clean_names(df_anura_cleaned)
  
  # Remove duplicate values, if any exist
  dim(df_anura_cleaned)
  
  df_anura_cleaned <- df_anura_cleaned %>%
    distinct() # interesting - before this there were approx. 3000 duplicate records
 
  dim(df_anura_cleaned)
  
  # Display the number of records that exists for qualitative columns
  View(df_anura_cleaned %>% 
         count(bin_uri, sort = TRUE))
  df_anura_cleaned %>% 
         count(order_name, sort = TRUE) # make sure only the order of interest is reported here, otherwise remove in cleaning steps
  View(df_anura_cleaned %>% 
         count(country, sort = TRUE))
  View(df_anura_cleaned %>% 
         count(continent, sort = TRUE))
  
  # Analyse numeric columns of interest by checking distribution and ranges
  hist(df_anura_cleaned$lat)
  hist(df_anura_cleaned$lon)
  summary(df_anura_cleaned$lat)
  summary(df_anura_cleaned$lon) # longitude range still seems to cover the entire globe
  
  # Look at plotted coordinates to visualize the spread of data
  ggplot(df_anura_cleaned, aes(x = lon, y = lat, colour = lat_band)) +
    geom_point(alpha = 0.1) +
    coord_fixed()
  
  write_tsv(df_anura_cleaned, "../data/df_anura_cleaned")
  
#_Data Exploration & Analysis --------
##_How does biodiversity between tropical and temperate regions differ? --------
  ####### PREVIOUS WORK
  # Spread the data
  df_bins_by_region <- df_anura_cleaned %>% 
    group_by(climate_region, bin_uri) %>% 
    count()
  
  df_bins_by_region_spread <- pivot_wider(data = df_bins_by_region, names_from = bin_uri, values_from = n, values_fill = 0)
  
  # Calculate the common sample size
  df_count_climate_region <- df_anura_cleaned %>% 
    count(climate_region)
  
  common_sample_size <- min(df_count_climate_region$n)
  
  # Create community matrix
  df_bins_community_matrix <- df_bins_by_region_spread[,-1] #removes region column
  
  # Calculate the number of observed species at the common sample size
  climate_regions <- df_bins_by_region_spread$climate_region # store this for combining with rarified richness data
  rarefied_richness <- rarefy(df_bins_community_matrix, common_sample_size)
  
  df_rarefied_by_climate_region <- data.frame(region = climate_regions,
                                              rarefied_richness = rarefied_richness)
  
  count_projects <- df_anura_cleaned %>% 
    mutate(project_id = str_extract(processid, "^[A-Za-z]+")) %>% 
    count(project_id)
  
  
  count_projects %>% 
    count(project_id) %>% 
    summary(count(project_id))
  hist(count_projects$n, breaks = 50, main = "Samples per project", xlab = "Number of samples")
  
##_What is the relationship between anuran biodiversity and latitude? --------
  
  # Test for normality on independent variable (latitude)
  qqnorm(df_anura_cleaned$lat)
  qqline(df_anura_cleaned$lat, col = "red")
  
  set.seed(123)
  sample_lat <- sample(df_anura_cleaned$lat, 5000) # Shapiro test doesn't allow for a sample larger than 5000
  shapiro.test(sample_lat) # Data is skewed, so cannot use mean as the value representing each latitude band
  rm(sample_lat)

  # Impute latitude band values with median instead of categories or mean for regression testing since latitude is skewed
  # Extract median latitude from all of the latitude points
  df_anura_lat_median <- df_anura_cleaned %>%
    group_by(lat_band) %>%
    summarize(lat_median = median(lat)) %>%
    mutate(lat_median = lat_median)
  
  df_anura_analysis <- df_anura_cleaned %>% 
    left_join(df_anura_lat_median, by = "lat_band")
  rm(df_anura_lat_median)
  
  # Transform data from cleaned data frame into a matrix for further analysis
  df_anura_analysis_counts <- df_anura_analysis %>% 
    group_by(lat_median, bin_uri) %>%
    count()
  
  # Pivot wider to make BINs column names
  df_anura_analysis_counts_spread <- pivot_wider(data = df_anura_analysis_counts, names_from = bin_uri, values_from = n, values_fill = 0)
  str(df_anura_analysis_counts_spread)
  
  # Label rows with median latitude because matrix has to ignore the latitude values
  rownames(df_anura_analysis_counts_spread) <- df_anura_analysis_counts_spread$lat_median
  
  # Create a matrix that removes the "lat_median" label
  mat_anura_abundance <- as.matrix(df_anura_analysis_counts_spread)
  str(mat_anura_abundance)
  mat_anura_abundance <- mat_anura_abundance[, -1]
  str(mat_anura_abundance)
  
  rm(df_anura_analysis_counts_spread, df_anura_analysis_counts)
  
  # Plot rarefaction curves for each median lat point to determine sampling completeness
  rarecurve(mat_anura_abundance,
            main = "Rarefaction Curves",
            label = FALSE)
  
  # Format labels to avoid overlap
  sample_sums <- rowSums(mat_anura_abundance)
  specimen_sums <- specnumber(mat_anura_abundance)

  offsets <- seq(-4, 4, length.out = length(sample_sums))
  
  text(x = sample_sums + 50,
       y = specimen_sums + offsets,
       labels = rownames(mat_anura_abundance),
       cex = 0.7,
       col = "black",
       xpd = TRUE)
  
  # Sampling not complete for some latitudes, therefore drop lowest site values and compute Hill numbers for each site 
  # Drop bands that do not fall within the 25th percentile
  site_totals <- rowSums(mat_anura_abundance)
  common_size <- quantile(site_totals, 0.25)
  
  # Combine with site names (latitude medians)
  df_site_totals <- data.frame(
    lat_median = as.numeric(rownames(mat_anura_abundance)),
    total_specimens = site_totals)
  
  
  
  
  
  
  common_size <- quantile(site_totals, 0.25)
  
  
  
  
  
  

