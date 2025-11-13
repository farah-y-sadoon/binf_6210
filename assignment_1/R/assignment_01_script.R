##***************************
## BINF-6210
##
## Assignment 1: Analysis of Anuran Biodiversity Across Latitude Bands
##
## Farah Sadoon
##
## 2025-10-10
##
## Main edits to code and suggestions for improvement by Eva Innocente on 12-11-25, with minor feedback from Fangyi Li and Lishita Rowjee
##
##***************************

rm(list = ls())

## _Setting up Libraries --------
library(tidyverse)
library(janitor)
library(sf)
library(vegan)
library(iNEXT)
library(countrycode)
library(maps)

## _Import Data --------

# Use BOLD API to extract data
#df_anura <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Anura&format=tsv")
# write_tsv(df_anura, "../data/df_anura.tsv")

getwd()

## _Inspect Data --------
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

## _Clean Data --------

# Filter fields create a new column for latitude bands
# Check if there are any records at the equator exactly
df_check_equator <- df_anura2 %>%
  filter(lat == 0) # no observations returned, no specimen were recorded on the equator exactly
rm(df_check_equator)

#df_anura_cleaned <- df_anura2 %>%
  #filter(!is.na(lat)) %>%
  #filter(!is.na(lon)) %>%
  #filter(!is.na(bin_uri)) %>%
  #filter(!is.na(country) & country != "Unrecoverable") %>%
  #mutate(hemisphere = ifelse(lat >= 0, "Northern", "Southern")) %>% # add hemisphere for later analysis of biodiversity differences
  #mutate(lat_band = case_when(
    #lat > 0 & lat < 10 ~ "0-10°N",
    #lat >= 10 & lat < 20 ~ "10-20°N",
    #lat >= 20 & lat < 30 ~ "20-30°N",
    #lat >= 30 & lat < 40 ~ "30-40°N",
    #lat >= 40 & lat < 50 ~ "40-50°N",
    #lat >= 50 & lat < 60 ~ "50-60°N",
    #lat >= 60 & lat < 70 ~ "60-70°N",
    #lat >= 70 & lat < 80 ~ "70-80°N",
    #lat >= 80 & lat <= 90 ~ "80-90°N",
    #lat < 0 & lat > -10 ~ "0-10°S",
    #lat <= -10 & lat > -20 ~ "10-20°S",
    #lat <= -20 & lat > -30 ~ "20-30°S",
    #lat <= -30 & lat > -40 ~ "30-40°S",
    #lat <= -40 & lat > -50 ~ "40-50°S",
    #lat <= -50 & lat > -60 ~ "50-60°S",
    #lat <= -60 & lat > -70 ~ "60-70°S",
    #lat <= -70 & lat > -80 ~ "70-80°S",
    #lat <= -80 & lat >= -90 ~ "80-90°S",
    #TRUE ~ "INVALID"
  #))


##### EDIT 1 by Eva:
# This is a simpler way to make bins for latitude, based on the highest and lowest latitude present in the dataset (this way you're not making bins with nothing in them). I used the palaeoverse package:

# install.packages("palaeoverse")
library(palaeoverse)

# checked what the min and max latitudes were

max(df_anura2$lat, na.rm=T)
min(df_anura2$lat, na.rm=T)

# filtered as before- I named this dataframe df_test

df_anura_cleaned <- df_anura2 %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(bin_uri)) %>%
  filter(!is.na(country) & country != "Unrecoverable") %>%
  mutate(hemisphere = ifelse(lat >= 0, "Northern", "Southern")) 

# First you create a dataframe of the bins

bins <- lat_bins_degrees(size = 10, min = -50, max = 70, fit = FALSE, plot = FALSE)

# Then run the bin_lat function to assign the bins you made and save it to a new object

df_anura_cleaned <- bin_lat(df_anura_cleaned, bins = bins, lat = "lat", boundary = FALSE)

# check how many specimens are in each latitude bin and making sure lat_bin is a factor, renaming to lat band to avoid confusion
df_anura_cleaned %>%
  group_by(lat_bin) %>%
  summarise(n())

df_anura_cleaned <- df_anura_cleaned %>%
  mutate(lat_band = as.factor(lat_bin))

#### end of edit

# Add continent for later analysis of biodiversity - NA and unrecoverable country codes are filtered out earlier for countrycodes()
df_anura_cleaned <- df_anura_cleaned %>%
  mutate(continent = countrycode(country, origin = "country.name", destination = "continent"))


# Check to see if any latitude bands were not captured by the case_when clauses
df_check_lat_band <- df_anura_cleaned %>%
  filter(lat_band == "INVALID") # no observations returned, no specimen in the cleaned data frame have invalid latitude values
rm(df_check_lat_band)

# Check that the new column categorizes lat_band as required, and that each band has enough records for further analysis
df_check_lat_band_range <- df_anura_cleaned %>%
  group_by(lat_band) %>% # Replace 'region' with your column name 
  summarise(min_lat = min(lat), max_lat = max(lat), med_lat = median(lat), n_records = n()) %>%
  arrange(min_lat)

#### EDIT 2 by Eva: 
# an effective way to show mumber of records per lat band could also be a bar chart of number of specimens in each lat band, or by grouping by the lat band variable and summarising

# making a quick rough bar chart to check numbers of specimens in each lat band

ggplot(data = df_anura_cleaned, aes(x = lat_band)) + geom_bar() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# summarising the data frame to check number of observations

df_anura_cleaned %>%
  group_by(lat_band) %>% # Replace 'region' with your column name 
  summarise(n())

#### end of edit

rm(df_check_lat_band_range)

#### EDIT 3 by Eva: 
# if you are just checking that something worked or looking at something in a dataframe (as above), it may be easier not to save it as an object so you don't have to keep removing objects, like I did in Edit 2. But this is up to personal preference! 

#### end of edit

# Clean names to make sure everything is in snake_case and filter duplicates
clean_names(df_anura_cleaned)

# Remove duplicate values, if any exist
dim(df_anura_cleaned)

df_anura_cleaned <- df_anura_cleaned %>%
  distinct() # approx. 3000 duplicate records, since there are multiple records for marker

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
#world <- map_data("world") # from maps() package

#(fig_distribution_specimen <- ggplot() +
  # Create world map outline
  #geom_polygon(
    #data = world, aes(x = long, y = lat, group = group),
    #fill = "gray100", color = "gray80"
  #) +
  # plot data points onto the map
  #geom_point(
    #data = df_anura_cleaned,
    #aes(x = lon, y = lat, colour = lat_band),
    #alpha = 0.6, size = 1
  #) +
  #coord_fixed(1.3) + # this ratio of 1.3 to 1 lat to lon to make the map look less stretched out
  #labs(
    #title = "Distribution of Anura Specimens Collected",
    #x = "Longitude", y = "Latitude", colour = "Latitude Band"
  #) +
  #theme_minimal())

#### EDIT 4 by Eva 

## The original figure with the datapoints coloured by lat band was slightly redundant as it showed a measure of latitude twice. I removed the legend specifying lat band and changed the y axis intervals to represent the lat bands instead. I changed point size and transparency to represent density of data. This way, we can see how many specimens were collected in a lat band by the transparency of points. If the colour is darker then it shows where a lot specimens were collected vs. fewer. Using the test dataframe I made above

### Changing figure:
world <- map_data("world")
(fig_distribution_specimen <- ggplot() +
     # Create world map outline
     geom_polygon(
       data = world, aes(x = long, y = lat, group = group),
       fill = "gray100", color = "gray80"
     ) +
     # plot data points onto the map
     geom_point(
       data = df_anura_cleaned,
       aes(x = lon, y = lat, colour = lat_band),
       alpha = 0.2, size = 0.8
     ) + scale_y_continuous(breaks = seq(-60, 80, by = 10)) + guides(colour = "none")
     + coord_fixed(1.3) # this ratio of 1.3 to 1 lat to lon to make the map look less stretched out
     + labs(
       title = "Distribution of Anura Specimens Collected",
       x = "Longitude", y = "Latitude", colour = "Latitude Band"
     ) +
     theme_minimal()) 

#### end of edit 

#ggsave("../figs/01_fig_distribution_of_anura_specimen.png", plot = fig_distribution_specimen, width = 12, height = 9, dpi = 800)

rm(world, fig_distribution_specimen)

write_tsv(df_anura_cleaned, "../data/df_anura_cleaned.tsv")

## _Data Exploration & Analysis --------
# What is the relationship between anuran biodiversity and latitude?
# Test for normality on independent variable (latitude)
qqnorm(df_anura_cleaned$lat)
qqline(df_anura_cleaned$lat, col = "red")


### EDIT 5 by Eva:
# Since the anura dataset has 13307 observations in it, and the Shapiro test only allows up to 5000, this test may not be the best for determining deviations from normality. Using a plot like the qqnorm() function is a better practice. The assumption of normality is important for larger datasets. Source: (Ghasemi and Zahediasl, 2012).

## commenting out these lines: 

#set.seed(123)
#sample_lat <- sample(df_anura_cleaned$lat, 5000) # Shapiro test doesn't allow for a sample larger than 5000
#shapiro.test(sample_lat) # Data is skewed, so cannot use mean as the value representing each latitude band
#rm(sample_lat)

## end of edit

# Impute latitude band values with median instead of categories or mean for regression testing since latitude is skewed
# Extract median latitude from all of the latitude points
df_anura_lat_median <- df_anura_cleaned %>%
  group_by(lat_band) %>%
  summarize(lat_median = median(lat)) %>%
  mutate(lat_median = lat_median)

df_anura_analysis <- df_anura_cleaned %>%
  left_join(df_anura_lat_median, by = "lat_band")
rm(df_anura_lat_median)

write_tsv(df_anura_analysis, "../data/df_anura_analysis.tsv")

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

# Write matrix to a table and specify parameters for delimiters and row/col names
write.table(mat_anura_abundance,
  file = "../data/mat_anura_abundance.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)


# Plot rarefaction curves for each median lat point to determine sampling completeness
png("../figs/02_fig_rarecurves.png", width = 3600, height = 2400, res = 300)
rarecurve(mat_anura_abundance,
  main = "Rarefaction Curves at Different Latitudes",
  label = FALSE
)

# Format labels to avoid overlap
sample_sums <- rowSums(mat_anura_abundance)
specimen_sums <- specnumber(mat_anura_abundance)
offsets <- seq(-4, 5, length.out = length(sample_sums))
text(
  x = sample_sums + 65,
  y = specimen_sums + offsets,
  labels = rownames(mat_anura_abundance),
  cex = 0.7,
  col = "black",
  xpd = TRUE
)

dev.off()

rm(sample_sums, specimen_sums, offsets)


# Sampling not complete for some latitudes, drop lowest site values and compute Hill numbers for each site
# Drop bands that do not fall within the 25th percentile
site_totals <- rowSums(mat_anura_abundance)
lat_cutoff <- quantile(site_totals, 0.25)
keep_sites <- site_totals >= lat_cutoff
mat_anura_filtered <- mat_anura_abundance[keep_sites, ]

rm(site_totals, lat_cutoff, keep_sites)

write.table(mat_anura_filtered,
  file = "../data/mat_anura_filtered.tsv",
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)

# Compare matrices before and after
dim(mat_anura_abundance)
dim(mat_anura_filtered)

# Split the matrix into lists where each list represents one latitude site and contains
mat_list <- split(mat_anura_filtered, rownames(mat_anura_filtered))
mat_list <- lapply(mat_list, as.numeric) # make each vector numeric

# Use iNEXT to calculate diversity indices
diversity_lat <- iNEXT(mat_list, q = c(0, 1, 2), datatype = "abundance")
class(diversity_lat$AsyEst)

rm(mat_list)

# Extract Asymptotic species richness from diversity_lat and create data frame for regression analysis
df_diversity_lat_analysis <- as_tibble(diversity_lat$AsyEst) %>%
  filter(Diversity == "Species richness") %>%
  dplyr::rename(
    latitude = "Assemblage", observed = "Observed", estimated_richness = "Estimator",
    standard_error = "s.e.", lcl = "LCL", ucl = "UCL"
  ) %>%
  select(-Diversity) %>%
  mutate(latitude = as.numeric(latitude)) %>%
  mutate(abs_lat = abs(latitude)) # Use absolute latitude for linear model

rm(diversity_lat)

write_tsv(df_diversity_lat_analysis, "../data/df_diversity_lat_analysis.tsv")

### _Statistical Testing ----

# Run linear regression to determine relationship between latitude and estimated asymptotic anura diversity
lm_diversity_lat <- lm(estimated_richness ~ abs_lat, data = df_diversity_lat_analysis)
summary(lm_diversity_lat)

# Assess residuals of the linear model
plot(lm_diversity_lat, which = 1) # Residuals of the model are not evenly distributed
shapiro.test(residuals(lm_diversity_lat)) # Residuals of the model are not normally distributed

# Log transform dependent variable to adjust for non-normal and skewed residuals
lm_log_diversity_lat <- lm(log(df_diversity_lat_analysis$estimated_richness) ~ df_diversity_lat_analysis$latitude, data = df_diversity_lat_analysis)
summary(lm_log_diversity_lat)

plot(lm_log_diversity_lat, which = 1) # Better homoscedasticity
shapiro.test(residuals(lm_log_diversity_lat)) # p > 0.05, this does not violate normal distribution

# Save output to a text file in output folder
capture.output(summary(lm_diversity_lat), file = "../output/01_output_lm_diversity_lat.txt")
capture.output(summary(lm_log_diversity_lat), file = "../output/02_output_lm_log_diversity_lat.txt")

rm(lm_diversity_lat, lm_log_diversity_lat)

# Plot log transformed linear regression
fig_log_linear_regression_diversity_lat <- ggplot(df_diversity_lat_analysis, aes(x = abs_lat, y = log(estimated_richness))) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(
    x = "Distance from the Equator (°)", y = "Log(Estimated Species Richness)",
    title = "Log Linear Model: Estimated Asymptotic Anuran Richness vs Absolute Latitude"
  ) +
  theme_minimal()
fig_log_linear_regression_diversity_lat

ggsave("../figs/03_fig_log_linear_regression_diversity_lat.png", plot = fig_log_linear_regression_diversity_lat, width = 16, height = 9, dpi = 800)