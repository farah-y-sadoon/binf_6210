##***************************
##  QUIZ #5 FUNCTIONS AND ITERATIONS
##
## Karl Cottenie
##
##
## 2025-11-20
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

## -------------------------------------------------------------------------------------------------------- ##

##_ Q1. READ PanTHERIA.tsv INTO SCRIPT ----
df_traits <- read_tsv("../data/PanTHERIA.tsv")

# Adjusting names in the data frame
names_original <- as.vector(names(df_traits))
names_edited <- str_replace(names_original, "^MSW05_|^[0-9]+-[0-9]+_", "")
names(df_traits) <- names_edited

##_ Q2. FUNCTION THAT CALCULATES WHETHER BERGMANN'S RULE FOR A GIVEN CUT-OFF VALUE OF NON-TROPICAL ----
cutoff_lat = seq(10, 50, 10)

fn_bergmann <- function(df, cutoff_lat) { 
  # filter data frame
  df_filtered <- df %>%
    filter(!is.na(GR_MidRangeLat_dd)) %>% 
    filter(!is.na(AdultBodyMass_g)) %>%
    filter(GR_MidRangeLat_dd >= abs(cutoff_lat)) %>%
    group_by(Genus) %>% 
    filter(n() >= 10) %>%
    ungroup()
  
  # create a list of data frames for each genus
  ls_traits_berg_genus <- df_filtered %>% 
    nest_by(Genus) %>% 
    pull(data, name = Genus)
  
  # run a linear model for body mass as a function of mid-range latitude for each genus
  ls_models_berg <- map(ls_traits_berg_genus, 
                        function(x) lm(log10(x$AdultBodyMass_g) ~ abs(x$GR_MidRangeLat_dd)))
  
  # extract the coefficients from each linear model - we are interested in the slopes for each data frame, aka each genus
  df_models = ls_models_berg %>% 
    map(coef) %>% 
    map(as_tibble_row) %>%
    list_rbind(names_to = "Genus")
  
  # conduct t-test on the slopes for each genus and extract p-values - this is the return value for the function
  p_value <- t.test(df_models[, 3], mu = 0)$p.value
  
  return(p_value)
}

##_ Q3. ITERATE OVER CUTOFF VALUES ----
# iterate over the cutoff values list to extract p-values from each mid-range latitude cutoff for each genus
p_values <- map_dbl(cutoff_lat, \ (x) fn_bergmann(df_traits, x))

# save results and cutoff values in a data frame to plot later
df_bergmann_results <- tibble(cutoff_lat = cutoff_lat, p_value = p_values)
df_bergmann_results

##_ Q4. CREATE A FIGURE THAT ILLUSTRATES WHETHER THE SIGNIFICANCE OF BERGMANN'S RULE DEPENDS ON THE DEFINITION OF NON-TROPICAL ----

ggplot(df_bergmann_results, aes(x = cutoff_lat, y = p_value)) + 
  geom_line() + 
  geom_point(size = 3) + 
  geom_text(aes(label = round(p_value, 3), nudge_x = 1.5)) + 
  labs(title = "Impact of Cutoff Values on p-Values for Bergmann's Rule", 
       x = "Cutoff Latitude Value",
       y = "p-Value of Slope from t-tests"
  )
             