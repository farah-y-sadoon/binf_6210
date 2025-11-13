##***************************
##  Software Tools - Functions
##
## Karl Cottenie
##
## 2025-10-27
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
#install.packages("usethis")

##_ Fork and Clone ------

# this will create a new project in RStudio in a new window

usethis::create_from_github(
  "https://github.com/karl-cottenie/clustering_2025.git",
  # replace the above github url with the green Code > Local > Clone using the web URL
  # you need to update this for each repository you want to contribute to
  destdir = "./", 
  # or whatever directory location you want
  # this destdir can be the same across different repositories if you want all your github projects to be in the same location
  # each repository will be in a different subdirectory in this location
  fork = TRUE
)