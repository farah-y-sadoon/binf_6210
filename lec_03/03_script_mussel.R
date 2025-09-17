##***************************
## mussel 0 problem
##
## Karl Cottenie
##
## 2024-02-08
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())
library(vegan)
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)

# Startup ends here

## _ Comment codes ------ %>% 
# Coding explanations (#, often after the code, but not exclusively)
# Code organization (## XXXXX -----)
# Justification for a section of code ## XXX
# Dead end analyses because it did not work, or not pursuing this line of inquiry (but leave it in as a trace of it, to potentially solve this issue, or avoid making the same mistake in the future # (>_<) 
# Solutions/results/interpretations (#==> XXX)
# Reference to manuscript pieces, figures, results, tables, ... # (*_*)
# TODO items #TODO
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

# _ data input ------

df_mussel = read_csv("TRIAGE_2021_data_mussel crew data_for Karl.csv")
df_mussel

#or use the below line, without the #

load("mussel.RData")

# Cleaning up column names
df_mussel = df_mussel %>% as_tibble %>% janitor::clean_names()

# Data Description
df_mussel %>% skimr::skim()

# Data validity - missing values at the data set level
df_mussel %>% na.exclude() %>% skimr::skim()

df_mussel %>% filter(!complete.cases(.)) %>% View()

# Univariate data validity - 0s, outliers, invalid values
df_mussel %>% dataMaid::visualize()
dev.off()

df_mussel %>% dataMaid::check()

save(df_mussel, df_mussel_wo_comb, file = "mussel.RData")

# _ analysis ---------

df_mussel = df_mussel |> 
  select(u_imbecilis:f_mitcheli) |> 
  rowSums() |> 
  bind_cols(df_mussel) |> 
  rename(total = ...1) |> 
  mutate(river = str_sub(site_id, 1, 2), 
         .after = site_id)
df_mussel  

df_mussel |> 
  ggplot(aes(mesohabitat, total^0.25)) +
  geom_boxplot() +
  facet_grid(river ~ .)
#higher densities in pools compared to riffles?

df_mussel |> 
  filter(mesohabitat != "Combined") %>%
lm(total^0.25 ~ mesohabitat * river, data = .) |> 
  anova()
# interaction not significant, mesohabiat almost sign, river sign!

df_mussel |> 
  filter(mesohabitat != "Combined") |> 
  select(u_imbecilis:f_mitcheli) |> 
  vegdist(method = "eucl") |> 
  stepacross(path = "extended")
#==> doesn't work, too many differences

df_mussel_wo_comb = df_mussel |> 
  select(u_imbecilis:f_mitcheli) |> 
  specnumber() |> 
  bind_cols(df_mussel) |> 
  filter(mesohabitat != "Combined") |> 
  rename(richn = ...1)

df_mussel_wo_comb |> 
  ggplot(aes(mesohabitat, richn)) +
  geom_boxplot() +
  facet_grid(river ~ .)

# higher densities in pools compared to riffles?

df_mussel_wo_comb %>%
  lm(richn ~ mesohabitat * river, data = .) |> 
  anova()

# _ second approach: dropping samples w/o ind------------

# __ check sample size ------
df_mussel_compl = df_mussel_wo_comb |> 
  filter(total > 0)
  #==> dropped from 182 to 87, huge loss of data!

df_mussel_compl |> 
  select(u_imbecilis:f_mitcheli) |> 
  metaMDS() |> 
  autoplot()
  
# __ multivariate analysis -------
ls_mussel = rda(decostand(select(df_mussel_compl,
                     u_imbecilis:f_mitcheli),
              method = "hellinger") 
    ~ mesohabitat + Condition(river), data = df_mussel_compl)
ls_mussel
autoplot(ls_mussel)
anova(ls_mussel)

RsquareAdj(ls_mussel)
