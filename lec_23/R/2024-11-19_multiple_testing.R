##***************************
## Problem of multiple comparisons
##
## Karl Cottenie
##
## 2024-11-19
##
##***************************

## _ Packages used -------
library(stats)
library(tidyverse)
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

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

## _ Starting point from Law et al. -----

# 3 samples from the basal and 3 samples from the MP group
# approximately 16000 potential genes
# average log-CPM = 5, spread -5 to 15
# => model log-CPM as a normal distribution with mean 5 and sd = 3

## _ One t-test, without any difference in gene expression ----

set.seed(5)

vc_basal = rnorm(3, mean = 5, sd = 3)
vc_mp = rnorm(3, mean = 5, sd = 3)

# mean is the same => no difference in gene expression 
vc_basal |> sort() 
vc_mp |> sort()

t.test(vc_basal, vc_mp, alternative = "two.sided")
# ===> p-value = 0.6
# no significant difference

## _ Repeat this analysis 16000 times

set.seed(5)
t.test(rnorm(3, mean = 5, sd = 3),
       rnorm(3, mean = 5, sd = 3), 
       alternative = "two.sided")$p.value
#this is the p-value of the t-test
#with the set.seed(5) it is always 0.5883

# now repeat the above analysis 16000 times

vc_pvalues = map_dbl(1:16000, # repeat my function 16000 times
                     function(x){
                       # I want different samples for each replication 
                       # => no set.seed!
                       t.test(rnorm(3, mean = 5, sd = 3),
                              rnorm(3, mean = 5, sd = 3), 
                              alternative = "two.sided")$p.value
                     })

vc_pvalues |> 
  boxplot() 
abline(h = 0.05)
# ===> there are p-values below the 0.05 cut-off line
# remember, in line 63-64, the mean is exactly the same!!

## _ How many are below 0.05?----

sum(vc_pvalues < 0.05)
# 519 significant p-values while samples are not significantly different!!!!!
# slightly different from the theoretically predicted 800 significant values

