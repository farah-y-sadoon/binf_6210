###### **********************
## 
## Lecture 12 - Machine Learning and Random Forest
## In class activity question 12
##
##***************************

##_Loading Libraries ----
library(tidyverse)
library(randomForest)
library(Biostrings)

#_Question 12 - Build own classifier to distinguish genera ----
###_ Load data and investigate ----
df_wasps<- read_tsv("../data/wasps.tsv")
names(df_wasps)
summary(df_wasps)
table(df_wasps$markercode)

df_wasps_cleaned <- df_wasps %>%
  filter(!is.na(nucleotides)) %>% 
  filter(!is.na(genus_name))

unique(df_wasps_cleaned$genus_name)
sum(is.na(df_wasps_cleaned$nucleotides))

###_ Create a new data frame ----
df_sequences <- df_wasps_cleaned %>%
  mutate(nucleotides2 = str_remove(nucleotides, "^[-N]+")) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  filter(str_count(nucleotides2, "N") <= (0.05 * str_count(nucleotides)))

df_compare_seqs <- cbind(df_sequences$nucleotides, df_sequences$nucleotides2)

q1 <- quantile(nchar(df_sequences$nucleotides2[df_sequences$genus_name == "Cotesia"]), probs = 0.25, na.rm = TRUE)
q1
q3 <- quantile(nchar(df_sequences$nucleotides2[df_sequences$genus_name == "Cotesia"]), probs = 0.75, na.rm = TRUE)
q3

table(df_sequences$genus_name)

df_sequences <- df_sequences %>%
  filter(((str_count(nucleotides2) >= q1 & str_count(nucleotides2) <= q3 & genus_name == "Cotesia") | genus_name == "Campoletis"))
dim(df_sequences)
unique(df_sequences$genus_name)
table(df_sequences$genus_name)
sum(is.na(df_sequences$nucleotides2))
summary(str_count(df_sequences$nucleotides2[df_sequences$genus_name == "Cotesia"]))
summary(str_count(df_sequences$nucleotides2[df_sequences$genus_name == "Campoletis"]))


###_ Create model ----
## Converting tibble to data frame for biostrings to use
df_sequences <- as.data.frame(df_sequences)
df_sequences$nucleotides2 <- DNAStringSet(df_sequences$nucleotides2)
class(df_sequences$nucleotides2)

df_sequences <- cbind(df_sequences, as.data.frame(letterFrequency(df_sequences$nucleotides2, letters = c("A", "C","G", "T"))))

## Calculating proportions of each nucleotide
df_sequences$Aprop <- (df_sequences$A) / (df_sequences$A + df_sequences$T + df_sequences$C + df_sequences$G)

df_sequences$Tprop <- (df_sequences$T) / (df_sequences$A + df_sequences$T + df_sequences$C + df_sequences$G)

df_sequences$Gprop <- (df_sequences$G) / (df_sequences$A + df_sequences$T + df_sequences$C + df_sequences$G)

df_sequences$nucleotides2 <- as.character(df_sequences$nucleotides2)

table(df_sequences$genus_name)

## Create validation data set
set.seed(217)

df_validation <- df_sequences %>%
  group_by(genus_name) %>%
  sample_frac(0.25)

table(df_validation$genus_name)

## Create training data set
set.seed(13)
df_training <- df_sequences %>%
  filter(!processid %in% df_validation$processid) %>%
  group_by(genus_name) %>%
  sample_n(1024)

table(df_training$genus_name)

names(df_training)

## Create final model
genus_classifier <- randomForest::randomForest(x = df_training[, 86:88], y = base::as.factor(df_training$genus_name), ntree = 50, importance = TRUE)

genus_classifier

genus_classifier$importance

## Validate
predict_validation <- predict(genus_classifier, df_validation[, c(70, 86:88)])
predict_validation

## Confusion matrix
table(observed = df_validation$genus_name, predicted = predict_validation)
