###### **********************
## CLASS 8 SOFTWARE TOOLS - INTRODUCTION TO BIOSTRINGS AND K-MERS
##
## Karl Cottenie
##
## 2024-09-29
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

#Script adapted from Sally Adamowicz

#Resources consulted include package vignettes and function documentation

####1 - OVERVIEW----

#In this tutorial, we will commence working with DNA sequences in R using packages available through Bioconductor, starting with the Biostrings package. Bioconductor is a curated repository with a focus upon bioinformatics tools.

#Bioconductor is separate from CRAN (Comprehensive R Archive Network), so we need to specify the source for obtaining packages through Bioconductor. I would recommend in general to use mainly CRAN and Bioconductor tools when we are beginning in R, as submission of R packages to these repositories requires software tests and documentation. Once you are more comfortable, you can feel free to explore tools under development that have been posted to GitHub (but you would have to be more careful about quality and ensuring the tool does what you need).

#This tutorial uses a DNA sequence dataset from a group of parasitoid wasps (genus Cotesia). We will first put our data exploration and filtering skills to use to create two data subsets, each containing sequences from a single marker gene (COI and 28S).

####2 - INSTALLING AND LOADING NEEDED PACKAGES----

#Packages available through CRAN:

#tidyverse suite of packages used for data manipulation in this script
#install.packages("tidyverse")
library(tidyverse)

#Packages available through Bioconductor. This is the code for installation for R version 3.5 or greater.
#Instructions from: https://bioconductor.org/install/
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install(c("Biostrings", "DECIPHER", "muscle"))

#For today, we are only using Biostrings
library(Biostrings)

#We can explore the functions available through Biostrings by accessing the vignettes. The point for now is so that you know the vignettes exist. Later, you can uncomment and run the below lines to check out the vignettes on your own.

#Seeing what vignettes are available.
#vignette()

#Checking what vignettes are available specifically for the package Biostrings
#vignette(package = "Biostrings")

#viewing the quick overview vignette
#vignette("BiostringsQuickOverview", package = "Biostrings")

#Viewing the Pairwise Alignments vignette
#vignette("PairwiseAlignments", package = "Biostrings")

#Viewing Multiple Alignment vignette
#vignette("MultipleAlignments", package = "Biostrings")

####3 - LOADING DATA FILE AND CHECKING DATA PARSING----

#Cotesia is an interesting genus of parasitoid wasps. This taxonomic group was the subject of a paper by University of Guelph researcher Alex Smith and colleagues in 2008 (published in PNAS), which is posted to CourseLink. We will use this dataset for exploring some DNA sequence properties. I would encourage you to complete your own explorations of sequences on other data sets of interest to you afterwards.

#This file was downloaded from BOLD on October 2, 2021. You would uncomment this line if you would like to run this or a similar API call yourself. Note if that if you run this yourself, the specific answers will vary. So, I suggest to use the same file as posted on CourseLink for now.
#Cotesia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cotesia&format=tsv")

#I used the following line to write the file to hard disk. You can download the file from CourseLink, under Class 8.
#write_tsv(Cotesia, "Cotesia_BOLD_data.tsv")

#Now, we can use the same commands, starting from here to read in the data. Again, note I am using a generic data frame name, to make this script easier to adapt for diverse taxonomic groups.
dfBOLD <- read_tsv("../data/Cotesia_BOLD_data.tsv")

#You should now have a data frame in your global environment containing 5162 observations of 80 variables.

#Summary of the data and structure of the data in Cotesia. It is always wise to check data parsing.
summary(dfBOLD)

#column names
names(dfBOLD)

#This is a familiar format by now! So, let's move on ...

####4 - REVIEWING PIPING, EXPLORING SAMPLE SIZE, FILTERING THE DATA----

#Checking class
class(dfBOLD)

#Next, we are seeing what markers are available in this dataset. Our rows in BOLD data are individual specimens sequenced for one marker. Often, a specimen only has one row, but a specimen may have more than one row if multiple genes have been sequenced from the same specimen. The below line will tell us what unique markers are in the dataset. Note that the "markercode" variable indicates which marker is housed in that row. (It is perhaps a bit confusing that there is a different variable, called "marker_codes", which contains data about all of the primers used for that specimen. So, ensure you use "markercode" variable here.)
unique(dfBOLD$markercode)

#Next, let's obtain a count of the sequences by markercode.
dfBOLD.by.marker <- dfBOLD %>%
  filter(!is.na(nucleotides)) %>%
  count(markercode) %>%
  arrange(desc(n)) %>%
  print()

#We can see that most of the data are COI-5P (the 5' end of the cytochrome c oxidase subunit I mitochondrial gene), which is the standard marker gene for DNA barcoding of animals. The next most, at a distant second, is 28s (the nuclear large ribosomal RNA subunit gene). Remember that BOLD focuses on standardized markers, so this is expected that complementary markers have lower representation.

#Creating COI dataset. Filtering to retain only those rows with the COI-5P marker and that have nucleotide data present. We're also trimming off trailing gaps (-) and Ns at the beginnings and ends of the sequences. We are also removing all gaps throughout because we don't want aligned sequences; here, we want k-mer profiles based upon biological nucleotide sequences (not gap characters that are introduced during alignment). For now, we will keep internal Ns, as those represent a nucleotide that is likely present but of unknown identity. Later, we will see that that choice has consequences, so you will need to revisit that issue of how to treat the Ns. We are retaining sequences 400 - 700 bp in length to constrain the analysis and avoid very long or very short sequences, which could have different properties compared to the standard DNA barcode region of COI. You will try the remaining of the script first while, first, commenting out line 97 and then, second, uncommenting and running the script again. Ensure that you carefully read and thoroughly understand each line of code.
dfCOI <- dfBOLD %>%
  filter(markercode == "COI-5P") %>%
  filter(!is.na(nucleotides)) %>%
  mutate(nucleotides2 = nucleotides) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "^[-N]+")) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  #filter(str_count(nucleotides2, "N") <= (0.02 * str_count(nucleotides2))) %>%
  filter(str_count(nucleotides2) <= 700) %>%
  filter(str_count(nucleotides2) >= 400)

#checking summary statistics for COI sequence lengths to ensure the above worked, this time using the base R function nchar().
summary(nchar(dfCOI$nucleotides2))

#showing what which.max() does. Gives index position for the maximum value.
which.max(nchar(dfCOI$nucleotides2))

#have a look at sequence with maximum length
dfCOI$nucleotides2[which.max(nchar(dfCOI$nucleotides2))]

#Checking that we have only COI sequences remaining in our new COI dataset and no records missing a sequence.
unique(dfCOI$markercode)
sum(is.na(dfCOI$nucleotides2))

#looking at distribution of sequence lengths for 28S, using our original data frame dfBOLD
summary(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "28S"]))
hist(nchar(dfBOLD$nucleotides[dfBOLD$markercode == "28S"]))

#Similarly, creating a dataframe containing 28S sequences, with constrained sequence length.
df28S <- dfBOLD %>%
  filter(markercode == "28S") %>%
  filter(!is.na(nucleotides)) %>%
  mutate(nucleotides2 = str_remove(nucleotides, "^[-N]+")) %>%
  mutate(nucleotides2 = str_remove(nucleotides2, "[-N]+$")) %>%
  mutate(nucleotides2 = str_remove_all(nucleotides2, "-+")) %>%
  filter(str_count(nucleotides2) <= 700) %>%
  filter(str_count(nucleotides2) >= 400)

#checking df28S
summary(nchar(df28S$nucleotides2))
unique(df28S$markercode)
sum(is.na(df28S$nucleotides2))

#Now we can remove unneeded objects.
rm(dfBOLD.by.marker)

####5 - CALCULATING DNA SEQUENCE FEATURES----

#Checking the data type of the nucleotides, character.
class(dfCOI$nucleotides2)

#However, character is not the correct data type for usage with Biostrings functions. We can see that when viewing the help for functions we want to try out from the Biostrings vignettes.

#Below, we are converting the format of the nucleotide sequences to a DNAStringSet. A DNAStringSet contains character (string) data but belongs to a stricter object class. If you wish to read further about this, Biostrings uses the "S4 object system" rather than the "S3" system. You can read more about that topic in various posts about Bioconductor as well as in "Advanced R" by Hadley Wickham, which I'd recommend to read after "R for Data Science": https://adv-r.hadley.nz/

#First, we convert to regular data frame. (We could have instead gone with base R style programming above and not used tibbles at all throughout this script. However, those functions are very handy for data filtering! Also, conversion from a tibble to a regular data frame can be easily completed as follows.)
dfCOI <- as.data.frame(dfCOI)
class(dfCOI)
df28S <- as.data.frame(df28S)
class(df28S)

#Changing the class of our trimmed nucleotide data to a DNAStringSet.
dfCOI$nucleotides2 <- DNAStringSet(dfCOI$nucleotides2)

#Checking the class of the new nucleotides column
class(dfCOI$nucleotides2)

#Also converting for the 28S marker dataset
df28S$nucleotides2 <- DNAStringSet(df28S$nucleotides2)

#checking that the class converted.
class(df28S$nucleotides2)

#Nucleotide composition is a fundamental feature of DNA sequences.

#For accessing the help in Biostrings, we leave off the normal () we have usually used so far after the function name. Here, we want to access the documentation from the function letterFrequency.
?letterFrequency

#This below line will calculate the nucleotide count of the specified letter, in this case A, in each sequence. We are assigning the results to a data frame so that we can easily look at the results.
dfCOI.AFreq <- as.data.frame(letterFrequency(dfCOI$nucleotides2, letters = "A"))

#In your environment window, click on the name of the new object, and have a look.

#Also, for comparison, let's try with the argument as.prob = TRUE. That would be important if we are comparing the composition of sequences of different lengths. See the difference? This will take into account sequence length and report decimal proportion of As out of the total sequence length, not the absolute count of As.
dfCOI.AFreq <- as.data.frame(letterFrequency(x = dfCOI$nucleotides2, letters = "A", as.prob = TRUE))

#We can also ask for the frequency for more than one letter, in this example A or T frequency (see documentation).
dfCOI.ATFreq <- as.data.frame(letterFrequency(dfCOI$nucleotides2, letters = "AT", as.prob = TRUE))

#We could instead do A and T separately if we want. Compare dfCOI.ATFreq and dfCOI.ATFreq2 using the viewer.
dfCOI.ATFreq2 <- as.data.frame(letterFrequency(dfCOI$nucleotides2, letters = c("A", "T"), as.prob = TRUE))
rm(dfCOI.ATFreq2)

#By contrast, the function alphabetFrequency calculates the frequency of all of the nucleotides separately. Have a look in the viewer. (Note there are various letter codes to denote uncertain nucleotides. See: https://www.bioinformatics.org/sms/iupac.html Most of these are unused in the present data set.)
dfCOI.NucFreq <- as.data.frame(alphabetFrequency(x = dfCOI$nucleotides2, as.prob = TRUE))

#Let's look at the summary and histogram of AT frequencies. (Note that we need quotation marks here because the column header has a special character.)
summary(dfCOI.ATFreq$"A|T")
hist(dfCOI.ATFreq$"A|T")

#Let's do the same for 28S.
df28S.ATFreq <- as.data.frame(letterFrequency(df28S$nucleotides2, letters = "AT", as.prob = TRUE))

#viewing data summary and histogram.
summary(df28S.ATFreq$"A|T")
hist(df28S.ATFreq$"A|T")

#Is there a relationship between sequence length and AT frequency? (Answer: not really here. Yes, there was before I added the filter to remove very short or very long sequence lengths, as different genetic regions can have different properties... on your own later, try playing with the filtering code above to see what happens!)
plot(nchar(dfCOI$nucleotides2), dfCOI.ATFreq$"A|T")
plot(nchar(df28S$nucleotides2), df28S.ATFreq$"A|T")

#Given this dataset, are the mean AT frequencies significantly different between COI and 28S in this data set? Answer: Yes.
t.test(df28S.ATFreq$"A|T", dfCOI.ATFreq$"A|T")

#But, let's check for outliers. Remember to use "Zoom" in RStudio to see plot. Looks like a more severe outlier in COI.
boxplot(dfCOI.ATFreq$"A|T")
boxplot(df28S.ATFreq$"A|T")

#let's look at the sequence of COI with the lowest AT frequency. Note that for the DNAStringSet, colour is used.
dfCOI$nucleotides2[which.min(dfCOI.ATFreq$"A|T")]

#same sequence using the regular character nucleotides column.
dfCOI$nucleotides[which.min(dfCOI.ATFreq$"A|T")]

#We see a lot of Ns in that sequence! That would likely be driving the AT frequency anomaly here. We should consider adding an additional filter above, to remove sequences with a lot of Ns. This is another example where we have to think carefully about how to treat missing data. Would we want to delete this entire sequence? Trim off the beginning part with a lot of Ns? Try to build a filter that removes sequences with, for example, 2% Ns (or another threshold). As species typically differ by more than 2% sequence divergence at COI, then we may not want sequences with more than 2% Ns, as that might confound identifications. Have a look again after adding the filter on line 97.

#Calculating k-mers. Below, we are calculating the counts of all of the short oligonucleotide sequences of length 4. Again, we could choose whether we want absolute counts or whether we want a probability (i.e. a decimal proportion considering total sequence length). Have a look at oligo4 in the viewer.
oligo4 <- oligonucleotideFrequency(x = dfCOI$nucleotides2, width = 4, as.prob = TRUE)

#Let's try a different k-mer length and then have a look at the object. Notice the relative size of those two matrices! The size of the data grows rapidly with increasing k!
oligo3 <- oligonucleotideFrequency(x = dfCOI$nucleotides2, width = 3, as.prob = TRUE)
