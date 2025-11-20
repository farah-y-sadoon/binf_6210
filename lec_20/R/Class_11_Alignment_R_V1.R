######SOFTWARE TOOLS CLASS 11 - DNA SEQUENCE ALIGNMENT V2----
###### **********************
## SOFTWARE TOOLS CLASS 11 - DNA SEQUENCE ALIGNMENT V1
##
## Karl Cottenie
##
## 2024-10-06
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

#last updated October 06, 2024

#Script adapted from Sally Adamowicz. Resources consulted include the documentation for the R package muscle as well as the original paper (Edgar 2004. Nucleic Acids Research) and online documentation for the MUSCLE algorithm.

#Overview: This script provides an introduction to DNA sequence alignment using tools available through the R package muscle. As you go, fill in the group activity sheet. Afterwards, on your own, try out these functions on data of interest to you for your Assignment #2. Also, on your own and for your assignment, try out the alignment functions available through the DECIPHER package. V2 contains example answers for you to consider (but note that there is typically more than one way to do things in R).

####PART 1 - LOAD REQUIRED PACKAGES----

#Packages from CRAN

#install.packages("tidyverse")
library(tidyverse)

#install.packages("stringi")
library(stringi)

#install.packages("ape")
library(ape)

install.packages("RSQLite") #a required package for DECIPHER
library(RSQLite)

#Tip: If you have trouble with package installation for CRAN packages, note that occasionally there can be a problem with a particular CRAN mirror site. I have found the following function helpful to choose a different CRAN mirror. I would suggest to try a different mirror within North America.

#chooseCRANmirror()

#Packages from Bioconductor.

#First, install and load the Bioconductor installation system, if needed.
#install.packages("BiocManager")
#library(BiocManager)

#Then, install any needed packages.
#BiocManager::install(c("Biostrings", "muscle", "msa", "DECIPHER"))

#Load libraries
library(Biostrings)
library(muscle)
library(DECIPHER)
#library(msa) - not used here, you may wish to use in the future

#Tip: be sure to READ error messages. Some packages have dependencies. You may need to install additional required packages. I suggest to "update all" if you have some out-of-date packages. And, remember to load libraries with every R session.

####PART 2 - IMPORT AND CHECK DATASET----

#Download from CourseLink and import the file.
#The below should work if you have set your working directory to the location where the below file is stored on your computer. We are using a familiar data set (from class 9) so that we can move on promptly.
setwd("./")
getwd()
dfBOLD <- read_tsv("../data/wasps.tsv")

#Check the data set size, data parsing and data summary, and variable names.
dim(dfBOLD)
summary(dfBOLD)
names(dfBOLD)

####PART 3 - FILTERING, DATA FORMATTING, AND PRELIMINARY DATA EXPLORATION----

#You may wish to select key columns to make your data frame smaller and easier to look at (e.g. the unique identifiers such as processid, taxonomy, markercode column, and nucleotides.)

#Tip: It is always important to keep a unique identifier. In this case, we will keep processid, which is unique for each specimen.

#Here, we are selecting the specified columns that we want, to be included in the new dataframe, using the select() function from the dplyr package. Note that this is an easy way to be able to use the column names.

#Coding challenge #1- Adapt the below line of code so that you also keep the country and lat and lon GPS coordinates.
dfBOLD1 <- dfBOLD %>%
  select(processid, bin_uri, genus_name, species_name, markercode, nucleotides, country, lat, lon)

#Have a look at object in viewer.

#Challenge #2: Filter data set to include only records having a nucleotide sequence present for the 28S marker. 
df28S <- dfBOLD1 %>%
  filter(markercode == "28S") %>% 
  filter(!is.na(nucleotides))

#check properties
dim(df28S)
names(df28S)
unique(df28S$markercode)
sum(is.na(df28S$nucleotides))

#Here, we are removing original data frame to clean up our workspace.
rm(dfBOLD)

#Build a summary of sequence lengths for data exploration.
summary(nchar(df28S$nucleotides))

#Challenge #3: make a histogram to display the distribution of sequence lengths.
df28S <- df28S %>% 
  mutate(sequence_length = nchar(nucleotides))

plot(hist(df28S$sequence_length))

#Tidyverse notes: For the summary, we could instead use the function str_length() from the tidyverse package stringr. For example:
summary(str_length(df28S$nucleotides))

#This will also count the number of characters and then we are passing that info to summary(). The empty pattern argument for str_count() will result in counting of characters.
summary(str_count(df28S$nucleotides))

#By contrast, here we are setting the pattern argument to something specific. This will count the number of "A"s, as an example. I would encourage you to read the document on stringr saved under the "R Cheatsheets" folder in CourseLink. This is a nice set of functions for string handling.
mean(str_count(string = df28S$nucleotides, pattern = "A"))

#Same thing but shorter (we can omit the names of arguments when using them in order; omitting argument names requires care and is most suitable when there are only one or a few arguments for a given function.)
mean(str_count(df28S$nucleotides, "A"))

#Tip: You may notice that what we just did above is very similar to the letterFrequency() function from Biostrings that we used the other day, but the letterFrequency() function will act upon a DNAStringSet (i.e. a Bioconductor style data object containing character data, following the S4 object-oriented programming system). So, as usual, there is more than one way to get to the same answer in R! I would suggest to use Bioconductor packages for large data sets, due to the more efficient way of representing large biological data.

#Conclusion about the sequence variability: So, there is some variability in sequence length in this dataset, which we may wish to consider during downstream analysis. However, we are NOT seeing VERY extreme values (e.g. a whole genome sequence or very long or very short sequences accidentally present in this dataset). If we mix sequences of very different length, then this would interfere with multiple sequence alignment. In cases where we genuinely have a short query sequence and long reference sequence (e.g. mapping short sequence reads to a reference genome), then we would need to choose an alignment method suitable for that (a LOCAL alignment option). Here, we will perform a multiple sequence alignment of a single gene, so we will choose a suitable alignment method for that purpose.

#Change the format of the nucleotide data to be suitable for analysis by functions in the muscle package, which is a Bioconductor package.
class(df28S)
df28S <- as.data.frame(df28S)
df28S$nucleotides2 <- DNAStringSet(df28S$nucleotides)
class(df28S)
class(df28S$nucleotides2)

#Setting processid as the sequence names, so that we are sure to carry forward a unique identifier for the nucleotides through the analysis steps.
names(df28S$nucleotides2) <- df28S$processid
names(df28S$nucleotides2)

#Viewing on screen
df28S$nucleotides

#We can use a great function from DECIPHER package to get a nicer view of the sequences in a browser.
BrowseSeqs(df28S$nucleotides2)

####PART 4 - MULTIPLE SEQUENCE ALIGNMENT----

#Using the muscle algorithm from the muscle package, let's perform a multiple sequence alignment.

#We will try different analysis settings to see what happens when we dramatically change the gap penalties. Note that the online manual for the muscle software has more detail than the muscle R package documentation. See: http://www.drive5.com/muscle/muscle_userguide3.8.html

#Here is an example. Here, we wanted a fast preliminary alignment to check if this runs, and so I am setting maxiters (maximum iterations) to 2. See Edgar paper on the algorithm. This will truncate the number of alignment refinement iterations. So, this is suitable for a fast, preliminary alignment, but I would remove that constraint for a final analysis. After alignment, we are converting to a DNAStringSet and assigning this object to df28S.alignment.
df28S.alignment <- DNAStringSet(muscle::muscle(df28S$nucleotides2, maxiters = 2), use.names = TRUE)

#Viewing subset on screen
df28S.alignment

#Can also view and browse our alignment in html format using the function BrowseSeqs() from the DECIPHER package. You may also consider the function msaPrettyPrint() from the msa package, but you need to install a LaTeX system and also ensure you have "safe" file and path names (i.e. without special characters). So, that function can be trickier to get going. By contrast, BrowseSeqs() is very user-friendly.
BrowseSeqs(df28S.alignment)

#Note that there is only partial help regarding the muscle package and function directly available through R. The following doesn't bring up all the options.
?muscle
#However, the help does link to the MUSCLE website, which explains the options in more detail.
#MUSCLE website: http://www.drive5.com/muscle/muscle_userguide3.8.html

#For example, for the maxiters option: "You can control the number of iterations that MUSCLE does by specifying the -maxiters option. If you specify 1, 2 or 3, then this is exactly the number of iterations that will be performed. If the value is greater than 3, then muscle will continue up to the maximum you specify or until convergence is reached, which ever happens sooner. The default is 16. If you have a large number of sequences, refinement may be rather slow." From: http://www.drive5.com/muscle/muscle_userguide3.8.html

#removing preliminary alignment
rm(df28S.alignment)

#Try different analysis settings....

#What happens when we dramatically change the gap penalties?

#For now, we won't alter the defaults, allowing the alignment to improve over multiple iterations if needed. We are requesting a log file so that we can see the defaults used in this R implementation of the muscle algorithm. We can open the log file, which was saved in the working directory. The default gap opening penalty is -400 and the default for maxiters is 8.
df28S.alignment1 <- DNAStringSet(muscle::muscle(df28S$nucleotides2, log = "log.tx", verbose = T), use.names = T)

#checking class of our resulting alignment so we know what we are working with. Our alignment is also a DNAStringSet.
class(df28S.alignment1)

#Viewing
BrowseSeqs(df28S.alignment1)

#let's have a look at the first element.
nchar(df28S.alignment1[1])

#How long is our DNAStringSet?
length(df28S.alignment1)

#Here, we are assessing the length of the first element (i.e. sequence length). We can use indexing with double square brackets here, similar to the case of lists and data framees, to dive into element #1 to get at the data.
length(df28S.alignment1[[1]])

#Coding challenge #4: Write a line of code to determine how many gaps there are in sequence 1 in the alignment.
####your code here.

countPattern("-", df28S.alignment1[[1]])

#Now, let's try making a different alignment, with a very different gap opening penalty. This is penalizing the opening of gaps more heavily for the calculation of the alignment score. i.e. Setting the gap opening penalty to be more negative means that it costs more to open a gap; i.e., more mismatches will be tolerated before a gap is opened.
df28S.alignment2 <- DNAStringSet(muscle::muscle(df28S$nucleotides2, gapopen = -10000), use.names = T)

#Viewing alignment
BrowseSeqs(df28S.alignment2)

#Now, checking length of the alignment (i.e. number of nucleotide positions in length), using sequence 1.
length(df28S.alignment2[[1]])

#So, we can see that the total alignment lengths from the same starting sequences can be different depending upon the gap penalty. We could also count the number of gaps or look at the distribution of gap counts.

#What if we want to look at the average number of gaps across all sequences in the alignment? Fortunately, our alignment (DNAStringSet object) is of a class that can be treated like a list. Therefore, we can use lapply() - hooray! Here, we want lapply() because we want to apply a function over the elements of a list. This is short way to code a loop-like structure in R. The entire apply family of functions is very helpful.

#Here, we are applying the function str_count(), with the argument set to the gap character (i.e. the string pattern for which we want to count occurrences), across all elements of our DNAStringSet object (our alignment).
lapply(as.character(df28S.alignment1), str_count, "-")

#lapply applies a function over the elements of a vector or list. Note that the arguments for the function we want to apply (in this case str_count) need to be placed after a comma, within the arguments to be passed to lapply. So, the syntax here is a little different than what we have seen so far. Check out the documentation for the base R function.
?lapply

#What if we want to find the mean gap count?

#If we want to find the mean of all of the gap counts across list elements, we could do the following. The below is first using lapply() to apply the function str_count() to all of the elements in our alignment. We are counting the number of gaps for each sequence in the alignment. Then, we can use unlist() to change the data type (from a list to a vector, a numeric vector in this case). This will enable us to pass the gap counts to the function mean() to calculate an overall mean. Note that the below doesn't save anything. Below, we are outputting the result to the screen. We could readily add a name and an assignment operator at the beginning if we wanted to save the results for further analysis. There are multiple ways to do this! You may notice I gravitate towards base R, but you don't have to.
mean(unlist(lapply(as.character(df28S.alignment1), str_count, "-")))

#Could do summary
summary(unlist(lapply(as.character(df28S.alignment1), str_count, "-")))

#The following is a way to do the same thing in tidyverse style, making use of piping.
as.character(df28S.alignment1) %>%
  lapply(str_count, "-") %>%
  unlist %>%
  mean

#Repeating the calculation of the mean number of gaps, for alignment 2
mean(unlist(lapply(as.character(df28S.alignment2), str_count, "-")))

#Looking at a distribution of gap counts in the two alignments, using a simple histogram. (Note that you might wish to use ggplot instead for preparing nicer figures for assignments, reports, websites, and publications.)
hist(unlist(lapply(as.character(df28S.alignment1), str_count, "-")))
hist(unlist(lapply(as.character(df28S.alignment2), str_count, "-")))

#So, we now have an idea about how changing the gap penalty influenced some features of the resulting alignments. Setting a heavier gap penalty resulted in a shorter, less-gapped alignment. We would also expect increasing nucleotide dissimilarity at certain positions in the alignment. We could examine this additional expectation through looking at the distribution of genetic distances among sequences under different alignment scenarios (we will cover genetic distance calculation next time).

#But, how would we choose a suitable method and gap penalty? For conserved genes from closely related species, we would typically get a good alignment with standard approaches because there are relatively few indels. However, for more distant relatives, it would be advisable to use more sophisticated alignment methods. For example, for rRNA we can use alignments that consider secondary structure. For protein-coding genes, it is preferable to align the sequences in amino acid format and then untranslate back to nucleotides (or use a codon-informed alignment method). Try out those options using DECIPHER functions!

#Note: we can write our alignment to a file to view in other software if you wish (such as MEGA: https://www.megasoftware.net/). There are multiple different packages that have a function for writing fasta files, e.g.:
writeXStringSet(df28S.alignment1, file = "Cotesia_alignment1_28S.fasta")
