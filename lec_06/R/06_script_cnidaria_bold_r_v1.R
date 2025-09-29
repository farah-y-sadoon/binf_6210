##***************************
## Cnidaria - BOLD
##
## Karl Cottenie
##
## 2025-09-17
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
####1- INTRODUCTION----

#This script, as many of our scripts, has two main objectives: biological and programming skills

#BIOLOGICAL OBJECTIVE:

#We will use the example of the phylum Cnidaria to help us to learn how to analyze data from BOLD. for this group, here were more species than BINs (whereas the reverse is commonly expected for many invertebrate groups due to the prevalence of undescribed or cryptic species). Are there really more species of Cnidaria than BINs (genetic clusters defined based upon the COI data alone)? Or, could something else be causing the higher count of species? For example: perhaps Cnidaria have a lot of interim species names rather than formal Linnaean (scientific) species names. We will explore this topic in this example script.

#INFORMATICS OBJECTIVES: The above topic provides a foray into several important programming skills you will need for Assignment #1 and beyond:

#First, we will briefly review the indexing of matrices. We will also clarify that we can treat matrices as a vector, and show that we can pass a matrix with numerical data to core mathematical functions in R.

#introduction to loading contributed packages (i.e. additional R packages that did not come with base R)

#introduction to data frame

#indexing data frame by position

#indexing data frame by variable names

#The following topics will be continued in class 5:

#introduction to logical expression

#indexing data frame by condition

#intro to helpful basic functions that you can start using right away: e.g. table(), plot(), hist(), length(), unique(), sum(), is.na()

#brief preview of the idea of regular expressions as a very helpful tool for analyzing string/character data (NOTE: This is NOT on quiz 2; we will have a future lesson on regular expressions. This is a preview only.)

####2- MATRIX: INDEXING AND MATHEMATICAL FUNCTIONS----

#Create matrix of integers 1 to 100, with 10 by 10 dimensions
y <- matrix(1:100, 10, 10)
y

#check class
class(y)

#Let's see what happens when we pass one index position:
y[2]
y[20]
y[50]

#Check out length
length(y)

#matrix is a vector with a dimensions attribute
attributes(y)
dim(y)

#Consequence: we can pass a matrix containing numerical data to core mathematical functions that take a vector. This is very useful. For example:
mean(y)
max(y)
min(y)
median(y)
sd(y) #standard deviation

#Let's review two-dimensional indexing for selecting rows, columns, and elements
y[2, ] #row 2
y[, 2] #column 2
y[1:5, ] #rows 1 through 5, all columns
y[, 1:5] #all rows, columns 1 through 5
y[3, 5] #individual element (row 3, column 5)
y[c(1, 3, 5), c(2, 4, 6)] #can use c() to create vector of index positions for rows (here rows 1, 3, 5), columns (here, columns 2, 4, 6) or both

#check out class for the following:
class(y[2])
class(y[2, 2])
class(y[2, ])
class(y[, 2])
class(y[1:5, 1:5])

#remove unneeded objects at the end of each code section
rm(y)

####3- READING IN TSV FILE, INTRO TO DATA FRAME----

#If you have not done so already, uncomment (remove #) and run the following line.
#install.packages("tidyverse")

#You only need to do the installation step once (unless you wish to update), but you need to load the library every R session. I suggest to code this into your scripts. Typically, I would recommend to load all of your libraries at the top of your script to make it very clear what the dependencies are in your script.
library("tidyverse")

#First, you need to set your working directory.

#You can check your current working directory:
getwd()

#create a folder with two subfolders, data and R
# download the Cnidaria_BOLD_data.tsv and copy it to the /data folder
# move this script file into the /R folder

#However, set working directory relative to script file location, see Rtips on how to set this up.

#Once the path is set exactly, we can read in the data file.

#The following line uses the function read_tsv() to read in our data file and assign the data to a new data frame. (Note that I chose to give the data frame a general name, as we may want to use this script in the future for other taxonomic groups.)

#TODO Create a directory with a folder structure such that this code works for everybody.

df_bold <- read_tsv(file = "../data/06_data_cnidaria_bold.tsv")
df_bold

#We should see our data frame coming up in the global environment. There are 21094 observations (rows) and 80 variables (columns). You can click on that object name in the environment pane to bring it up and view it. Each row contains information about an individual specimen.

#We can see the variable names using the names() function:
names(df_bold)

#We can obain a summary of the data, such as what mode of data each variable contains. (Note that we can switch the data mode among compatible types for analysis, if we need to do so down the line.)
summary(df_bold)

#compare with how length() performed on a matrix
length(df_bold)

#Notice that our length is 80. Our dataframe contains 80 elements. Each element (variable) is a vector of the same length. Data frames house rectangular data. Even if there are NAs in the data, each vector (i.e. column in the data viewer) is still the same length. In this way, a data frame is a special case of a list, and indexing works the same for both. (Lists can be also be used if the data are more complex, such as containing vectors of variable lengths. Lists can also have multiple hierarchical levels... we will see examples of such lists later on.)

#Congratulations! You now have your data read in!

####4- INDEXING OF DATA FRAMES BY POSITION AND NAME----

#We may want to reduce the number of variables to create a smaller dataframe containing only the key information for analysis at this time. For example, here we are selecting 8 variables. processid is a unique identifier for each specimen. These codes are globally unique on BOLD. (Note that one specimen might have two rows, if there are different genes sequenced for one specimen. But most of the time it is one row per processid.) bin_uri is the alphanumeric code for that BIN (sequence clusters). markercode is the gene/marker name. lat is short for latitude. lon is short for longitude.

#This is an example of indexing by name. We are using the names of the variables, and this is quite a friendly way to do this.
df_bold.sub <- df_bold[, c("processid", "bin_uri", "species_name", "country", "lat", "lon", "markercode", "nucleotides")]
df_bold.sub

#let's compare to indexing by position.

#We can get an output of the position of the different variables
names(df_bold)

#Here is the same indexing task, with the positions of the variables, rather than the names:
df_bold.test <- df_bold[, c(1, 8, 22, 55, 47, 48, 70, 72)]

#Were those two approaches the same?
all.equal(df_bold.sub, df_bold.test)

#remove an interim object we were just using for that test; don't need further
rm(df_bold.test)

#Let's narrow into one column. Let's see what kind of data structure we get using single square brackets [] to subset the data frame down to one variable of interest.
test <- df_bold[, "lat"]
test
class(test)

#If we index using single square brackets, we get a data frame again. Note that "tbl" means "tibble", which is a tidyverse style of data frame. (Differences between a regular data frame and a tibble include how the data are displayed to the screen, for example. However, they typically very similar in practice.)

#Knowing we get a data frame back with single square brackets has practical implications for our analyses. Can we pass a dataframe to the mathematical function mean()?
?mean()
mean(df_bold[, "lat"], na.rm = TRUE)

#No, we cannot!

#Single square brackets are used for subsetting; this keeps the object's class.

#Now, let's compare with double square brackets, as follows:
class(df_bold[["lat"]])

#Let's try now:
mean(df_bold[["lat"]], na.rm = TRUE)

#Yes, this should work now! Answer: 30.7 is the mean latitude for specimens in this data set.

#The double square brackets allow us to extract the elements and act upon them. Each element of our data frame is a vector.

#comparing indexing by position.
mean(df_bold[[47]], na.rm = TRUE)

#Note that double square brackets are the SAME as $. Many people find the $ easier to use (as the syntax is shorter). However it is important to be aware of the long version, such as for reading example scripts online and also for writing loops. Here, again calculating the mean of latitude:
mean(df_bold$lat, na.rm = TRUE)

#Let's compare all three indexing styles for accessing a single element (variable latitude) of our data frame. Here, we are using "lat" (which houses the latitude data).
all.equal(df_bold$lat, df_bold[["lat"]], df_bold[[47]])

#I find the first style easiest for many cases.

#We can use functions that take a vector on single variables. For example, to build a histogram of latitude values:
hist(df_bold$lat)

#Recap: You have seen how to subset a data frame [] and how to select the elements [[]] you want to analyze using both position-based and name-based indexing. These are foundational skills. So, please do ensure you are comfortable with this content.

####5- EXPLORING DATA FRAME----

#In this short section, we will explore the data frame a little bit more by using a few further helpful functions.

#For example, let's see what countries are present in this dataset:
table(df_bold.sub$country)

#Sort is another helpful function. Note that we can nest functions, using parentheses. These will be executed from the inside out. Here we set decreasing to TRUE, which is an argument to the function sort(), so that we can see the countries ordered from most data to least data. (TIP: It is convenient to be able to nest functions, but I find it's good not to go overboard. We can break things up into steps (e.g. create interim objects) if we have a lot of functions we want to nest.)
sort(table(df_bold.sub$country), decreasing = TRUE)

#We may also be interested in looking at only, for example, the top 10 countries. It's easier to view on the screen quickly rather than all of the countries:
sort(table(df_bold.sub$country), decreasing = TRUE)[1:10]

#Let's create a plot showing the amount of data for the top 5 countries. (Note we won't spend time customizing the plot ... we will work on making nicer figures in ggplot2 another day!)
plot(sort(table(df_bold.sub$country), decreasing = TRUE)[1:5])

#Your turn - play with the dataframe and see what you can discover! For example, try using the table() function on variables from df_bold to see how many taxonomic classes or orders are present in this dataset.

#NEXT TIME: indexing by condition and exploring the taxonomic names in this data set.

