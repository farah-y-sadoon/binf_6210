##***************************
## SOFTWARE TOOLS 2024 - CLASS 5
##
## Karl Cottenie
##
## 2024-09-15
##
##***************************

## _ Packages used -------
library(stats)
library(tidyverse)
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

# Startup ends here

####1- INTRODUCTION----

#This script, as many of our scripts, has two main objectives: biological and programming skills

#BIOLOGICAL OBJECTIVE:

#We will use the example of the phylum Cnidaria to help us to learn how to analyze data from BOLD. Cnidaria is a group that includes jellyfish, sea anemones, and hydrozoans. When we observe the summary statistics for Cnidaria on the BOLD web portal, an interesting observation can be made: there are more species than BINs (whereas the reverse is commonly expected for many invertebrate groups due to the prevalence of undescribed or cryptic species). Are there really more species of Cnidaria than BINs (genetic clusters defined based upon the COI data alone)? Could this be due to biological factors, such a low variability in COI among closely related species? Or, could something else be causing the higher count of species? For example: perhaps Cnidaria have a lot of interim species names rather than formal Linnaean (scientific) species names? We will explore this topic in this example script.

#INFORMATICS OBJECTIVES: The above topic provides a foray into several important programming skills you will need for Assignment #1 and beyond:

#First, we will briefly review the indexing of matrices. We will also clarify that we can treat matrices as a vector, and show that we can pass a matrix with numerical data to core mathematical functions in R.

#Introduction to loading contributed packages (i.e. additional R packages that did not come with base R)

#introduction to data frame

#indexing data frame by position

#indexing data frame by variable names

#introduction to logical expression

#indexing data frame by condition

#intro to helpful functions that you can start using right away: e.g. table(), plot(), hist(), length(), unique(), sum(), is.na(), names(), dim()

#brief preview of the idea of regular expressions as a very helpful tool for analyzing string/character data (NOTE: This is NOT on Quiz #2; we will have a future lesson on regular expressions. This is a preview only.)

####2- MATRIX: INDEXING AND MATHEMATICAL FUNCTIONS----

#This section includes a review about what we learned to date about matrices and extends this by showing how we can apply mathematical functions on numerical matrices.

#Create matrix of integers 1 to 100, with 10 by 10 dimensions
y <- matrix(1:100, 10, 10)

#see documentation
?matrix

#remember, this is the same as:
y1 <- matrix(data = 1:100, nrow = 10, ncol = 10)
all.equal(y, y1)

#check class
class(y)

#Let's see what happens when we pass one index position:
y[2]
y[20]
y[50]

#Check out length
length(y)

#matrix is a vector with a dimensions attribute! That is why we get full length of the vector when we used the length() function on our matrix y.
attributes(y)
dim(y)

#Consequence: we can pass a matrix containing numerical data to core mathematical functions that take a vector. This is very useful! For example:
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
rm(y, y1)

####3- READING IN TSV FILE, INTRO TO DATA FRAME----

#If you have not done so already, uncomment (remove #) and run the following line.
#install.packages("tidyverse")

#You only need to do the installation step once (unless you wish to update), but you need to load the library every R session. I suggest to code this into your scripts. Typically, I would recommend to load all of your libraries at the top of your script to make it very clear what the dependencies are in your script.
library("tidyverse")

#Here is the line of code that I used to obtain the Cnidaria data download from BOLD, in tsv format. I like the read functions from the readr package (part of the tidyverse suite of packages) due to improved data parsing (i.e. ability to correctly "guess" data types and read in data) functionality compared to the base R file-reading functions. Please fee free to run this yourself or choose a different taxon to try. I suggest not to use a huge taxonomic group (as download will take a while).
#dfBOLD <- read_tsv(file = "http://www.boldsystems.org/index.php/API_Public/combined?taxon=Cnidaria&format=tsv")

#I then wrote to disk so that we can all start with exactly the same file
#write_tsv(dfBOLD, "Cnidaria_BOLD_data.tsv")

#Next, we will all read in the same file. What is tsv? File with tab-separated values.

#First, you need to set your working directory.

#You can check your current working directory:
getwd()

#If that's not where your file is, you would uncomment the below and set this to YOUR PATH. The below path is an example only! You need to specify your path.
#setwd("C:/Dropbox/Bioinformatics/SoftwareTools/Class4")

#You can also set your working directory using the graphical user interface (GUI) of RStudio. At the top, click, on "Session", then "Set Working Directory", then "Choose Directory". You would navigate to the folder where your data file is stored.

#Once the path is set exactly, we can read in the data file.

#The following line uses the function read_tsv() to read in our data file and assign the data to a new data frame, actually a tibble (a tidyverse-style data frame). Note that I chose to give the data frame a general name, as we may want to use this script in the future for other taxonomic groups.

# For this exercise, make sure that you have set up a folder structure such that the below line works
dfBOLD <- read_tsv(file = "../data/Cnidaria_BOLD_data.tsv")

#We should see our data frame coming up in the global environment. There are 2,4535 observations (rows) and 80 variables (columns). You can click on that object name in the environment pane to bring it up and view it. Each row contains information about an individual specimen.

#We can see the variable names using the names() function:
names(dfBOLD)

#We can obtain a summary of the data, such as what mode of data each variable contains. (Note that we can switch the data mode among compatible types for analysis, if we need to do so down the line.)
summary(dfBOLD)

#compare with how length() performed on a matrix
length(dfBOLD)

#Notice that our length is 80. Our dataframe contains 80 elements. Each element (variable) is a vector of the same length. Data frames house rectangular data. Even if there are NAs in the data, each vector (i.e. column in the data viewer) is still the same length. In this way, a data frame is a special case of a list, and indexing works the same for both. (Lists can be also be used if the data are more complex, such as containing vectors of variable lengths. Lists can also have multiple hierarchical levels... we will see examples of such lists later on.)

#Congratulations! You now have your data read in!

####4- INDEXING OF DATA FRAMES BY POSITION AND NAME----

#We may want to reduce the number of variables to create a smaller dataframe containing only the key information for analysis at this time. For example, here we are selecting 8 variables. processid is a unique identifier for each specimen. These codes are globally unique on BOLD. (Note that one specimen might have two rows, if there are different genes sequenced for one specimen. But most of the time it is one row per processid.) bin_uri is the alphanumeric code for that BIN (sequence cluster). markercode is the gene/marker name. lat is short for latitude. lon is short for longitude.

#This is an example of indexing by name. We are using the names of the variables, and this is quite a friendly way to do this.
dfBOLD.sub <- dfBOLD[, c("processid", "bin_uri", "species_name", "country", "lat", "lon", "markercode", "nucleotides")]

#Let's compare to indexing by position.

#We can get an output of the position of the different variables
names(dfBOLD)

#Here is the same indexing task, with the positions of the variables, rather than the names:
dfBOLD.test <- dfBOLD[, c(1, 8, 22, 55, 47, 48, 70, 72)]

#Were those two approaches the same? Yes!
all.equal(dfBOLD.sub, dfBOLD.test)

#remove an interim object we were just using for that test; don't need further
rm(dfBOLD.test)

#Let's narrow into one column. Let's see what kind of data structure we get using single square brackets [] to subset the data frame down to one variable of interest.
test <- dfBOLD[, "lat"]
class(test)

#If we index using single square brackets, we get a data frame again. Note that "tbl" means "tibble", which is a tidyverse style of data frame. Differences between a regular data frame (base R) and a tibble include how the data are displayed to the screen, for example. Also, when subsetting a tibble using single square brackets, we will get the same object class back. By contrast, for a base R style data frame, the default is to simplify... you may get a vector if you select one column, for example. For tibbles, as the default, the object class is not changed on us, which is nice and can make downstream analysis smoother.

#Knowing we get a data frame back with single square brackets has practical implications for our analyses. Can we pass a dataframe to the mathematical function mean()?
?mean()
mean(dfBOLD[, "lat"], na.rm = TRUE)

#No, we cannot!

#Note that na.rm() is an argument to the function mean(), and this means to remove the NAs (NA stands for "not available"). Data that are not available cannot be used in the calculation of the mean.

#Single square brackets are used for subsetting; this keeps the object's class.

#Now, let's compare with double square brackets, as follows:
class(dfBOLD[["lat"]])

#Let's try now:
mean(dfBOLD[["lat"]], na.rm = TRUE)

#Yes, this should work now! Answer: 31.78 is the mean latitude for specimens in this data set.

#The double square brackets allow us to extract the elements and act upon them. Each element of our data frame is a vector.

#comparing indexing by position.
mean(dfBOLD[[47]], na.rm = TRUE)

#Note that double square brackets are the SAME as $. Many people find the $ easier to use (as the syntax is shorter). However it is important to be aware of the long version, such as for reading example scripts online and also for writing loops. Here, again calculating the mean of latitude:
mean(dfBOLD$lat, na.rm = TRUE)

#Let's compare all three indexing styles for accessing a single element (variable latitude) of our data frame. Here, we are using "lat" (which houses the latitude data).
all.equal(dfBOLD$lat, dfBOLD[["lat"]], dfBOLD[[47]])

#Yes - they are all the same!

#I find the first style, using $, easiest for many cases, as it's easier to remember the variable names than the column number, but all are important to know.

#We can use functions that take a vector on single variables. For example, to build a histogram of latitude values:
hist(dfBOLD$lat)

#Removing unneeded test object
rm(test)

#Recap: You have seen how to subset a data frame [] and how to extract the elements [[]] you want to analyze using both position-based and name-based indexing. These are foundational skills. So, please do ensure you are comfortable with this content.

#Practise subsetting the BOLD data in a way of interest to you!

####5- EXPLORING DATA IN DATA FRAME----

#In this section, we will explore the data frame more by using a few further helpful functions.

#For example, let's see what countries are present in this dataset. We are passing the variable called "country" from tibble dfBOLD.sub to the function table(). We will obtain a count of records, by country.
table(dfBOLD.sub$country)

#Let's read about the function table()
?table
#Note after you run the above line, pick the option from the base R package (cross tabulation and table creation)

#Let's assign our table to a named object and store it in working memory so that we can explore it further.
my.table <- table(dfBOLD.sub$country)

#What class of object is this? a table, which is an array of a specific kind (table).
class(my.table)

#What are the dimensions of our table? Note that attributes() will give us more information than dim() alone, which gives the dimensions. 
dim(my.table)
length(my.table)
attributes(my.table)

#We have a one-dimensional array here (which is similar conceptually to a vector, but the dimensions attribute is assigned), because we only asked for the counts by one variable (county). There are 114 countries in our data set. (Note that oceans, e.g. Atlantic Ocean, are also stored in the "country" variable in BOLD data; so, ocean names will appear as well.)

#What names do the elements have? We see the names of the countries (or oceans). Each element (count of the records) is an integer, and each is associated with the name of the country.
names(my.table)

#Let's have a look at the first element.
my.table[1]
class(my.table[1])

#Let's pass the names in my.table to the function class(). The names are character data.
class(names(my.table))

#Let's see, on average, how many records there are per country. Pretty high! 143.7
mean(my.table)

#What about median? 36.
median(my.table)

#What does this information tell us, when we see a large difference between the mean and median? For Cnidaria, there are lot of countries with a low count of records and a few with a high count. Doing exploratory plot building and summary stats is an important step in data analysis.

#We can also look at a histogram to see the distribution of frequencies. Along the x axis we see the counts of records. Along the y-axis, we have the frequency (i.e. number of countries) with that amount of records.
hist(my.table)

#We can improve the labeling by setting optional arguments to the function hist().
?hist
hist(x = my.table, xlab = "Count of BOLD Records per Country", ylab = "Frequency (No. Countries")

#We can also adjust how the data are broken into binned sets of counts, if we wish. See the documentation for options. There are several ways you can choose to set the breaks if you don't like the default. First, I checked what is the maximum value in the data set. Then, I decided to bin in increments of 250, thus showing many countries have a count of Cnidaria records on BOLD in the bin 0-250. This is an example. Play around to see what happens.

max(my.table)

hist(x = my.table, xlab = "Count of BOLD Records per Country", ylab = "Frequency (No. Countries", breaks = c(seq(0, 3000, by = 250)))

#This is perhaps better, in showing that there are a lot of countries with low counts of records.

#TIP: Keep the base R cheat sheet handy. The seq() function and many other helpful functions are on there, as well as an overview of indexing.

#sort() is another very helpful function. Note that we can nest functions, using parentheses. These will be executed from the inside out. Here we set decreasing to TRUE, which is an argument to the function sort(), so that we can see the countries ordered from most data to least data. (TIP: It is convenient to be able to nest functions, but I find it's good not to go overboard. We can break things up into steps (e.g. create interim objects) if we have a lot of functions we want to nest, to help to cut down on errors and keep the code readable. However, nesting 2-4 functions is readable, in my opinion.)
sort(table(dfBOLD.sub$country), decreasing = TRUE)

#We may also be interested in looking at only, for example, the top 10 countries. It's easier to view on the screen quickly rather than all of the countries. We are creating a table, sorting it in decreasing order (by the counts of records) and displaying the first through 10th elements (i.e. the countries with the most BOLD records):
sort(table(dfBOLD.sub$country), decreasing = TRUE)[1:10]

#Let's create a plot showing the amount of data for the top 5 countries. (Note we won't spend time customizing the plot today ... we will work on making nicer figures in ggplot2 another day!). The plot() function is still very helpful for quick plots. It is for x, y data (i.e. where you want two axes). There are many arguments you can use for plot() if you wish.
plot(sort(table(dfBOLD.sub$country), decreasing = TRUE)[1:5])

#You may need to increase the size of your plot pane to see the plot properly.

#Your turn - play with the dataframe and see what you can discover! For example, try using the table() function on various variables from dfBOLD to see how many taxonomic classes or orders, for example, are present in this data set.

####6- INDEXING DATA FRAME BY CONDITION----

#Interesting observation regarding the Cnidaria data summary from BOLD: there were more species than BINs. Why? Before looking for a biological explanation, let's explore potential technical explanations. For example, did a lot of the records have no sequence? Or, are there a lot of interim species names, which could be causing the species count to look artificially high?

#unique() is helpful function. Here, this will give us all of the unique BIN codes. This will return to us a vector of the unique values of the bin_uri column of dfBOLD.sub (i.e. removing duplicates).
unique(dfBOLD.sub$bin_uri)

#length() is also a very helpful function, which I use regularly. How many elements do we have? Reminder: functions will be run from inside out, starting with the most nested pair of parentheses. So, this is saying how MANY unique BINs there are.
length(unique(dfBOLD.sub$bin_uri))
#1396.

#How many unique species are there?
length(unique(dfBOLD.sub$species_name))
#3253. But is this count influenced by interim names?

#Note we can create a new variable on the end of our dataframe using the assignment operator. Here, we are using a function from the stringr package (part of tidyverse) to have a look at the character data and to count specific characters. This is a BRIEF preview of regular expressions (regex). We will return to this in a future class (this isn't on Quiz #2!). Here, we will count the number of spaces, dots, and digits (i.e. 0-9) in each species names. Do we have unusual or interim species names? The expected number of spaces in a standard scientific (Linnaean) species name is 1.
dfBOLD.sub$spaces <- str_count(string = dfBOLD.sub$species_name, pattern = "[\\s\\.\\d]")

#Have a look at our data frame in the viewer. We should have a new variable, called spaces, added onto the end of our dataframe, containing the counts of spaces, dots, and digits in the species names for each BOLD record of Cnidaria.

#Let's explore the distribution of our data. How many of the records have a species name containing just one space/dot/digit, and how many have more than that?
hist(dfBOLD.sub$spaces)

#Next, we want to apply a condition to our data. In a minute, we will want to keep and examine the records (rows) with exactly 1 space in the species name (as an indicator of a formal, scientific species name).

#== means "exactly equal to", and here produces a logical vector. TRUE means that the spaces variable was exactly equal to 1 for that record, and FALSE means that anything OTHER THAN 1 was present in the spaces variable for that record. NA means that there was NA in the species colunn (check for yourself). Let's have a look...
dfBOLD.sub$spaces == 1

#logical vector also is numeric, as FALSE means 0 and TRUE means 1. Therefore, we can use logical vectors with mathematical functions! cool, eh? So, by summing 1s (TRUEs), we can readily obtain a count of the BOLD records in our dataset (i.e. rows in our dataframe) that met out condition. We are removing NAs (data not available).
sum(dfBOLD.sub$spaces == 1, na.rm = TRUE)
#16772

#Now that we've seen the logical vector alone, let's try indexing by condition. By placing a logical expression to be evaluated in the rows indexing position, a logical vector will be created (TRUE/FALSE or 0/1). As a consequence, we will KEEP all rows meeting that condition. Placing nothing after the comma means we are keeping all columns.
df3 <- dfBOLD.sub[dfBOLD.sub$spaces == 1, ]

#Let's check the properties. We should have a smaller dataframe that contains the original number of columns but fewer rows. Yes.
dim(df3)

#Let's also check if we got rid of all rows having something other than 1 in the "spaces" variable.
unique(df3$spaces)

#Almost! We now just have two values in the spaces variable: NAs and 1s. For further analysis, we only want the rows where the species name has only 1 space and is NOT missing.

#Therefore, we can add MORE conditions. For today, we will perform this in classical style using & (logical operator). Using & means it is used as vectorized operator. In addition to keeping those rows for which column spaces is equal to 1, we are also keeping ONLY those rows for which species_name and bin_uri is NOT missing (i.e. is present). Missing data for these variables will not help us with our question, so we are removing them. Note that ! means "not". While it seems like a double negative, in R !is.na() is actually one of the easiest ways to get data that are present!
df3 <- dfBOLD.sub[dfBOLD.sub$spaces == 1 & !is.na(dfBOLD.sub$spaces) & !is.na(dfBOLD.sub$species_name) & !is.na(dfBOLD.sub$bin_uri), ]

#Again: we are KEEPING rows where there is exactly 1 in the spaces column and for which spaces, species_name, and bin_uri are all NOT missing (i.e. there is data there).

#Comment: Isn't it getting annoying to type the dataframe name over and over again? We will see an easier way to work with the variables in data frames once we start using more tidyverse functions and piping next week. However, it is important to understand base R style syntax. For today, an alternative approach that is shorter is to use the subset() function (thanks Jacqueline):
df3_alt1 <- subset(dfBOLD.sub, spaces == 1 & is.na(spaces) == F & is.na(species_name) == F & is.na(bin_uri) == F)

all.equal(df3, df3_alt1)
rm(df3_alt1)

#Now, let's compare now that we've filtered the dataset to the more useful portion for our question.

#The only value remaining in the spaces variable should be 1. Yes.
unique(df3$spaces)

#How many unique species names remain? 1663.
length(unique(df3$species_name))

#How many unique BIN names remain? 912.
length(unique(df3$bin_uri))

#We can see that there remains a difference between the number of BINs and the number of species names in Cnidaria. So, this pattern is not only driven by the interim species names. Let's discuss why! This could possibly be caused by a slower rate of molecular evolution than in many other animal groups. Other interesting possibilities could include hybridization. Another possibility could be phenotypic plasticity, whereby the same genotype can produce different phenotypes, which could be described under different taxonomic names, under different environmental conditions. From the literature, slower rate of molecular evolution has been put forward as an explanation for lower resolution of DNA barcode sequence data in Cnidaria.  However, further investigation of the literature, and possibly other data types, would be needed to draw more substantial conclusions.

#We have seen a variety of indexing tools and helpful functions to let us explore the data in our data frame.

#YOUR TURN - Play around with the code in this tutorial to explore the data and functions on your own. Then, you can try it out on your own data for your Assignment 1! You can go ahead and download data using the API tool (by changing out "Cnidaria" from the example provided above). Have fun!

