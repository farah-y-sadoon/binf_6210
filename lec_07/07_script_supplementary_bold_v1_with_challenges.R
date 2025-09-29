##***************************
## SOFTWARE TOOLS 2024 - CLASS 5 - SUPPLEMENTARY SCRIPT V1
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
#+ scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())

# Startup ends here

#V1 - includes coding challenges at the bottom. See V2 for example answers.

#Purpose: The purpose of this script is to provide another example of working with data from BOLD. This script provides a review of concepts we have covered so far in other scripts.

#Overview: In part 1 of this script, we will go over how to obtain data (sequence and specimen data together in one file) from BOLD, the Barcode of Life Data Systems. We will do this directly within R, without typing search terms into the web tool. In parts 2 and 3, we will next cover some examples of how to work with data stored in a data frame. Part 2 addresses indexing, while part 3 covers examples of a basic yet very useful type of plot (histogram), dives further into indexing, and covers examples of how to find data that meets a condition. Part 4 includes a few coding challenges to test and improve your skills.

#Here, we will all use the same input data file, for consistency. However, I do suggest to try out the API tool on your own. BOLD's API tools are described here:
#http://www.boldsystems.org/index.php/api_home

#For today, we will use some functions available through the tidyverse packages. We will do a more extensive introduction to the tidyverse another day. In brief, near the top of the script, we are using the function read_tsv(). That function will allow us to read in a .tsv file, i.e. a file with tab-separated values. I suggest to have a look at our .tsv file first in a text editor to see how it looks. As well, that function provides data parsing (i.e. "guessing" the type of data in each column). In general, the tidyverse file reading/parsing functions tend to perform better than the base R functions, in terms of both quality of data parsing and also speed. However, note that the base R function read.table, for example, is also commonly used by many R users. Similarly, for .csv files (comma-separated values), I would recommend generally to use the function read_csv from the tidyverse package readr, rather than the base R function read.csv(). You can read more about the readr package in the book: R For Data Science.

##If you haven't installed the tidyverse suite of packages, do so by uncommenting (i.e. removing the "#") and running the following:
#install.packages("tidyverse")

#We need to load tidyverse using library(). It is important to load the library, or else you will not be able to use the functions from those specific packages. This needs to be done EVERY new RStudio session. I always recommend to load all needed libraries near the top of your script.
library(tidyverse)

#What is the tidyverse? It is a set of R packages that work together. This system is designed to make data handling, manipulation, reformatting, clean-up, visualization, etc. more friendly and functional. You can learn more here:
#https://www.tidyverse.org/

#TIP: Remember that you can put your cursor anywhere on the line and hit Control-Enter to run that line. You can also highlight blocks of code and hit Control-Enter to run the whole block. You can also use Control-A to select all and then run an entire script at once.

######PART 1 - OBTAINING DATA FROM BOLD----

#First, we will show you how we obtained the data. We used the API tool to scrape data directly from the BOLD database.

#If you want to obtain data yourself, you would uncomment and run this line. Please do run this or similar lines on your own later and also experiment with the BOLD API tools.

#It is important for your methods to record how and when you obtained data. I ran the below line on Sept 24, 2021. You should include such information about data acquisition in your reports. Here, we are accessing publicly available data from BOLD.

Daphnia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Daphnia&format=tsv")

#Next, we can set the working directory.

#We can set the working directory for a session using the RStudio GUI.
#Session -> Set Working Directory
#You can then go to "Choose Directory". Navigate until you find the folder where you want to work. 
#You can also code in your working directory for this script into your code with setwd().

#I would recommend NOT to work out of your computer's default Downloads folder, as that folder can become very cluttered. Rather, I would recommend to set up organized folders for your courses, thesis work, personal projects, etc. For example, you might have a folder called "SoftwareTools" and then a sub-folder for Class1, Class2, Class3, etc. You would then keep your notes, data file, R scripts, and results/outputs/plots together in one folder for each class and each project.

#You can check your working directory using to be sure it is set to what you want:
getwd()

#Here is the line we can use to write the BOLD data download to hard disk. After performing your own API call, you would uncomment and run the below to gain practise in writing to hard disk. Check your folder on hard disk. Verify that the file was saved in the expected location.
#write_tsv(Daphnia, "Daphnia_BOLD_data.tsv")

#Now, we all have the save file. Read the data in from the working directory you have set, where your data file should be stored.
#You can load in the provided data file, which was uploaded to CourseLink under class 5. Here, I am reading in the saved .tsv data file from hard disk.
Daphnia <- read_tsv("../data/Daphnia_BOLD_data.tsv")

#If everything worked, you should now have an object named Daphnia in your Environment pane with 2285 observations of 80 variables. You can look at the data by clicking on "Daphnia" in the environment pane. This will open up the data in the viewer. You would need to toggle between the script editor and viewing the data. You can have a look at the variables in this data. This dataset is organized in tidy data format. Variables are columns, samples are rows, and values are in cells. For data from BOLD, each row is typically an individual biological specimen, and each row is a consistent type of sample. (Occasionally, a specimen will have more than one row, if the same specimen was sequenced for additional genetic markers.)

#Getting the data read in is an important first step - well done!

######PART 2 - WORKING WITH DATAFRAMES: INDEXING----

#In this section, we will go through some examples of indexing, whereby we pull out specific components of a dataframe. Subsetting data is a very commonly-used skill in almost any project. So, you should spend some time becoming comfortable with indexing.

#Click on Daphnia in environment pane in RStudio (generally upper-right). You should see a tab called Daphnia in the same window as your script editor. Let's explore. See that you can sort the data in the viewing pane. That is very helpful. Also, you may notice that many of the columns are sparsely populated. There are many optional data fields in BOLD that data providers can fill when they submit data to BOLD, but many fields may also be left blank. There are a few core fields that are nearly always filled in. So, data completeness could be a consideration for analyzing some variables. Note that DNA sequence data mined from GenBank into BOLD often contain fewer metadata (e.g. such records more often lack GPS coordinates, detailed collector information, etc.). BOLD also publishes data to GenBank.

#We can use a combination of numerical summaries and visualizations to explore the data and check for sample size, outliers, etc.

#First, we are checking the class of the data object.
class(Daphnia)

#This is a tidyverse-style dataframe called a "tibble". So, yes, we do have the expected data structure. Note that if we just type "Daphnia" we will get a more reasonable summary to the screen for a tibble than using a base R style dataframe.
Daphnia

#Remember that a data frame is a rectangular data object, which may contain multiple types of data. By contrast, a matrix should contain just one type of element. Data frames are a very common and useful data structure. They can be considered analogous to the familiar Excel spreadsheet. Typically, we want "tidy data", where samples are in rows, variables are in columns, and each cell contains a value (e.g. the latitude value for a particular biological specimen would be housed in one cell). Each column should contain just a single type of data.

#We want to explore our data before we start doing analysis or hypothesis testing. i.e. We want to be able to answer questions such as: What data do we have? What format are the data in? Do we have any potential errors or outliers in the data that need further examination or omission prior to analysis? We will go through a few EXAMPLES of data exploration. What would be suitable for a specific study could vary. It is important to employ critical thinking! Throughout the semester, we will see various examples of data checking and quality control; what we are looking for may vary by data type. Our brain is an important tool to keep activated.

#We can get a lot of information about data structure with this function: str()
str(Daphnia)

#Perhaps this is a little overwhelming at first, although the info is very helpful! Try to read the outputs from this for a few variables. For example, "double" is R's way of representing numerical data (accurate to 15-16 significant digits). We also have a lot of character data.

#A summary that is easier to read can be obtained using summary()
summary(Daphnia)

#This will tell us if our data are being parsed as expected. e.g. We expect nucleotide sequence data to be considered character data. We expect GPS coordinate data, i.e. latitude (column named lat) and longitude (column called long), to be parsed as numerical data. If correctly parsed, we will get a numerical summary of the distribution of the numerical data when using summary().

#If data aren't parsed correctly, we could check for formatting errors in the original data or reset the data type, if compatible.

#Getting the names and position numbers of the variables (our columns, which are the elements of our data frame).
names(Daphnia)

#For ease of working with these data further for this tutorial, we can extract a subset of variables to work with, as 80 is a lot for our purposes for now. As well, this common task of subsetting is a helpful general skill.

#So, in the below line we are creating a new tibble(tidyverse data frame) called Daphnia2. Our starting data is Daphnia. We are pulling a subset of Daphnia out and assigning it to the new object. What data are we pulling out? Specifically, we are pulling out columns 1, 8, etc., as listed below. The square brackets are for indexing; i.e. these numbers refer to specific positions. This contrasts to our usage of round parentheses to pass arguments (inputs) to functions.
Daphnia2 <- Daphnia[, c(1, 8, 10, 12, 14, 16, 18, 20, 22, 25, 26, 32, 47, 48, 55:59, 70, 72)]

#Let's go over what this line is doing in more detail. We can read from the inside out. First, we are using the c() function to create a vector. What is in this vector? The positions of the columns from Daphnia that we want to include in Daphnia2. It is helpful to work through the different components. See what this does:
column.vector <- c(1, 8, 10, 12, 14, 16, 18, 20, 22, 25, 26, 32, 47, 48, 55:59, 70, 72)

#show on screen
column.vector

#What class is this? We have created a numeric vector.
class(column.vector)

#Then, working outward, we surrounded our new vector with square brackets. Out of Daphnia, we want to pull out the specific columns we want (using position-based indexing). Then, we assigned this subset of the data to a new tibble/dataframe, to make it easier to work with downstream.

#What is the comma doing there in front of the vector of columns numbers? We didn't give a row number in front of the comma, and so we are including ALL rows (i.e. all specimens in the dataset). We are keeping all rows but only keeping a subset of the columns. We can pull out a subset of rows instead or a subset of columns and rows. Let's run through a few simpler examples and then we will return to Daphnia2 for analysis.

#What does this do? Create this object, then click on df.test in the environment pane to view the data.
df.test <- Daphnia[4, 2]
df.test
class(df.test)

#Answer: This is giving us a dataframe that actually contains just a single cell. This pulled out row 4, column 2. Don't just take my word for it. Look at the original Daphnia data in the viewing window to confirm this.

#Another example. This line pulls out row 8, column 3. Again, check for yourself by looking back at  Daphnia.
df.test2 <- Daphnia[8, 3]
df.test2

#What do you predict this will do? Make your prediction, then check.
df.test3 <- Daphnia[c(1:5), c(1:10)]
df.test3

#What do you predict this will do? Make your prediction first, then check.
df.test4 <- Daphnia[c(1:8), c(1:4)]
df.test4

#See? We can readily pull out a subset of the data that we want by position. (Later, we will see that we can also subset the data by criteria applied to the data, i.e. logical indexing also known as indexing by condition.)

#Note that we can also refer to the names of the variables, rather than using position (i.e. indexing by name). 
df.test5 <- Daphnia[, c("processid", "bin_uri")]
df.test5
dim(df.test5) #checking dimensions. should be all rows and just 2 columns.

#What will this do?
df.test6 <- Daphnia[c("processid", "bin_uri")]
df.test6

#Check if these are the same?
all.equal(df.test5, df.test6)

#Yes, but why? How can these be the same?

#Basically, if only one index position or vector is used for a dataframe, the DEFAULT is that it will be considered as columns. So, the comma isn't strictly necessary in this case, but it is better to start out with more "complete" coding practices prior to making shortcuts! So, I would recommend to you to include the leading comma for pulling out columns, i.e. using two-dimensional indexing. If we wanted to pull out only specific rows, but use all columns we would use:
df.test7 <- Daphnia[c(1:4), ]
dim(df.test7)

#What did that do? This kept the first 4 rows but all the columns... as always, check for yourself! So, here, the comma at the end is needed if we want to pull out only rows but keep all the columns.

#Please do take the time to read about indexing. I recommend that, for critical concepts, it is helpful to review the same information from different sources and formats. For example:
#http://www.hep.by/gnu/r-patched/r-lang/R-lang_45.html

#Let's tidy up our environment by clearing the interim dataframes created for our exploration of indexing.
rm(df.test, df.test2, df.test3, df.test4, df.test5, df.test6, df.test7)

#Now we are ready to continue exploring our smaller dataframe of interest, Daphnia2.

######PART 3 - WORKING WITH DATAFRAMES: BASIC PLOTS, MORE ON INDEXING, AND FINDING DATA THAT MEET A CONDITION----

#We will return to our tibble dataframe Daphnia2, which contains a selected subset of columns that were present in the original data download from BOLD. We will first create a few simple plots, introduce the $ indexing operator, and then work through a few examples of how to find data that meet a condition.

#See summary of what is in tibble dataframe Daphnia2 and the names of the variables included in this subset of our original data.
summary(Daphnia2)
names(Daphnia2)

#Plotting the data. We want to explore the relevant variables we are interested in for outliers and potential errors. Below, we will do a few example data explorations of numerical data. We will have a look at GPS data first, as location of collection is important for biodiversity research, biogeography, epidemiology, tracking invasive species, etc.

#What values do you expect for latitude? Let's formulate our expectation BEFORE we look at the data. For example, would we expect a latitude of 2000? No, the North pole is at 90 degrees (North), and the South pole is at -90 (i.e. 90 degrees South). Then, we can make a histogram plot and check whether we have unexpected values outside that range, which would likely represent data entry errors.
hist(Daphnia2$lat)

#What is this doing? Some notes about the $ sign ... This is also an indexing operator. On the left of the $ sign, we gave the dataframe name, here Daphnia2. After the $ sign, we provided the name of the variable that we wanted to look at (in this case latitude, which has the column name lat). Note that here we don't need to use quotation marks around the text of the column headers. So, $ is efficient to use for extracting a specific column in a dataframe, by column name: df$lat is equivalent to df[['lat']]. We can access that specific element and analyze it as an atomic vector.

#Let's look at the column names and numbers again:
names(Daphnia2)

#Which is the column number for latitude? 13.

#Let's further look at indexing and ways to pull out data

#Let's look at the latitude column of the dataframe Daphnia
Daphnia2$lat

#What type of data structure is this? Atomic vector of type numeric. Therefore, we can pass Daphnia2$lat to functions that take a numerical vector, such as histogram, mean, etc. hooray!

#Checking class
class(Daphnia2$lat)

#Building histogram
hist(Daphnia2$lat)

#Finding mean.
mean(Daphnia2$lat, na.rm = TRUE)

#Uncomment and run the below line to see what happens
#mean(Daphnia2$lat)

#what happened? Note we need to remove missing values here. What is the default? na.rm = FALSE. If we check the document, we can see defaults.
?mean()

#So, in summary, by using the indexing operator $ and the column name, we are getting a numerical vector in the case of latitude, which we can use in various ways we may be interested in.

#Let's contrast this with the case of indexing using square brackets. It is important to understand that, when using a single set of square brackets, we will have returned to us an object of the SAME class as the original (for tidyverse style data frame). We are looking at the latitude column as our example:
Daphnia2[13]
class(Daphnia2[13])

#See, here we are returned a dataframe! Can we pass a dataframe to a function that takes a vector?
#Does this work? Uncomment this and try to run this to see what happens.
#hist(Daphnia2[13])

#NO - that gives an error! TIP: It is very common to get errors of this type, but we can easily fix cases of such errors if we think about the data structure we need. Such errors returned to us are actually helpful, as they can keep us from executing nonsensical analyses.

#What if we want to index by column number and pass the data to a mathematical function? We instead can use DOUBLE square brackets. This allows us to dive into the data within that column, to access and analyze the numerical data in this case. Does this work? Yes, this should work.
hist(Daphnia2[[13]])

#check class between these two options:
class(Daphnia2[13])
class(Daphnia2[[13]])

#See the difference? So, we can index by either the $ indexing operator or by using square brackets. Single square brackets applied to a tibble dataframe will return us a dataframe data object, i.e. subsetting and returning us a smaller data frame. I think of double square brackets as letting us dive deeper into the data. Double square brackets will let us access the numerical data as a vector, which can be passed to functions that take a numerical vector. $ also does this! and is very helpful when we want to index by names of columns, rather than by positions. We will see more cases of when we want [] vs. [[]] later when working with lists and also when writing loops.

#Have a look at the histogram. Are there any unexpected data? Or, do latitudes lie betwen -90 and 90? Which latitudinal zone is the best represented in our dataset? Answer: This taxonomic group has been better studied for DNA barcoding in northern-hemisphere temperate and sub-Arctic environments. It is actually very common in temperate and polar areas, but sampling bias expected also to play a role in this distribution of latitude values.

#NOTE: positive latitudes are north of the equator; negative latitudes are south of the equator. BOLD uses decimal degrees.

#What range of values do you expect for longitude? (Answer: -180 to 180). Positive numbers refer to eastern longitudes, and negative numbers refer to western longitudes. What do you see when you explore the data graphically and numerically?
hist(Daphnia2$lon)
summary(Daphnia2$lon)

#Now, we can move on to asking more types of questions about the data.

#An example... A simple question we might which to ask: What is the northern-most sample? We could sort our data in the viewer. To code this, we can use:
Daphnia2$processid[which.max(Daphnia2$lat)]

#Explanation of the above line. We are starting with the lat column of dataframe Daphnia2. We are passing that vector to the function which.max(), which will tell us the index position of the maximum value. Then, we use the index position to find the processid (a unique identifier for that record on BOLD) associated with the largest latitude value.

#So, what is the northern-most sample in our dataset? We can do this in one step with the above line. Answer is: ACHAR345-18.

#check out the documentation
?which.max()

#check out what this component alone returns to us: index position. In this case, the maximum value is element 717 in that vector.
which.max(Daphnia2$lat) 

#We can also look up the actual maximum latitude value in this dataset using either of the two below methods:

#This method treats latitude as a numerical vector and asks for the maximum value using the max() function.
max(Daphnia2$lat, na.rm = TRUE)
#So, the maximum latitude is 69.1

#This approach pulls the value at index position 717 for latitude to check the above.
Daphnia2$lat[717]

#Remember that Daphnia2$lat is a vector and so we use indexing in a single dimension (not two-dimensional like for a dataframe). Again, verifying class, so that we know we have a vector.
class(Daphnia2$lat)

#Another question we might ask is: Who collected that sample? Here, we are finding the value of the variable collectors that has the maximum value of latitude. The syntax is similar to the above case.
Daphnia2$collectors[which.max(Daphnia2$lat)]

#Answer: CBG Team 1. This mean "Centre for Biodiversity Genomics" field team 1.

#What else is interesting about this sample? We can pull out that one row, row 717, and let's have a look at some of the locality information (in columns 15, 16, and 18, for example).
Daphnia2[717, c(15, 16, 18)]
#Interesting! Sample is from Cambridge Bay in Arctic Canada, the location of the Canadian High Arctic Research Station. That is a site of increasing research activity and establishment of a DNA barcoding baseline biodiversity inventory, for further monitoring.

#What is the southern-most sample?
Daphnia2$processid[which.min(Daphnia2$lat)]

#Who collected that sample?
Daphnia2$collectors[which.min(Daphnia2$lat)]

#Answer: D. Lukashanets.

#TIP: Note that you can use tab completion to save typing. This is worth our while to get used to this!

#How can we find values that meet a particular criterion?

#There are various ways to do that. Here is one example.

#How many records are collected north of the 60 degrees N?
sum(Daphnia2$lat > 60, na.rm = TRUE)

#Note that, here, we are passing a logical argument to the function sum(). "Logical true values are regarded as one, false values as zero" (from the help). So, summing the 1s will give us a count of the values that exceed 60 in this example. We can also confirm this by clicking on the object name in the environment window and then viewing the data and sorting by the lat column (to get that column in descending order would be easiest, as then the NAs will go to the bottom).

#We can check out the logical vector alone:
Daphnia2$lat > 60
#We get TRUEs and FALSEs (and NAs). TRUE is considered as a "1", and so summing the TRUEs will give us a count using the sum() function above. Logical vectors are very helpful, and we will use them a lot in various ways. It is common to want to apply a condition to our data, for exploring and analyzing data subsets.

#About now, you may be wondering why not to use the attach() command, as you may have seen that before, so that you can just use the variable names and not the full dataframe name. It is NOT advisable to do so. You may have multiple objects in your R environment, and you are less likely to create an error if you clearly specify the dataframe name and the variable name. (However, we will start to be able to avoid typing the dataframe name over and over again as we start to use more tidyverse functions and piping... stay tuned!) So, please do NOT use attach and note that you should specify the dataframe name when using the base R programming style.

#Another question we might want to answer is: How many taxonomic species names are in the dataset? i.e. What is the species richness of this dataset? So, remember we need to read from the inside out. We want to look at the column called species_name in the dataframe Daphnia2. We first use the function unique() to pull out 1 copy of each unique species name (i.e. removing duplicates). Then, we pass the results to the function called length(), which will give us the number of elements in that vector.
length(unique(Daphnia2$species_name))
#Answer: 107

#How many BINs are in the sample? i.e. What is the taxon richness of this dataset as defined by BINs (which may be considered proxies for species units in many taxonomic groups).
length(unique(Daphnia2$bin_uri))
#Answer: 137 unique BINs in this sample.

#What is the ratio of BINs to species? What does the result mean?
length(unique(Daphnia2$bin_uri)) / length(unique(Daphnia2$species_name))
#Result: 1.28. This suggests that there are more distinctive clusters, or mitochondrial lineages, than named species in this dataset. This can be an indicator of cryptic species diversity or populations that have been geographically isolated for some time.

#NOTE: To receive a BIN, a sequence must be at least 500 bp long and have < 1 % internal Ns (undetermined nucleotides). Sequences lacking a sequence, or that have a short sequence, won't have a BIN.

#How many specimen records bear a BIN (a Molecular Operational Taxonomic Unit) identifier? There are various ways to do this. One way:
sum(!is.na(Daphnia2$bin_uri))
#Answer: 1981

#Explanation: One of the easiest ways to explore data that is missing is the function is.na(). It will generate a logical vector telling us which elements within the variable bin_uri of data frame Daphnia 2 is missing. Then ! is "not". This tells us which ones are NOT missing (by flipping the TRUEs and FALSE's). The ones not missing are available! Then, remember that since TRUEs are 1s and FALSe is zero, we can sum the logical vector to get the count of records for which BIN is present (i.e. NOT missing).

#Summary: In this section we have outlined ways of working with data in a dataframe, from subsetting by column name or index position in a dataframe to seeking data that meet a specific condition (logical indexing). These are fundamental skills, so we will continue to work on these skills over the coming weeks.

######PART 4 - CODING CHALLENGES----

#Your turn... in this section in V1, you will find only commenting with instructions, no code. Please write code to solve the tasks requested. I also encourage you to play with the data in your own way. Be sure to CHECK each of your solutions to see if it works. Then, check out V2 with example answers.

#Create a new tibble called HundredDaphnia, which contains rows 1 to 100 from Daphnia.

#check the object to confirm.

#Create a new tibble called TestingIndexing, from original data Daphnia, which contains rows 50 to 100 and columns 10 to 15.

#Using HundredDaphnia for the next task...How many unique collectors collected these samples?

#Who are they?

#Now, let's see if there are any additional markers beyond COI. How many unique values of the variable markercode in Daphnia2 data frame are there? And, how are the records distributed among the markers?

#What are they?

#How are data distributed among the genes as indicated in the markercode column? What is the count for each gene?

#Create a separate tibble for the data from alternative markers (i.e. not COI-5P) and for which nucleotide data are present.

#How many unique countries are represented in  Daphnia2? 

#How are the samples distributed among countries? i.e. What is the count for each country?

#How many records are missing geographic coordinate data?

#How many records have a BIN assignment?

#Explore another element of the data that interests you!
