##***************************
## Software Tools - Class 6 - Intro to Biodiversity Analysis Using tidyverse
##
## Karl Cottenie
##
## 2024-09-20
##
##***************************

## _ Packages used -------
library(tidyverse)
conflicted::conflict_prefer("filter", "dplyr")
library(viridis)
# + scale_color/fill_viridis(discrete = T/F)
theme_set(theme_light())
library("vegan")


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

#Version 2 - includes example answers to the coding challenges.

####1 - INTRODUCTION----

#This script has two main sets of objectives:

#Biological objectives - We will conduct analysis of biodiversity using BOLD data. Specifically, we will investigate sampling completeness through building species or BIN accumulation curves. Additionally, our data will be ready for further analysis using the R package vegan, such as for calculating community diversity and composition metrics. For example, we may wish to calculate the dissimilarity of BIN composition among different geographic regions, using the vegdist() function.

#Informatics objectives - This script provides an introduction to:

#a) using piping to perform multiple steps in sequence and using tidyverse functions to summarize data in a convenient fashion that is both easy to code and easy to read

#b) using function filter() from dplyr package to filter data set easily by a condition

#c) using the pivot_wider() function from the tidyr package to reshape our data set for analysis in vegan package

#NOTE: Filtering, summarizing, and reshaping data are core skills that you will use widely beyond this specific example using BOLD data. Therefore, I recommend to take the time to become comfortable with these functions.

####2 - LOAD PACKAGES AND CHECK CONFLICTS----

#Uncomment and install package, if needed.
#install.packages("vegan")
#Then, load the library so that you will have access to the functions by running the following line. The vegan package provides substantial functionality for ecological analysis.
library("vegan")

#uncomment and install package suite, if needed. Load package. The tidyverse suite of packages provides diverse tools for reading and writing data, reshaping data sets, filtering data, facilitating analysis, and visualizing data.
#install.packages("tidyverse")
library("tidyverse")

#Notice "conflicts" when you load tidyverse. These are cases where there is already a function with that name loaded. These functions from the base R set of packages, including basic stats package, are masked when we load tidyverse. However, this functionality is not lost to us. Rather, if we want to use one of those functions, we would need to use the full function name, which includes the package name at the start. For example: stats::filter(). If we only use filter() alone, we will be accessing the dplyr package version of that function name. Actually, writing full function names is good practice and is especially helpful if you want to write a fully reproducible script that can be run correctly by anyone, regardless of which packages they have loaded. In practice, for our own analysis, it is more common simply to use the function name and be aware that these conflicts can exist (and can be easily remedied by using the full function name when needed). Using full function names is especially important if we are using lots of packages in a more complex script.

####3 - READ DATA----

#We will continue with analyzing the Daphnia data (which we analyzed in the class 5 supplementary script). This is the line I used to obtain the Daphnia data from BOLD. Data downloaded Sept. 24, 2021.
#Daphnia <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Daphnia&format=tsv")

#This is the line I used to write the data to disk, so that we would all have the same data.
#write_tsv(Daphnia, "Daphnia_BOLD_data.tsv")

#Remember to set your working directory.
#setwd()
#Check working directory
#getwd()

#Load data. We will continue to work with our Daphnia data from BOLD, for now.
dfBOLD <- read_tsv("../data/Daphnia_BOLD_data.tsv")

class(dfBOLD)

dim(dfBOLD)

#Note that I am using a more generic data frame name (rather than Daphnia), as we may want to adapt the below code for other taxonomic groups. That way, we can change our input data but don't have to change the data frame name throughout. It would be confusing to be analyzing, for example, echinoderm data under a data frame named Daphnia. Therefore, for code reusability, I suggest to consider using more general object names.

#If this is successfully read in, you should see a dataframe in your global environment containing 2285 observations and 80 variables.

####4 - PIPING, SUMMARIZING DATA BY GROUPS, AND FILTERING BY CONDITION----

#In this section, we will briefly introduce R programming using piping. Piping allows us to perform multiple actions in sequence, without cluttering our environment by creating many interim data objects. When using piping, the OUTPUT of one step serves as the INPUT for the next step. The code can also be easier to read, avoiding long sets of nested parentheses.

#If you are working on an analysis involving a lot of steps, it can be helpful to build up your code (using small datasets to check your work along the way) and then use piping to reduce the creation of interim objects that clutter your environment. Also, piping allows a logical flow of steps. However, for a more clear flow and for manual verification, it can still sometimes be helpful to create interim objects. So, don't be afraid to go ahead and create an interim object (e.g. a subset of your original data) if you think it will help you.

#Which countries have the most barcode data? Here is just the code alone. Note you can put your cursor anywhere in this block and hit Control-Enter and it will all run.
dfBOLD %>%
  count(country, sort = TRUE)

#Explanation: What is going on here? First, we specify what tibble (a tidyverse-style data frame is called a tibble) we are working with (dfBOLD). Then we put the pipe (%>%) at the end of that line. That means that the output from the first line is the input to the next line. One GREAT thing about piping is that we can now just use the variable names! We don't have to type the name of the data frame over and over again. Next, we are passing the variable country to the function count(). What count() will do for us is to group the rows by country and count the number of rows (observations) within each country. I added the argument sort = TRUE due to personal preference. I like to see such data sorted.

#Now, let's compare against what we have seen so far, to generate such information using base R functions:
sort(table(dfBOLD$country), decreasing = TRUE)

#Very similar! Some differences: with piping, we can avoid having multiple sets of nested parentheses. Also, the tidyverse function gives us explicit information about the NAs in our data set, which is actually really helpful! What is the most common country? NA. (Those are likely records mined from GenBank, which are lacking geographic information in BOLD. Data submitted directly to BOLD, by researchers engaged in DNA barcoding research, often have geographic metadata.)

#Can we calculate other summary statistics for groups? Yes! The function summarize() is very helpful for that.

#For example, perhaps we wish to find the mean latitude for each province and state. (This is an example, but why might anyone wish to do such a task? Perhaps we want to create a global-scale map about sample sizes. We might wish to calculate the mean latitude (and, afterwards, longitude) for the data in this dataset for each province or state, so that we can impute missing coordinates by using the mean for this data set for provinces and states.)

#First, base R style example solution ...

#Using aggregate() function. What are we aggregating? latitude data (argument x). What are we grouping by? province_state (we need to make a list for the by argument). What function do we want to apply to the groups? mean, with na.rm set to TRUE. We are creating a new dataframe, dfAvgLat1.
dfAvgLat1 <- aggregate(x = dfBOLD$lat, by = list(ProvinceState = dfBOLD$province_state), FUN = function(x) mean(x, na.rm = TRUE))

summary(dfAvgLat1)

#If we wanted to use simple means, we can just pass mean() alone to the FUN argument. However, if we do a summary, we will see that we get fewers NAs in the outcome if we remove NAs while taking the means. This is because NA's will be dropped, enabling the calculation of the mean all of the provinces and states that have any numerical data available in the latitude column. You can uncomment and run the below lines and see the difference.
#dfAvgLat.test <- aggregate(x = dfBOLD$lat, by = list(ProvinceState = dfBOLD$province_state), FUN = mean)
#summary(dfAvgLat1)
#summary(dfAvgLat.test)
#rm(dfAvgLat.test)

#We may want to order the data frame. We can do this, for example, by using the order() function. Easier to do this in separate step.
(dfAvgLat1 <- dfAvgLat1[order(dfAvgLat1$x), ])

#Now, have a look at the data frame dfAvgLat1 in the viewer.

#Now, let's try this again using dplyr functions. Here is the code all at once. We will then go through line by line and explain.

dfAvgLat2 <- dfBOLD %>%
  group_by(province_state) %>%
  summarize(avgLat = mean(lat, na.rm = TRUE)) %>%
  arrange(desc(avgLat)) %>%
  print()

#The same code is repeated below, with commenting. Let's walk through this line by line:

#dfAvgLat2 <- dfBOLD %>%
#Name of data frame we want to end with, then assignment operator, then name of data frame we will be manipulating. Line ends with pipe, %>%, telling R there is more. i.e. dfBOLD is our input to line 2.

#group_by(province_state) %>%
#Here, we are using the function group_by() to group the data by province_state. Note we don't have to specify the data frame name! That was already specified at the beginning in line one, so we don't have to type it again. So, using group_by(), we are grouping rows of data by provinces and states, e.g. grouping all records from Ontario together, etc. I find that group_by() is one of the most helpful functions in existence in R. We end the line with the pipe, indicating that the grouped data are the input to the next line.

#summarize(avgLat = mean(lat, na.rm = TRUE)) %>%
#The function summarize() is a general function. We can summarize the data in many ways. Here, we are creating a new variable called avgLat. We are setting this variable to be the mean of latitude. The mean will be calculated for each GROUP (i.e. each province and state). We can add the argument na.rm = TRUE to pass to the function mean. This means that it will be possible to calculate the average latitude for any provinces or states that have any data, even if there are NAs for some records.

#arrange(desc(avgLat)) %>%
#We are sortng the results of averaging latitude in descending order. 

#print()
#ending with a print statement, if we want to print to screen right away.

#Also, note than when you print a tibble (a type of data frame) to the screen, only a small amount of the data prints to screen, unless you specify otherwise. This can make larger data easier to look at on the screen. Many R users prefer to use tidyverse structures and functions, but there is a lot of variability in personal preference.

#We can also look at dfAvgLat2 in the viewer in RStudio.

#Now, what if we also want to apply a condition to the data first? For example, maybe we first want to filter out records that are lacking province_state information. Then, perhaps we are only interested in sub-Arctic and Arctic records for our project. We can easily apply filtering steps first. For example:
dfAvgLatArctic <- dfBOLD %>%
  filter(!is.na(province_state)) %>%
  filter(lat >= 55) %>%
  group_by(province_state) %>%
  summarize(avgLat = mean(lat, na.rm = TRUE)) %>%
  arrange(desc(avgLat)) %>%
  print()

#Here, we have added two filtering steps. First, we retain only records for which province_state is NOT NA. Next, we filter to retain records for which latitude is equal to or greater than 55. Note, that we use filter() to specify which records to RETAIN. So, we need to structure our code accordingly. Try removing the not operator ! and see what happens. Not much! Because we would be only keep rows for which province_state is NA. These are EXAMPLES. For your own analyses, you would need to think carefully about what records you want to retain or filter out to answer your question.

#Want to compare with base R style indexing by condition? We could do this in one step, but it's easier to do in multiple steps:

#First, we are indexing by condition to filter the data set. We are applying two conditions to our rows. We want to keep rows for which province_state is NOT missing, and we are keeping rows with latitude at or north of 55 degrees. There is nothing after the comma, as we are keeping all columns. (Please note we could add more conditions if we want, such as if you want to keep only rows that have a numerical value present for latitude.)
dfBOLDfiltered <- dfBOLD[!is.na(dfBOLD$province_state) & dfBOLD$lat >= 55, ]

#The rest is the same as further above, except we are using the filtered data frame.
dfAvgLat3 <- aggregate(x = dfBOLDfiltered$lat, by = list(ProvinceState = dfBOLDfiltered$province_state), FUN = function(x) mean(x, na.rm = TRUE))

#Ordering our data frame by latitude. And, by putting parentheses around the whole thing, this means that we will print to screen right away.
(dfAvgLat3 <- dfAvgLat3[order(dfAvgLat3$x), ])

#Play around with the code above to understand what the different components do. Summarize the data in a different way, of interest to you.

#Which programming style do you like better?

#Let's finish off this section by removing objects we don't need any more downstream in this script.
rm(dfBOLDfiltered, dfAvgLat1, dfAvgLat2, dfAvgLat3, dfAvgLatArctic, dfAvgLat.test)

####5 - REPHAPING DATA AND INTRODUCTION TO BIODIVERSITY ANALYSIS USING VEGAN----

#For this part, we will use the BINs as a molecular proxy for species to explore the completeness of sampling in the genus Daphnia for the DNA barcoding campaign. Barcode Index Numbers (BINs) are unique identifiers assigned on the basis of clustering patterns in the DNA barcode sequence data. The BIN algorithm was originally calibrated against traditional species identifications for selected, well-studied animal groups. The BIN algorithm is also run for sequences lacking a traditional taxonomic identification. BINs are very helpful for biodiversity analysis, especially of invertebrates, where obtaining morphological species-level identifications is very difficult or sometimes impossible, due to the high diversity, many undescribed species, and sometimes cryptic species (evolutionarily distinct species that look the same to humans).

#We will use several functions from the vegan package. For example, the function specaccum() will build a species accumulation curve. Remember we will use BINs as a proxy for species. In this example, we want to ask: how well sampled is this taxon? Note that this function requires a community object (comm) to act upon. So, first we need to get our data ready.

#We will be getting the data ready for downstream analysis in the package vegan. We need an object of the type comm, as can be seen by looking at the documentation for specaccum() and also in some online tutorials about this package. That object needs to have species (here, BINs as proxies for species) as column headers and then counts as the data. Sites are in rows. First, we are analyzing all the data together as one big site (the whole world).
?specaccum()

#Here, we are grouping the data by BIN and counting the number of records in each BIN. We will walk through this together. We are creating a new tibble called dfCount.by.BIN, by taking the data from dfBOLD, grouping by BIN (a proxy for species), and counting the number of specimen records (which are rows in our BOLD data) per BIN.
dfCount.by.BIN <- dfBOLD %>%
  group_by(bin_uri) %>%
  count(bin_uri)

#Note that count() will ungroup the groups afterwards. See documentation for count(). Here, this is good for what we want. However, we might frequently use summarize() instead, even for counts, in the case we want to keep working with the groups further down in a set of piped expressions

#Have a look at dfCount.by.BIN in the viewer and be sure that you understand what the above lines of code did.

#Now, we are using a very useful function called pivot_wider() from the tidyr package (within tidyverse suite of packages) to get the data into the community data object format. (Note that this is similar to an older function called spread(), which is no longer under development, and so usage of pivot_wider() is now recommended.)
?pivot_wider()
dfBINs.spread <- pivot_wider(data = dfCount.by.BIN, names_from  = bin_uri, values_from = n)

#Explanation: We are creating a new data frame called dfBINs.spread. We are using the function pivot_wider() to reshape our data. The data argument is dfCount.by.BIN. Those are the counts of records, by BIN. Have a look at that object. All the data are there that we need to make our community object. For the names_from argument, we are specifying what we want our new columns to be called. We want our new columns to be the BINs. So, what does key do? names_from argument will result in pulling pulling out each VALUE from the bin_uri column (i.e. each unique BIN code) and put those as column headers in our new data frame. For the values_from argument, we are specifying what we want the data to be in the body of the new data frame. We are specifying "n", which are the counts in the dfCount.by.Bin data frame.

#What does this look like when we look at dfBINs.spread? Now, rows are sites (here, we are treating all samples as belonging to one big site), columns are BIN identifiers (i.e. a different column for each BIN), and values in cells are the counts of individuals per BIN. So, this is a different format for our data, which is needed for biodiversity and community ecology analysis using the vegan package.

#I suggest to consult the data wrangling cheat sheet and to keep that handy to help you (noting that pivot_wider() and pivot_longer() are the new function names).

#Now, let's do something with the data! There are many options. The vegan package alone has a manual that is nearly 300 pages long! If you are interested in community ecological analysis and biodiversity, this is an invaluable package for you to explore further. For this course, we will just explore a few functions, but please do explore further on your own outside of class towards your Assignment #1.

#Let's try out the function rarecurve() to build a rarefaction curve. First, see documentation:
?rarecurve()

x <- rarecurve(dfBINs.spread, xlab = "Individuals Barcoded", ylab = "BIN Richness")

#Have a look at the plot. You can click on "Zoom" if you want to get a better view.

#Have a look at the documentation and the defaults to understand what this function is doing. This function is randomly sampling the individuals available in the total data set, at various levels of completeness (shown on the x-axis). The default is to step in increments of 1. So, at sample size 1, one individual is chosen. At sample size 2, two individuals are chosen, etc. Shown on the y-axis is the number of BINs that were detected at a given level of sampling.

#And, what do you notice about the shape of the curve? the rate of discovery of BINs is slowing with increasing sampling size, but the curve is still increasing. We would predict that we would continue to add BINs to the database with further sampling.

#Now, we will consider the country data. The question we are addressing through the below analysis is: As we add countries, do we add a lot of new BINs? i.e. Do unique BINs tend to be added as we sample additional countries? First, as often the case, we need to get the data into the correct format. First, we are counting the number of specimen records per BIN per country. Note that we can pass two grouping variables to group_by(), as we want a count of the number of specimens per BIN per country.
dfBINs.by.country <- dfBOLD %>%
  group_by(country, bin_uri) %>%
  count(bin_uri)

#First, make a prediction about what this did. Then, have a look at the resulting object to check your prediction and understanding of the syntax here.

#Next, we can remove the rows where country is NA, as such records do not contribute to answering our question. We do not know the country for those records. We also remove rows where the BIN is NA, as such records also do not contribute to answering our specific question in this case. (You may ask: why are such records present in this dataset? There can be specimen records on BOLD where sequencing hasn't been performed yet, where sequencing failed, or where the sequence was too short to be assigned to a BIN. Sequences must be at least 500 base pairs long and contain fewer than 1% Ns (undetermined nucleotides) to receive a BIN assignment.)

#Example using base R style indexing by condition:
dfBINs.by.country.na.rm <- dfBINs.by.country[!is.na(dfBINs.by.country$country) & !is.na(dfBINs.by.country$bin_uri), ]

#Have a look at this carefully to see the syntax. This is an example of combining two logical conditions. We want to retain rows having country present (i.e. country is not missing) and also to have bin_uri present (i.e. bin_uri is not missing). The & (and) and | (or) logical operators are applied element-wise to vectors. So, here, with &, both conditions will be applied. 

#You can read more about operators in R, including logical operators, here:
#https://www.datamentor.io/r-programming/operator/
#https://cran.r-project.org/doc/manuals/R-intro.html#Loops-and-conditional-execution

#If we wish, we could instead do this filtering of missing data using tidyverse functions. Which do you find nicer?
dfBINs.by.country.na.rm2 <- dfBINs.by.country %>%
  filter(!is.na(country)) %>%
  filter(!is.na(bin_uri))

#checking that these are the same - yes! So, can remove one.
all.equal(dfBINs.by.country.na.rm, dfBINs.by.country.na.rm2)

rm(dfBINs.by.country.na.rm2)

#Remember that we don't have to give the name of the data frame more than once when using piping. The filter() function retains rows matching our condition. So, we can very easily pipe together multiple filtering criteria in a logical flow that can be easier to code and to read. yeah!

#By contrast, the select() function is used to keep the mentioned columns (variables) and is another helpful function to be aware of. See the dplyr cheat sheet.

#Then, again, we need to spread the data so that we create a comm object (dataframe of a specific data type and formatting). Have a look at that object.
dfBINs.spread.by.country <- pivot_wider(data = dfBINs.by.country.na.rm, names_from = bin_uri, values_from = n)

#check class.
class(dfBINs.spread.by.country)

#Have a look at dfBINs.spread.by.country to be sure you follow what the above line did. Our first argument (data) to pivot_wider() was our dfBINs.by.country data.na.rm (have a look at that object to see its structure). Our second argument (names_from) is the column we want to spread by and for naming the colunns. i.e. We want to take our column bin_uri and turn each unique BIN code in that column into its own column header in our new data frame. Our third argument (values_from) is what we want the data to be. We want the cells to contain counts.

#I would also recommend to read Chapters 5 and 12 in "R for Data Science" for more background on tidy data, data transformations, and data reshaping (e.g. spreading and gathering).

#Next, we are converting the NAs to zeroes, as required for downstream analysis. Our zeroes in dfBINs.spread.by.country represent that a given BIN hasn't been barcoded in that country (at least not yet!).
dfBINs.spread.by.country[is.na(dfBINs.spread.by.country)] <- 0

#Have a look again at that data frame in the viewer. Our NAs have been converted to zero.

#check out documentation for function
?specaccum()

#In the following section, I am deliberately leaving a line of code that generates an error. You can uncomment the next line and try to run it.
#specaccum(dfBINs.spread.by.country)
#What error message do you notice? Why? Tip: Rather than become frustrated, it is helpful to take the time to READ error messages. Upon reading the message, we realize the function is looking for numeric data.

#We can investigate the structure of our data using str().
str(dfBINs.spread.by.country)

#And, summary() will generate an output summary that some users may find more readable.
summary(dfBINs.spread.by.country)
#See... the first variable is of type character. We need our data all to be in numerical format to have a correct comm object type to pass to the function specaccum.

#We can see that our first variable contains the country information (which is what we would want for some types of analyses, tidy data format with all data in columns, but not for the specaccum() function in vegan).
names(dfBINs.spread.by.country)

#Now, to get the data ready for vegan, we are going to do something that is typically AGAINST tidyverse principles. Typically, we want all of our data to be in the columns and cells, not as a names attribute of the rows. However, we need to go with the data format wanted by vegan, if we want to access the powerful set of functions available through that package.

#The below line is setting the rownames as country, rather than having country as a data column. Note that we don't want to lose the country information.
dfBINs.spread.by.country <- dfBINs.spread.by.country %>%
  remove_rownames %>%
  column_to_rownames(var = "country")

#Have a look at dfBINs.spread.by.country before and after. Before, country values were in their own column. After, country is to the left margin as row names (not data inside the cells). This is sort of analogous to putting a label on the outside of a parcel, rather than as contents inside the parcel.

#Note that tidyverse wouldn't like that we set the row names as special row attributes. In older versions of the tidyverse packages, we would receive a warning message that "Setting row names on a tibble is deprecated." what does "deprecated" mean? This means that the tidyverse disapproves! In the tidyverse style, all data should be in data columns. However, here, we want to use the R package vegan, as it is very useful and has many useful functions for ecology and biodiveristy analysis. And, we need to have the data in the right format for that package. We do want to keep the country information with the rows, and so we set the row names as special attributes. Tidyverse says right on its homepage that its authors are opinionated about how things should be done! There are many useful things about the tidyverse, but we will frequently need to step outside of that framework for bioinformatics. For bioinformatics, we will need to learn tools to deal with "non-tidy" data, including non-rectangular data and complex strings.

#Now we can run a species accumulation curve analysis. Resampling of sites is performed to see how BINs accumulate as sites (in this case countries are treated as sites) are added.
AccumCurve <- specaccum(dfBINs.spread.by.country)

#We can make a plot of the model. click on "Zoom" to get a better view of the plot.
plot(AccumCurve, xlab = "Countries Sampled", ylab = "BIN Richness")

#You can see the help for the function specaccum() for more information. You can also customize the plot (e.g. add different axis labels, change the colours, etc.)

#What is that plot telling us? As countries are sampled, additional BINs are added to the dataset. So, not all BINs of Daphnia are found everywhere. This plot gives us a somewhat different answer compared to the individual-based rarefaction curve. The shape of the individual-based curve is indicating that some regions have likely been sampled quite a bit, and there are few new BINs to be sampled, as the curve is leveling off. However, the accumulation curve by country is telling us that as we sample different countries, new BINs are being encountered at a steeper rate. So, from a global perspective, there are likely new BINs to be sampled as additional parts of the world are sampled.

#Note that getting the data into the correct format was half the trouble! Probably more than half the trouble! That is often the case. Outside class, you can now explore more vegan functions now that the data are formatted. Below, you will find several coding challenges to help you develop your ability to work with data in R. Have fun and see you next time!

#Let's clear out objects not needed any more.
rm(AccumCurve, dfBINs.by.country, dfBINs.by.country.na.rm, dfBINs.spread, dfBINs.spread.by.country, x, dfCount.by.BIN)

####6 - CODING CHALLENGES WITH EXAMPLE ANSWERS----

#Write code to solve these tasks. Also, we encourage you to manipulate and play with the data in your own way.

#CHALLENGE #1: Start with the original data called dfBOLD. Create a subset as a new object, retaining only records found in the northern AND western hemispheres. Let's do this task based upon GPS coordinates only. We could also recover records that meet this criterion on the basis of country names, if we wanted to go further with a particular analysis, but that is beyond our intended scope here. Then, write lines of code to check that your solution works as intended.

#Potential answer in base R indexing style. It is helpful to remove NAs as well. In general, careful attention is needed as to how we treat NAs. Recall that west longitudes are negative using decimal degree GPS coordinates format, and east latitudes are positive. North latitudes are positive, and south latitudes are negative. Example solution:
dfBOLD.northwest <- dfBOLD[dfBOLD$lat > 0 & dfBOLD$lon < 0 & !is.na(dfBOLD$lat) & !is.na(dfBOLD$lon), ]

#The following lines of code check that we did this task correctly. It is always wise to check! Anyone can make a typing error. So, we recommend to get into good practices always to check your work at EVERY step.

#First, checking dimensions. We should NOT have lost columns. column number should be the same in original and Daphnia.northwest datasets. Yes, both have 80 variables (columns).
dim(dfBOLD)
dim(dfBOLD.northwest)

#Next, checking distribution of latitude and longitude values. Yes, the min and max for each is as expected. Latitudes should be positive and between 0 and 90, and longitudes should be negative and between 0 and -180.
summary(dfBOLD.northwest$lat)
summary(dfBOLD.northwest$lon)

#Checking for NAs. No NAs in lat and lon.
sum(is.na(dfBOLD.northwest$lat))
sum(is.na(dfBOLD.northwest$lon))

#CHALLENGE #2: Try to code the above task using both base R style and tidyverse style (using piping). Soon, you will see which you prefer and can choose just to use one. However, it is to your benefit to be able to read and understand code written in both styles.

#Second solution with tidyverse style.
dfBOLD.northwest2 <- dfBOLD %>%
  filter(lat > 0) %>%
  filter(lon < 0)

dim(dfBOLD.northwest2)

#Can check other properties too, as above.

#CHALLENGE #3: How many unique BINs (bin_uri column) are found in Argentina in this dataset, dfBOLD? And, what are those BINs?

#One possible base R solution. answer = 20. Note I am omitting rows for which either country or bin_uri is NA, as such rows don't contribute to answering our question.

length(unique(dfBOLD$bin_uri[dfBOLD$country == "Argentina" & !is.na(dfBOLD$country) & !is.na(dfBOLD$bin_uri)]))

#what are they?
unique(dfBOLD$bin_uri[dfBOLD$country == "Argentina" & !is.na(dfBOLD$country) & !is.na(dfBOLD$bin_uri)])

#Let's try again with tidyverse style. So, first we filter for rows where country is exactly equal to Argentina and where BIN is not missing. Then, we summarize the data. We create a variable called count.unique.bins and we assign to count.unique.bins the length of the vector of unique BINs. Answer = 20.
dfBOLD %>%
  filter(country == "Argentina") %>%
  filter(!is.na(bin_uri)) %>%
  summarize(count.unique.bins = length(unique(bin_uri)))

#CHALLENGE #4: How many records in this data set are missing country information? TIPS: Remember that sum() is helpful here, because it will count the number of TRUEs (i.e. 1s), if we pass a logical vector to sum(). So, the is.na() function is creating a logical vector, which can be passed to sum(). I suggest always to check the components for yourself if you are still getting used to this structure! Answer: 872
sum(is.na(dfBOLD$country))

#CHALLENGE #5: How many records are missing a BIN assignment? Answer: 304.
sum(is.na(dfBOLD$bin_uri))

#CHALLENGE #6: Let's prepare an overall summary of the data by country. Let's summarize the number of unique BINs per country.

#Example solution: the below will output the results to the screen. Note I am first filtering the dataset to remove rows from consideration that lack country or BIN data. (Also, remember, we could assign the output to a data object, if we wanted to work with the output further downstream. In this example solution, we are just printing to screen.)
dfBOLD %>%
  filter(!is.na(country)) %>%
  filter(!is.na(bin_uri)) %>%
  group_by(country) %>%
  summarize(count.unique.bins = length(unique(bin_uri))) %>%
  arrange(desc(count.unique.bins)) %>%
  print()

#CHALLENGE #7: What proportion of records are missing country information among records mined from GenBank vs other records? (Hint: see institution_storing column)

#Example answer #1 using base R indexing and functions: As usual, there are multiple ways to solve this! Here is one way. I am first creating a subset for records from GenBank, which I will use further below to obtain the answer to the question posed. See the documentation for the function grep(). I am using the function grep() to search for a specific string in the institution_storing variable.
dfBOLD.GenBank <- dfBOLD[c(grep("GenBank", dfBOLD$institution_storing)), ]

#Let's check our work. We can check if all records in this dataset indeed are mined from GenBank. If the above line was written correctly, there should now be just one unique value in the variable institution_storing. (Note that for some other records, the institution storing is where the actual specimen sequenced is housed, i.e. a biological repository.) It is wise practice always to check; it is very easy to make a syntax error! We don't want to assume our code is correct and move on without checking.
unique(dfBOLD.GenBank$institution_storing)

#Example solution #2: We can practice doing the same thing as above in tidyverse style instead if we wish, assigning the output to a dataframe (tibble style) with a different name, so that we can compare. Notice I am now using the tidyverse function str_detect(), from the stringr package, rather than grep() from base R.
dfBOLD.GenBank2 <- dfBOLD %>%
  filter(str_detect(institution_storing, "GenBank"))

#Are these solutions the same? i.e. Do our resulting objects contain the same processid's (i.e. the same samples)? Yes.
all.equal(dfBOLD.GenBank$processid, dfBOLD.GenBank2$processid)

#Full objects are not quite the same. Attributes not quite exactly the same. so, we need to be careful what information we are actually interested in comparing when using all.equal(). Above, I was interested to see if the same processid's were recovered.
all.equal(dfBOLD.GenBank, dfBOLD.GenBank2)
class(dfBOLD.GenBank)
class(dfBOLD.GenBank2)
attributes(dfBOLD.GenBank)
attributes(dfBOLD.GenBank2)

#Now, to calculate a proportion, we can find the number of records with country missing and divide this by the total number of records. Answer: 0.6266, i.e. 62.7%. (Remember not to report results in more significant digits than what makes sense given the amount of input data.)
sum(is.na(dfBOLD.GenBank$country)) / length(dfBOLD.GenBank$processid)

#Look at the above again, and ensure you understand why we need sum() in the numerator and length() in the denominator.

#Now ... repeating our analysis for records OTHER than those from GenBank. By using a minus sign "-", we are omitting those index positions having GenBank in the institution_storing column. This could be solved in one line, but this would be harder to read, in my view. So, if we create an object that we don't need for downstream analysis, we can remove it right after, using rm()

#This line is finding records NOT from GenBank.
dfBOLD.not.GenBank <- dfBOLD[-c(grep("GenBank", dfBOLD$institution_storing)), ]

#This line is calculating the proportion of records that are missing country data from among those NOT from GenBank.
sum(is.na(dfBOLD.not.GenBank$country)) / length(dfBOLD.not.GenBank$processid)

#Answer: 0.0011, or 0.11%. So, records mined by BOLD from GenBank often are missing the country information. That doesn't mean the information is necessarily unknown; it could be in the source paper associated with the record or somewhere in the GenBank record, which could be processed. We will later have lessons on accessing NCBI data. By contrast, researchers participating in the Barcode of Life project most commonly upload such data to BOLD. Geographic data in a standardized format available for download is very helpful for large-scale geographic analyses!

#We could do the above task in other ways, including using piping and tidyverse style R programming.

#CHALLENGE #8 - Create a filtered dataset that only contains rows that have country data and that have latitude and longitude data and that have a BIN.

#Example solution using tidyverse style: 
dfBOLD.subset <- dfBOLD %>%
  filter(!is.na(country)) %>%
  filter(!is.na(lat)) %>%
  filter(!is.na(lon)) %>%
  filter(str_detect(bin_uri, "BOLD:"))

#Again, there are a multiple ways to solve this. For applying multiple conditions, I would tend to use piping to make a logical flow that is easy to read, even if this does end up with code that is longer in terms of number of line. Piped code tends to be longer but narrower. Often, it can be easier to write and to read.

#checking our data subset
sum(is.na(dfBOLD.subset$country))
sum(is.na(dfBOLD.subset$lat))
sum(is.na(dfBOLD.subset$lon))
sum(is.na(dfBOLD.subset$bin_uri))

#CHALLENGE #9 - Summarize the count of records by country in the data subset created in step #8.

#Example solution:
table(dfBOLD.subset$country)

#We can readily sort our table using the following line, which outputs the results to the screen. Note that I changed a default argument of sort() (the default is decreasing = FALSE), as I thought it would be more helpful to have the counts in descending order to address our question.
sort(table(dfBOLD.subset$country), decreasing = TRUE)

#We can add onto the above line to ask for the first five elements only if we wish, i.e. the countries with the highest counts. Answer: Canada, Mexico, Russia, United States, Antarctica.
sort(table(dfBOLD.subset$country), decreasing = TRUE)[1:5]

#CHALLENGE #10 - Using this subset, practise spreading the data by country using pivot_wider() and use specaccum() to create a species accumulation curve, using BINs as a proxy for species identification. Following similar code as above example in part 5.

#Example solution:

#Getting counts by country and BIN. See object to confirm.
dfBOLD.subset.by.country <- dfBOLD.subset %>%
  group_by(country, bin_uri) %>%
  count(bin_uri)

#Spreading so that countries are in rows, BINs are in columns, and counts are in cells. Again, check the object.
dfBOLD.subset.spread.by.country <- pivot_wider(data = dfBOLD.subset.by.country, names_from = bin_uri, values_from = n)

#Changing NAs to zeroes in the spread data set.
dfBOLD.subset.spread.by.country[is.na(dfBOLD.subset.spread.by.country)] <- 0

#Removing country names from a column and setting country names as row attributes.
dfBOLD.subset.spread.by.country <- dfBOLD.subset.spread.by.country %>%
  remove_rownames %>%
  column_to_rownames(var="country")

#Building accumulation curve.
dfBOLD.subset.accum <- specaccum(dfBOLD.subset.spread.by.country)

#We can make a plot of the accumulation curve.
plot(dfBOLD.subset.accum)

#We can fix the axis labels to be more accurate:
plot(dfBOLD.subset.accum, xlab = "No. of Countries Sampled", ylab = "BIN Richness")

#Further investigation: Look at the document. What are the vertical bars? Answer: confidence interval.

#Concluding remarks: I hope that you can imagine that such data subsetting and exploration steps would be very commonly-used tasks! I hope that you are now more comfortable with subsetting and filtering data, using both base R style and tidyverse style. For more information, I suggest to read "R for Data Science" chapters as we go.

