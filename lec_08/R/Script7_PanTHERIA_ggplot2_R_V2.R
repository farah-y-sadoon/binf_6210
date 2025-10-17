###### **********************
## CLASS 7 V1 - SOFTWARE TOOLS 2021 - INTRODUCTION TO GGPLOT2----
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

#V2: last updated September 25, 2021.  Version 2 contains some tips and example answers for part 6.

####1 - OVERVIEW----

#The purpose of this script is to provide a brief introduction to plotting in the ggplot2 package, which is among the most popular packages in all of R. It is valuable for enabling users to produce appealing, informative graphics in a flexible framework.

#This script also provides a brief introduction to using regular expressions for cleaning up character data. (We will be returning to that subject later.)

#We will also work with a new biological data set today. The PanTHERIA dataset was a major collaborative effort by Jones et al. (2009) to compile species-level biological information for mammal species.
#https://ecologicaldata.org/wiki/pantheria
#http://esapubs.org/archive/ecol/E090/184/#data

#I chose the following file for us to use: "PanTHERIA_1-0_WR05_Aug2008.txt -- 5416 records, not including the header row, ASCII text, tab-delimited, no compression scheme was used."

#Here was the code I used originally to acquire the data (data download performed on Oct 28, 2019):
#df <- read_tsv("http://esapubs.org/archive/ecol/E090/184/PanTHERIA_1-0_WR05_Aug2008.txt")

#Code I used to write the data to hard disk so that we have a copy of the original data saved. I have also uploaded this file to CourseLink.
#write_tsv(df, "PanTHERIA.tsv")

####2 - READING DATA----

#Remember to set your working directory
#getwd()
#setwd()

#Reading in data
dfTraits <- read_tsv(file = "../data/PanTHERIA.tsv")

#If this was successful, you will see an object in your environment with 5416 observations of 55 variables.

#Have a look at this object in the viewer.

####3 - EXPLORING DATA AND CHECKING FORMATTING----

#check class
class(dfTraits)

#see variables
names(dfTraits)

#see summary
summary(dfTraits)

#Let's look at an example histogram. Let's use element 8 (adult forearm length). Remember we can use double square brackets followed by the index position to extract that single element (column), and thus pass that atomic vector to a function such as hist() that takes numerical data.
hist(dfTraits[[8]])

#compare - uncomment and try this. Remember, here for a tibble (tidyverse data frame) when we use single square brackets and and numerical indexing we get a smaller dataframe
#hist(dfTraits[, 8])
#class(dfTraits[, 8])

#Let's look at the data in the viewer and also look at the summary. What do you notice? Lots of -999. That's commonly used to indicate missing data! The reason is because there is almost no real measurement that would give that value; thus, -999 is a placeholder showing missing data. However, for data analysis, that doesn't work for us, as we don't want -999 mixed in with our real numeric data. Let's clean this up.

#Making a copy so we keep original df
df1 <- dfTraits

#Let' replace -999 values with NA
df1[df1 == -999] <- NA

#Let's compare original and edited data side by side, by looking at column 6. We will bind column 6 from the original tibble with column 6 from the new tibble. Then, we will look at the first 100 rows.
check <- cbind(dfTraits[, 6], df1[, 6])
check[1:100, ]

#If a comma isn't used, the default when using numeric indexing of data frames is for columns to be used. So, this will be the same. Check for yourself. However, I generally recommend to be more formal and use two-dimensional indexing when using single square brackets for data frames. I find this can help code to be more explicit and reduce errors.
#check2 <- cbind(dfTraits[6], df1[6])
#check2[1:100, ]
#rm(check2)

#You may notice that I often spend about half of my code checking that the other half is working as intended! I do recommend always to check that the code did what you mean it to do.

#Let's look at a summary. We can see that there are a lot of NAs. So, thoughtful decisions about the treatment of NAs would be needed for statistical analysis. Do you need complete cases for selected variables? Will you remove NAs only for single variables of interest? (In a future lesson, we will also discuss imputation strategies.)
summary(df1)

#Let's have a look at some of the data. We would want to investigate the distribution of the data prior to any statistical analysis.

#Typically, in R, variables are not supposed to start with a number. They should start with a letter or a . (which I would recommend to save only for special cases). In base R, we can use single or double quotation marks to work with variable names starting with numbers, e.g. df1$`5-1_AdultBodyMass_g`. However, I find this clunky to work with, to always have to use ticks for these column headers. Instead, we can readily clean up the variable names before working with the data. We want to be able to use the variable names easily for working in ggplot2.

#Example of using base R functions with variable names starting with a number.
hist(df1$`5-1_AdultBodyMass_g`)
mean(df1$'5-1_AdultBodyMass_g', na.rm = TRUE)

#However, we'd like to be able to access ggplot2 graphics easily. So, our specific task is: clean up the column headers for dataset PanTHERIA. What is the GENERAL kind of problem we are trying to solve? R doesn't care that these are biological variables. What we have here is a STRING problem. We need tools for cleaning up the strings: regular expressions.

#Let's look at the names. We can see these are structured names.
names(df1)

#Let's work with the names separately for now, in their own character vector.
names.original <- as.vector(names(df1))
names.original
class(names.original)

#check documentation
?str_replace()

#First, we can readily clean up the first 5 column names by replacing the specific starting string with nothing. The hat, "^", means we are looking at the beginning of the string (in case that string shows up later in one of the column headers, but that actually isn't needed here... I include the ^ for an example).
names.edited <- str_replace(string = names.original, pattern = "^MSW05_", replacement = "")

head(names.edited)

names.edited

#Now, what pattern do you see in the names of the remaining variables? They all start with a 1 or 2-digit number. Then, there is a hyphen. Then, there is another digit, and then there is an underscore. So, these names have a very specific structure, making them easier to work with.

#There are a couple of ways we could fix these. Here, we introduce the character "^" which means to start at the beginning of the string. This isn't strictly needed here, but this usage can be helpful to ensure we didn't overlook something at the end of the strings that would conflict with the edits we are trying to make. We know we are trying to get rid of the numbers and punctuation at the beginning so that our column names start with letters. In a regular expression, "+" means one or more. So, we are looking for one or more digits. Also, in regex, [] means OR. So, we are looking for 0 or 1 or 2, etc. 0-9 means all digits.
names.edited <- str_replace(string = names.edited, pattern = "^[0-9]+-[0-9]+_", replacement = "")

names.edited

#Now, we can see that we have cleaned up the beginnings of the variable names, and the column names now start with a character and not number.

#assignment of new names of columns to df1
names(df1) <- names.edited
names(df1)

#Looking good! OK, now these names are easier to work with.

#Please note that, in a research project, we would do more extensive data exploration before heading into statistical testing. Examples would include examining distributions of all variables. Outliers (through both visualizations and statistical metrics) would be explored. This does NOT mean we would automatically delete outliers. We would not want to throw away real biological data, but we might check for data entry errors. As well, we can consider whether data transformations (e.g. log, ln) may be relevant for our data set. But let's move on for now into our plotting lesson. We will conduct plotting of just a small number of the variables as part of our introduction to ggplot2.

#ending this section with removing objects not needed downstream.
rm(check, names.original, names.edited)

####4 - INTRODUCTION TO GGPLOT2----

#Let's make some plots! Please note that "+" is used at the end of the line for ggplot2, rather than the pipe we have seen for other tidyverse packages "%>%". If you are interested in this, you can read more about that directly from Hadley Wickham, with that post part way down the following link: https://community.rstudio.com/t/why-cant-ggplot2-use/4372/5

#Let's make a plot of type histogram for body size
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = AdultBodyMass_g)) +
  labs(title = "Histogram of Mammal Body Mass", x = "Adult Body Mass (g)", y = "Count")

#Look at very large x-axis! We can also try looking at the data in boxplot format, and separating out the orders to see where the very large values are coming from. What do you predict?
ggplot(data = df1) +
  geom_boxplot(mapping = aes(x = Order, y = AdultBodyMass_g), outlier.color = "red", outlier.fill = "red", outlier.size = 2) +
  labs(title = "Boxplot of Mammal Body Mass", x = "Mammal Orders", y = "Adult Body Mass (g)") +
  coord_flip()

#We can see that very large sizes are present in one order (Cetacea), which makes it hard to see the distributions of body sizes in other orders.

#Let's create a new dataframe omitting Cetacea
df1.xCetacea <- df1 %>%
  filter(!Order == "Cetacea")

#repeat plot without the Cetacea and this time with the orders in alphabetical order.
ggplot(data = df1.xCetacea) +
  geom_boxplot(mapping = aes(x = Order, y = AdultBodyMass_g), outlier.color = "red", outlier.fill = "red", outlier.size = 2) +
  labs(title = "Boxplot of Mammal Body Mass", x = "Mammal Orders", y = "Adult Body Mass (g)") +
  coord_flip() +
  scale_x_discrete(limits = c(sort(x = unique(df1.xCetacea$Order), decreasing = T)))

#So, we have large variability in body size. To look at the relationship between body size and other parameters, we may wish to perform a data transformation.

#Note! log() used alone means ln (natural logarithm) in R. If you want log base 10, then we need to use log10(). Compare:
log(10)
log10(10)

#Let's have a look at the distribution when log10-transformed. Still right-skewed, and reflecting there are some very large values in the dataset, but this is much more reasonable to work with for plotting and statistical analysis.
hist(log10(df1$AdultBodyMass_g))

#If we want, we could create a new variable onto the end of our dataframe. Or, we can specify log in our ggplot function. We could do this if we wanted to work with the log values a lot more downstream:
#df1 <- df1 %>%
#  mutate(logBodyMass = log10(AdultBodyMass_g))

#We can build a ggplot
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = log10(AdultBodyMass_g))) +
  labs(title = "Histogram of Mammal Body Mass", x = "Log10 of Adult Body Mass (g)", y = "Count")

#OK. Now looking at maximum longevity
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = MaxLongevity_m)) +
  labs(title = "Histogram of Longevity in Mammals", x = "Maximum Longevity (m)", y = "Count")

#Looking at ln-transformed longevity. A histogram is for examining the the distribution of values for ONE quantitative variable.
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = log(MaxLongevity_m))) +
  labs(title = "Histogram of Longevity in Mammals", x = "(ln) Maximum Longevity (m)", y = "Count")

#Now, plotting body mass vs. longevity as a scatterplot. This is a suitable type of plot when you are looking at the relationship between two quantitative variables.
ggplot(data = df1) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m))) +
  labs(title = "Body Size vs. Longevity in Mammals", x = "(log10) Adult Body Size (g)", y = "(ln) Maximum Longevity (m)")

#Same plot, but now we are colouring by orders of mammals. Using colour and/or symbol can be an effective way to add a third variable to a two-dimensional plot. In general, I find such strategies more successful than trying to make "three-dimensional" plots, which often end up hard to read.
ggplot(data = df1) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Order)) +
  labs(title = "Body Size vs. Longevity in Mammals", x = "(log10) Adult Body Size (g)", y = "(ln) Maximum Longevity (m)")

#Very interestingly, we can see that bats (order Chiroptera) have high longevity in relation to their body mass. However, this is a little tricky to see due to the many groups in this plot.

#We could drop taxonomic groups with few species, focusing on the more diverse orders, to see relationships more clearly. Here, we are creating a filtered dataset, whereby only orders with 200 or more species are retained. Note that this is very easy to do using tidyverse functions.
df1.large.orders <- df1 %>%
  group_by(Order) %>%
  filter(n() >= 200)

#repeating our scatterplot for orders containing >= 200 species. Now we have fewer groups to look at.
ggplot(data = df1.large.orders) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Order)) +
  labs(title = "Body Size vs. Longevity in Species-Rich Mammal Orders", x = "(log10) Adult Body Size (g)", y = "(ln) Maximum Longevity (m)")

#Repeating our scatterplot for orders containing >= 200 species, this time using both colour and symbol shape to show the orders. We can use colour and symbols with redundancy to make plots more accessible.
ggplot(data = df1.large.orders) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Order, shape = Order)) +
  labs(title = "Body Size vs. Longevity in Species-Rich Mammal Orders", x = "Log10 Adult Body Size (g)", y = "Ln Maximum Longevity (m)")

#Very interestingly, we can see that bats have high longevity for their body mass compared to other mammals in this plot. An effect of order upon the body mass/longevity relationship could be tested statistically through multiple linear regression, i.e. with order and body mass as predictor variables and longevity as the response variable.

####5 - FIX THAT PLOT!----

#The purpose of this section is to help you to get familiar with a few of the key functions and arguments involved with making a few of the most common types of plots. I have created a few plots that have something WRONG with them. Figure out what is wrong and fix the code to correct the figure.

#Remember to Zoom in your view of the plots to see better.

#PLOT1. Body size vs. longevity for diverse mammalian orders (species with >= 200 species). I want to colour by a third variable, order. What is wrong with this figure? Can you figure out what has gone wrong and fix it?
ggplot(data = df1.large.orders) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Order, shape = Order), size = 20) +
  labs(title = "Body Size vs. Longevity in Species-Rich Mammal Orders", x = "Log10 Adult Body Size (g)", y = "Ln Maximum Longevity (m)")

#Answer: I set the "size" argument to be very large! this makes it very hard to read this figure. Play around with the size setting to get a plot that is easy to read.

#PLOT 2. This is the same general plot as above. We want to show body size vs. longevity, with order also conveyed through the use of colour and shape, but there is a different problem. What is it and how can I fix it?
ggplot(data = df1.large.orders) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Family, shape = Order), size = 2) +
  labs(title = "Body Size vs. Longevity in Species-Rich Mammal Orders", x = "Log10 Adult Body Size (g)", y = "Ln Maximum Longevity (m)")

#Answer: I set the colour argument to be Family instead of Order. This leads to a confusing plot with way too many colours. Lesson: Check match between code and what we are trying to do. Also, we need to ensure that figures are readable.

#PLOT3. Final version of this plot. What is the problem here? Look carefully at everything: what data are we plotting? How is the figure labeled?
ggplot(data = df1.large.orders) +
  geom_point(mapping = aes(x = log10(AdultBodyMass_g), y = log(MaxLongevity_m), colour = Order, shape = Order), size = 2) +
  labs(title = "Body Size vs. Longevity in Species-Rich Mammal Orders", x = "Log10 Adult Body Size (g)", y = "Litter Size")

#Answer: I have incorrectly labeled the y-axis as "Litter Size". Again, we have to make sure everything matches up. This is supposed to be a plot between body mass and longevity. Labeling is a critical component of figures! Very important!

#PLOT 4. Let's return to a histogram figure. This is intended to be a histogram of (ln) maximum longevity. There are TWO problems with this code. Can you find both?
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = log(LitterSize))) +
  labs(title = "Boxplot of Longevity in Mammals", x = "(ln) Maximum Longevity (m)", y = "Count")

#Answer: The plot title is wrong. The title says boxplot, but the type of figure is actually a histogram. Also, note the data argument to aes(). This is actually a plot of LitterSize. I referred to the wrong variable. Again, the data being plotted need to match the plot title and axis labels!

#PLOT 5. Let's actually build a histogram of litter size, but what is weird with this plot? It's rather hard to read. Can you fix it?
ggplot(data = df1) +
  geom_histogram(mapping = aes(x = LitterSize), binwidth = 0.1) +
  labs(title = "Histogram of Litter Size in Mammals", x = "Litter Size", y = "Count")

#Answer: I set the binwidth to be very small, yielding too fine-grained a view that is harder to interpret overall. We can increase this to a reasonable size. This is a very useful plot component to play with. Sometimes, the binwidths that are selected if you leave this argument as default might not be what we want. For example, it can be common to make binwidth a little smaller than default to show more details in the distribution of the data, but not too small.

#PLOT 6. Let's compare the distribution of litter size for the most diverse orders. If you want to compare distributions of data across groups, a great option is a frequency polygon. Fix the plot title.
ggplot(data = df1.large.orders) +
  geom_freqpoly(mapping = aes(x = LitterSize, colour = Order, fill = Order)) +
  labs(title = "Histogram of Litter Size in Mammals", x = "Litter Size", y = "Count")

#Answer: To be accurate, the title of the plot should be changed to something like: "Frequency Polygon of Litter Size in Mammal Orders".

#PLOT 7. Let's return to our boxplot of mammal body sizes, by order, excluding cetaceans. Remember to use "zoom" to see your plots. I have left out one function that I used above in section 5. What was it? And, what did it do? Put it back in and compare! Which is more readable?
ggplot(data = df1.xCetacea) +
  geom_boxplot(mapping = aes(x = Order, y = AdultBodyMass_g), outlier.color = "red", outlier.fill = "red", outlier.size = 2) +
  labs(title = "Boxplot of Mammal Body Mass", x = "Mammal Orders", y = "Adult Body Mass (g)") +
  scale_x_discrete(limits = c(sort(x = unique(df1.xCetacea$Order), decreasing = T)))

#Answer: I left out the flipping of coordinates using coord_flip(). Instead of that, another solution to this cluttered, unreadable x-axis would be to make the order names in vertical print. However, it is easier for the reader to be able to read left to right. So, in this case, flipping of the x and y axes could be a good option for readability.

#PLOT 8. Same plot. I've changed the colour and symbol size for displaying outliers. I've also REMOVED another function that I used above in section 5. What was it? What did it do? Which plot did you like better?

ggplot(data = df1.xCetacea) +
  geom_boxplot(mapping = aes(x = Order, y = AdultBodyMass_g), outlier.color = "turquoise", outlier.fill = "turquoise", outlier.size = 1) +
  labs(title = "Boxplot of Mammal Body Mass", x = "Mammal Orders", y = "Adult Body Mass (g)") +
  coord_flip()

#Answer: In this case, I removed the reordering of the order names. In the above example in part 5, the order names were placed in alphabetical order for readability. Data could be ordered in various ways, depending upon the message you want to convey. In some cases, ordering by value (e.g. smallest to largest) can make sense. In other cases, you might order by another variable to show relationships in your data. We have to think about what is the message we are trying to convey.

####6 - YOUR TURN----

#This section is your chance to express your curiosity and creativity!

#Have another look at what variables are available in this dataset:
names(df1)

#What variables are you interested in? Create plots to explore the relationships among variables and share!

#You could create plots with other variables. You can also explore creating new kinds of plots. Some links below may be helpful.

####7 - FURTHER RESOURCES, TUTORIALS, AND READINGS----

##R documentation, of course!

#Chapters 1, 3, 5, and 7 from "R for Data Science (2)" for an introduction to handling data, plotting in ggplot2, and data exploration.
#https://r4ds.hadley.nz/

#Overview of ggplot2 and cheat sheet
#https://ggplot2.tidyverse.org/

#Histograms and frequency density plots:
#https://ggplot2.tidyverse.org/reference/geom_histogram.html

#Barplots:
#https://ggplot2.tidyverse.org/reference/geom_bar.html
#http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization

#Scatterplots:
#https://ggplot2.tidyverse.org/reference/geom_point.html

#Boxplots
#http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization

#Violin plots
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

#setting themes for colour and font size
#https://ggplot2.tidyverse.org/reference/ggtheme.html

#Regular expressions: See below link for good information about regex. It is worth it to keep working on regular expressions, and incorporate usage of regex into your coding practices, when suitable. Over time, you will SAVE time overall if you spend some time on this now. Also, see the stringr cheat sheet, saved under "Cheat Sheets" on CourseLink.
#https://rstudio-pubs-static.s3.amazonaws.com/74603_76cd14d5983f47408fdf0b323550b846.html