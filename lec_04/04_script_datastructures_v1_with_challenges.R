#######SCRIPT 2 FOR SOFTWARE TOOLS----

#Introduction to RStudio and R Data Structures (Part 2)

#Purpose: This script is designed to help us walk through the basics of using RStudio. RStudio is a friendly yet powerful software for conducting analyses in R and also has tools to facilitate software development. We will also go over 4 fundamental data object classes in R. Understanding data structures is essential for being able to use R. Many functions only operate on selected types of data, and understanding the structure of data will help us to avoid errors.

#Tip: If you've not done so already, please be sure first to review the uploaded script and example answers to the coding challenges for script #1, posted under class 1 on CourseLink.

#Tip: Please let me know in class if you have trouble seeing due to font size, and I will adjust.

#Tip: For beginning programmers, please don't worry if you find coding challenge #4 too difficult to solve immediately. We will get to many of those functions and options over time. Many of our example scripts will include a more advanced challenge as our class consists of learners who are coming in with various levels of experience. 

#Today, we will start using RStudio, which is an IDE (Integrated Development Environment) for R. I've posted the RStudio "cheat sheet" to CourseLink so that you may explore other features. However, we can productively get started by using a few key features of the software.

#For your reference again, here is the link to the formal language definition:
#https://cran.r-project.org/doc/manuals/R-lang.html
#Definitions provided in this script are from this resource.

###Part 1 - Introduction to Working in RStudio----

#Set RStudio to wrap commenting. Tools -> Global Options -> Code -> Soft-wrap R source files. This will mean you will be able to read the comments, without them running off the right side of the script display window.

#Also, notice that you can click on the arrows on the left that are beside the section headers. You can collapse and expand each section of the script, making it easier to navigate and to focus on one section at a time. You can use this formatting in your own scripts too. You need to add four or more dashes at the end of line to designate it as a header.

#While getting used to RStudio, we are going to work with some simple functions and commands and build upon class 1.
#4 panes or windows: script editor, console, environment, display window (e.g. for plots and documentation)
#In the script editor, use # to indicate comments
#We type our script in the script editor and send it to the console using Ctrl-Enter. I suggest to get used to using common keyboard shortcuts. You can put your cursor anywhere on a line and send it to the console. You can also highlight multiple lines and run or even select all (using control-A) and run the entire script. Lines of commenting will be skipped over.

#Let's try running a few lines of code similar to the class 1 script. What are we doing here? Predict what this line will do before you run it.
y <- c(1:10)
y
z <- 1:10
z

#YOUR TURN: Create a vector of integers from 1 to 20 and assign your data to a vector named my.vector

my.vector <- 1:20
my.vector

#See V2 for example answers.

#What do you notice about the colour of the comments vs. the code? I find the colour difference helpful for working with scripts.

#RStudio helps with finding errors. Try leaving out the closing parenthesis and see what happens.
#Uncomment (i.e. remove the #) and try the below
y <- c(1:10) #added comma here because it was missing! Rest of the code shows that I am missing a bracket if I do not resolve this issue.

#See the red symbol that appears to the left? pretty helpful. What was missing? Yes, the closing parenthesis.
y <- c(1:10)
y

#Also, check out how RStudio helps you keep track of parenthesis pairs by scrolling through this line and also writing your own line that is similar.
z <- c((3:10), (y + 2))

#Predict what the above line will do before you uncomment and run the below (i.e. remove the # and run the below line to check).
z

#Look in your global environment (typically, the pane to the upper right, although you can re-arrange your view). What do you see? RStudio will help you to keep track of your objects. Note these are in working memory. In other lessons, we will also see how to write data to hard disk.

#To remove one object:
rm(y)
rm(z)

#Or, remove two:
rm(y, z)

#You can also use the sweep symbol in your environment pane if you want to clear out everything with your mouse. For a larger script that you want to be reproducible, rather than interactive, code in the removal of your objects. You can remove individual objects or a subset or all. To remove all:
rm(list = ls())

#Now, let's create some data. It is very helpful to be able to create data quickly to work with. Additionally, at some point, you may wish to generate data from various statistical distributions for real investigations (e.g. you might wish to compare real data to a null distribution). R has many different distributions available (e.g. random, uniform, poisson, etc.), but we will use just one for today.
x <- rnorm(1000)
y <- rnorm(1000)

#See x
x

#See the top of x
head(x)

#See the end portion of x
tail(x)

#See how many elements are in x
length(x)

#What does rnorm() do? As always, check out the documentation. What are the defaults for this function?
?rnorm

#YOUR TURN: generate 100 numbers under a random normal distribution, with a mean of 10 and standard deviation of 5. Assign your result to a vector called my.rnorm
my.rnorm <- rnorm(100, mean = 10, sd = 5)


#See V2 for example answer.

#Let's learn how to use some more simple, yet commonly used, functions for mathematics. Check out the mean, median, and standard deviation of your results.
mean(my.rnorm)
median(my.rnorm)
sd(my.rnorm)
max(my.rnorm)
min(my.rnorm)

#What do you think these functions do? See! There are many useful functions in R with simple, intuitive names that are easy to remember.

#Do the answers you are getting look reasonable? On your own, have a look at the documentation for these functions.

#R can generate data under various distributions. This can be very helpful for generating null models to compare against real, observed data or for generating random numbers to help you with experimental design. Some other examples of distributions you can sample from: random uniform using the function runif(), random numbers under a poisson distribution using rpois(), binomial distribution using rbinom(), exponential distribution using rexp(), etc. Remember that R has its origins in statistical computing and has a lot of functionality in stats, although it's now used for many types of analyses and visualizations.

#Now, let's create and view some simple figures in the lower-right window. Later we will see ways to customize figures and create more visually appealing figures using the ggplot package. hist() is used to build a frequency distribution of values, a histogram, while plot is used when you have two variables to plot.
hist(x)
hist(y)
plot(x, y)

#On your own: How do our plots look if we generate just 10 values?
z <- rnorm(10)
a <- rnorm(10)

plot(z,a)

#let's build a linear model, summarize the model, and plot the regression line.

#We are creating a linear model and assigning it to the name model. We are using the function lm() to fit the model. The response variable (or "dependent variable") comes first, here y. This is followed by the tilde, i.e. ~. The predictor variable (or "independent variable") come after the tilde. (There can be multiple predictor variables in the case of multiple regression, but for now we will use just one predictor variable.)
model <- lm(y ~ x)

#Do we expect a significant trend? Why or why not?

#ANSWER: No, we do not expect a significant trend. We generated x and y independently under the random normal distribution, with no association between them.

#Let's see. Summarizing the model:
summary(model)

#adding the regression line for the model to the plot.
abline(model)

#Is there a significant trend? What do we find in the model summary and in the plot?

#Now let's view the model diagnostic plots
plot(model)

#Note you then need to hit enter multiple times within the console window in order to scroll through the plots.
#What are we looking for? We want to see if our model assumptions are plausible. Is a linear model suitable for these data?
#Please review from your prior stats course(s) the meaning of these plots and regression assumptions. I suggest that the following online tutorial is helpful for understanding diagnostic plots for linear models:
#https://data.library.virginia.edu/diagnostic-plots/

#Help will also show up in the lower right window. Try looking at the help and check out the default arguments.
?plot()

###Part 2 - Matrices and Coding Challenge!----

#Remember that you can work in small groups if you wish to write code to perform the below tasks. Or, you could complete this individually and then compare with your peers and with the posted example answers. I will typically post example answers within a few days, but please do remind me if I ever forget!

#Generate the integers from 1 to 100 and assign to a vector named hundred.
hundred <- 1:100

#Show the contents of hundred on the screen
hundred

#Check out the documentation for the function called matrix()
?matrix()

#Take the contents of hundred and place into a 10 X 10 matrix called ten.by.ten
ten.by.ten <- matrix(data = hundred, nrow = 10, ncol = 10)

#Check out the contents of ten.by.ten on the screen
ten.by.ten

#Indexing is used very frequently for analyzing subsets of data. Here are examples of indexing by position for matrices... 

#We can access rows from matrices as follows:
ten.by.ten[2, ] #row 2
ten.by.ten[4, ] #row 4
ten.by.ten[1:4, ] #rows 1 through 4
ten.by.ten[c(1, 3, 5), ] #rows 1, 3, and 5

#Above, by putting nothing after the comma, this means we want all columns.

#We can pull out columns from matrices as follows:
ten.by.ten[, 2] #column 2
ten.by.ten[, 4] #column 4
ten.by.ten[, 1:4] #columns 1 through 4
ten.by.ten[, c(1, 3, 5)] #columns 1, 3, and 5

#By putting nothing before the comma, this means that we want all rows.

#We can also select one individual elements. This can be considered similar to selecting a single cell in an Excel sheet.
ten.by.ten[2, 3] #row 2, column 3
ten.by.ten[3, 5] #row 3, column 5

#Examples of selecting more than one element:
ten.by.ten[1:5, 1] #rows 1 through 5, column 1
ten.by.ten[1, 1:5] #row 1, columns 1 through 5

#Rewrite your line for creating your matrix two ways: with vs. without naming the arguments. (What do you remember about specifying the arguments and argument order? Tip: You don't have to name the argument if you will use them in the default order. However, I recommend to name your arguments when you are using a function that takes lots of optional arguments, so that everything goes as planned. Also, naming your arguments makes your code explicit and increases its future utility.)
matrix.test <- matrix(hundred, 10, 10)

#Are they the same?
all.equal(ten.by.ten, matrix.test)

#Now, depending upon what we are doing, we may want to fill by row instead of column. The default for byrow = FALSE. Try again, setting byrow = TRUE. See how those two options differ.
matrix.test2 <- matrix(data = hundred, nrow = 10, ncol = 10, byrow = TRUE)
matrix.test2

#Create your own example of using the sample() function.
?sample()
my_sample <- sample(200, 100)
my_sample
#Code it out!

#Also, here are 3 specific examples to try, plus a more advanced challenge:

#Example #1: In your hypothetical study of environmental DNA (eDNA), there are 200 potential ponds you could include. You can have numbered the ponds 1-200. The scope of your study is to sample 40 of them. Pick 40 of 200 ponds randomly using R.
ponds <- sample(x = 1:200, size = 40)
ponds
length(ponds)

#Tip: When sampling or generating random numbers, if you want your script to be completely reproducible, you can use set.seed() first. First, try generating several sets of samples without doing that. Then, try set.seed() and do it again.

ponds_2 <- sample(x = 300, size = 23)
ponds_2

set.seed(45)
ponds_2 <- sample(x = 1:200, size = 40)
ponds_2


#See V2 for example answers.

#Example #2: In this hypothetical example, you have 10 friends but only 5 extra cookies. Create a vector containing the names of your friends. Tip: You should be creating an atomic vector of data mode character. Randomly pick 5 to whom you will deliver a cookie. (The remainder would presumably receive a cookie next week!)
vc_friends <- c("sabrina", "farah", "sharon", "eva", "rebekah", "jimmy", "brandon", "remi", "liona", "laura")
who_gets_a_cookie <- sample(vc_friends, 5)
who_gets_a_cookie

#See V2 for example answer.

#Example #3 (hardest). Imagine you have collected Trichoptera (caddisflies) larvae from your ponds. Caddisflies are important organisms for biomonitoring of freshwater quality and ecosytem health, due to their sensitivity to pollution. From each pond, you collected 100 caddisflies, which you will number 1 through 100. From each of your 40 ponds, you will randomly choose 30 of the 100 caddisflies for DNA sequencing. You can build up to this by doing the different components in different steps, but then also try to do this in ONE line of code. Tip: check out the function replicate(). Tip: Note that you can nest functions inside other functions. See examples at the bottom of documentation.

#My multi-step answer:
ponds <- (1:40)

for (i in ponds) {
  fly_seq <- sample(x = 100, size = 30, replace = FALSE)
  print(fly_seq)
}

#Answer with replicate(): 
options(max.print = 4000)
replicate(40, sample(x = 100, size = 30, replace = FALSE))

#See V2 for example answer

# creating Objects as arguments for the list function
a = c(1, 2, 3, 4, 5)
b = c("A", "B", "C", "D", "E")
c = c("A", "B", "c", "1", "2")

# implementing the list() function 
list(a, b, c)

?list()

#Check out class and dimensions of that object.
class(a)
dim(a)
max(a)
min(a)
mean(a)

#TIP: See V2 for tips on using viewer.

#Challenge #4 and friendly competition! Come up with as many (reasonable) ways as possible to solve example 3. In R, there are many ways to do the same thing! Sometimes, different solutions are a matter of personal preference. Sometimes, there are better and worse ways to do things (e.g. considering ease of reading the code and computational time).

#1: for loop
ponds <- (1:40)

for (i in ponds) {
  fly_seq <- sample(x = 100, size = 30, replace = FALSE)
  print(fly_seq)
}

#2: replicate()
options(max.print = 4000)
replicate(40, sample(x = 100, size = 30, replace = FALSE))

#3: while loop
count <- 0
while (count >= 0 & count < 40) {
  count <- count + 1
  fly_seq <- sample(x = 100, size = 30, replace = FALSE)
  print(count)
  print(fly_seq)
}

#I will put example answers to challenge #4 in a separate script, so that this script doesn't get too long. Everyone: Feel free to have a look once you are comfortable with this script.

#Code out these examples, and then create your own example that is relevant for your own life or for research you wish to pursue!

#As a review, here is a very brief post about the main data structures in R. It is often helpful to review similar information in different ways to help us to solidify our understanding of the most important concepts.
#https://jamesmccaffrey.wordpress.com/2016/05/02/r-language-vectors-vs-arrays-vs-lists-vs-matrices-vs-data-frames/


###Part 3 - Common Object Classes: Atomic Vector, Matrix, and List----

#We will review a few concepts from last time and expand our toolkit of object classes. By the end of the lesson you should be able to explain the difference between a vector and a list and a matrix and a dataframe.

#Revisiting object class definitions. Vectors are "contiguous cells containing data". An atomic vector contains one type of element. You can play around.

vect.numb <- c(1, 8, 12, 15, 18)
vect.numb
class(vect.numb)

vect.char <- c("Analysis", "in", "R", "is", "fun")
vect.char
class(vect.char)

vect.mixed <- c(1, 3, 5, "analysis", "fun")
vect.mixed
class(vect.mixed)

#What class did R determine as the class for this vector? Character. Note that numerical data can be treated as character. See data types and their hierarchy in the base R cheat sheet.

#Indexing allows us to pull out selected data elements. Note that we use square brackets for this. For example, check out:
vect.char[1]
vect.char[2]
vect.mixed[3]
vect.mixed[c(1, 3)]

#Now, enter numeric data but use the function list(). So, note that "a list" is a data object of the class list, but the function list() is a function that creates a list! A list is a generic vector. It can still be considered a linear data object. It may (but doesn't have to) contain more than one type of element.

listA <- list(1, 45, 3, 5, 9)
listA #How do the data look on the screen compared to a vector?
#it shoes you the order of the values / index of the list
class(listA) #What do you notice? what class of data object is returned?
#list is it's own data object, whereas vectors hold elements of the same object class.
# Q: DOES THIS MEAN THAT VECTORS CAN HOLD LISTS? <- NO BECAUSE THEY CAN ONLY HOLD ATOMIC ELEMENTS


#What does the function list do?
?list()
#"Functions to construct, coerce and check for both kinds of R lists." (from help)

#What is a list? "generic vector". It is a more general class of data object. The elements in a list could each be a matrix, for example. Thus, a list can contain hierarchical data and can also contain different types of data. Because lists can house hierarchical data, indexing can be a little more complex for lists than atomic vectors or matrices, for example, but we will get there in due course.

#It is always important to check what class of object you have. Many functions only act upon a specific class of object. Objects can be converted from one form into another, if compatible. For example:
mylist <- as.list(vect.char)
mylist
class(mylist)

#Lists are an extremely useful data object class in R! We will be returning to lists in future lessons and learn how to conduct indexing using lists and hierarchical data structures.

#Now moving on to rectangular data objects. This is a very familiar type (think an Excel sheet, which contains rows and columns). Here, we are repeating code similar to that from part 2 above.

#Using the function matrix(), we will create a new object, of the class matrix. A matrix is a two-dimensional rectangular data object which should contain one type of element. We will do this by sampling from the vector called students. Our matrix will have 5 rows. Note that I didn't specify the "size" argument as I want to include ALL of the students in a group. By reading across the rows, we will see which hypothetical students are assigned to which groups.
students <- c("Xin", "Jane", "Emma", "Maria", "Joe", "Kim", "Bev", "Li", "Sal", "Francois", "Marie", "Jen", "Vic", "Ang", "Jose", "Len", "Viv", "Don", "Kev", "Rick")
length(students)
groups1 <- matrix(sample(x = students), nrow = 5)
groups1
class(groups1)

#Think: In this case, would we want to sample with or without replacement? See V2 for answer.
#without, because you don't want to add the same people to multiple groups

#Reminder about viewing the help.
?matrix

#We can add arguments if we want to when using the matrix() function. In this case we will specify the argument dimnames to add names for the rows and columns.
groups2 <- matrix(sample(x = students), nrow = 5, dimnames = list(c("group.Awesome", "group.Brave", "group.Clever", "group.Diligent", "group.Enthusiastic"), c("Col1", "Col2", "Col3", "Col4")))
groups2

#Why CAN'T we do below? i.e. Why can't we use c() rather than list()? What does c() do? Would it preserve the dimensions? Try this. Note we are expecting an error so that you can see what happens. So, uncomment and try to run the below.

groups2 <- matrix(sample(students), nrow = 5, dimnames = c(c("group.Awesome", "group.Brave", "group.Clever", "group.Diligent", "group.Enthusiastic"), c("Col1", "Col2", "Col3", "Col4")))

#Try to explain why we can't do this. WHY do we get an error? Would we want all of our column and row names in one vector?
#c() creates a vector, whereas list() creates a list (that can hold multiple dimensions, like both the names of the rows and the columns). We want both the column and row names in one vector so that we can translate this list back into a matrix if needed, and reproduce our work.

#To follow the above code, read from the inside out. Start with the most nested set of parentheses. The function c(), or concatenate, is used to create two vectors, one for row names and one for column names. The function c() is a "generic function which combines its arguments" (from the help). "The default method combines its arguments to form a vector" (from the help). These vectors are then combined into a list object using the function list(). The list contains two elements, each of which is a vector (one for row names and one for column names). So, we need to use list(). 

#Tip: Try placing the parentheses in different places and see what happens. Parenthesis placement is very important and always needs to be in pairs! Try playing with vectors vs. lists.

#Now, let's look at our new and improved matrix.
groups2

#Much better! We will next use data generated by the class to demonstrate dataframes.

###Part 4 - Dataframe and Our Class's Background in R----

#In class - check time.

#The last type of data object we will look at today is a dataframe. A dataframe is like a matrix in that it is a rectangular data object. However, it can contain more than one type of element. It is one of the most common types of data object classes that we will work with.

#Typically, samples are in rows, and variables are in columns. Each column should contain just one type of data (numerical, character, categorical, etc.) to facilitate analysis. That is called "tidy data", and we will discuss that idea further in a future class.

#First, we will create a vector using familiar code, using data generated by and about the class. (These numbers now updated for 2021 class from Quiz Zero.)
background <- c(10, 17, 3, 1)
?names()
names(background) <- c("Beginner", "Novice", "Intermediate", "Advanced")
background
class(background)

?names

#For fun, let's make a pie chart.
pie(background)

#You can change the colours if you want, in various ways.
pie(background, col = rainbow(length(background)))

#Notice that the above is more general than the following, which also works.
pie(background, col = rainbow(4))

#Tip: In general, it is wise to make our code more general. We want to build code that can be easily rerun if the data change. Often, we will run the analysis many times as we develop a project and also as we add data.

?pie() #See the note about pie charts being less reliable to interpret visually than dot plots or bar charts.

#Now, let's create a data frame instead. This is a rectangular data object. Here, "level" will be treated as data (i.e. have its own column). For analysis purposes, that is generally preferable to using the names attribute (like a label "outside" the data).
df <- data.frame(c("Beginner", "Novice", "Intermediate", "Advanced"), c(9, 17, 3, 1), stringsAsFactors = TRUE)
df
names(df) <- c("Level", "No.Students")
df
class(df)

#Click on "df" in your environment window. See the dataframe using the tab next to your script in the script editor. This is very handy to be able to view the data in this way.

plot(df, xlab = "Level", ylab = "Number of students")

?plot()

#Note: We will learn how to make much nicer-looking figures with ggplot2 soon.

#Have a look at the $ operator. Separates the name of the data frame from the name of the variable. We can look at specific variables in a dataframe this way.
class(df)
class(df$Level)
class(df$No.Students)

#What is the class of each of the variables in our dataframe df?
#Level is a factor
#No.Students is numeric

#Parting note: ls() is very different from list()
ls()
#What does this do?
#ls() lists all the values in the environment
?ls()



#This was a brief overview... We will work with dataframes much more and become comfortable with them over the next few classes.


#Before next class, please install needed packages.
#If you don't have the below packages, you would uncomment (i.e. remove the # symbol) and run the following lines
install.packages("tidyverse")
install.packages("vegan")
install.packages("iNEXT")
