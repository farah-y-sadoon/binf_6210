######Software Tools  - Getting Started with R----

###PART 1 - INTRODUCTION----

#Welcome! We will get started through learning a few fundamental features of R.

#Note that the # symbol is used to designate comments. Everything after that on a given
#line will be treated as commenting and not run.
#You can use # at the beginning of a line or after the code you want to run.
#It is good practice to comment your code smartly. Your collaborators, advisors, any future researchers who wish to build on your project, and your future self will thank you!
# In one of the next classes we will discuss what "smartly" means. In these scripts that I share with you, I often err on providing too many comments to explain the reasons behind the code. When you are coding, this will often be too many comments!

#If you wish to use the R console, you should read the below lines of commenting and then you can type in the lines of code (i.e. lines without a # at the beginning) into the R console and hit <Enter> to run that line of code. We will also begin using RStudio. In RStudio, we can easily load an entire commented script and run the lines of code.

#Topics for today: some basic syntax, assignment operator, the idea of functions and arguments.
#We will also create two types of data object, atomic vector and matrix. 
#We will try out the sample() function and a few other helpful functions and finish by performing arithmetic.

###PART 2 - VECTORS----

#First, we will generate integers from 1 to 20 using the following line of code:
1:20

#The above line of code resulted in the integers being printed to the screen. However, commonly, we want to assign the data to a data object to work with further. We use the assignment operator <- for this.
numbers <- 1:20

#Next time, once you have it installed, try running the above line again in RStudio.
#What do you notice in your environment window?

#Notice the use of spaces around the assignment operator.
#The code will work without the spaces (try it!), but it is good practice to use the spaces to make your code more readable. Readability of your code to humans is also a consideration!
#Please use good formatting practices for your own practice sessions and for your assignments.

#We can print the contents of our new object to the screen by simply typing numbers.
numbers

#Now, we want to do something with numbers. First, we will check what type of data we have created.
#The function class() will tell us what class of object we have. Try it.
class(numbers)

#What is returned is the word "integer". This is telling us that 
#we have created integer data. Our object numbers is actually a "vector" containing integer data.

#Vectors can be considered "contiguous cells containing data", i.e. a linear object.

#This quote is from the R Language Definition document, which can be accessed at the link below. I suggest it is very helpful to check this document to read the precise definitions of the language components that are new to you.
#https://cran.r-project.org/doc/manuals/r-release/R-lang.html

#And here is a bit more about vectors:
#https://cran.r-project.org/doc/manuals/r-release/R-lang.html#Vector-objects

#In the case of our numbers object, we have an "atomic" vector object (i.e. it contains just one type of element).i.e. Every element in the vector numbers is an integer.
#From the R language definition: "R has six basic ('atomic') vector types: logical [i.e. TRUE/FALSE], integer, real, complex, string (or character) and raw." For our Software Tools course, we will commonly be working with numeric data (real numbers or integers) and string data (e.g. DNA sequences and words) and sometimes logical data.

#Note that the storage mode for numeric data (real numbers) is double precision floating point. This is way of encoding a real number (decimal format) in binary. The "floating point" phrasing refers to the fact that there can a variable number of digits before and after the decimal point. Double precision floating point is a way of representing real numbers to at least to 15-16 significant digits. Beyond that, there can be small rounding errors. This level of precision is nearly always suitable, and note that it isn't possible to store infinite digits. We typically don't notice the minor rounding errors in daily practice as the level of precision is far beyond measurement error and scientifically relevant digits for most types of biological research, but it is worth being aware of this. You can occasionally get unexpected results due to this feature of data storage, e.g. when comparing whether two numbers are exactly equal. But most of the time, we will not notice this at all!

#A final note on vectors for now... While "atomic vectors" contain just one data of data, lists are "generic vectors" and can contain more than one type of element and can be more complex. Lists may even contain hierarchical elements (e.g. element 1 could be a number, element 2 could be a character string, element 3 could be a vector, element 4 could be a matrix, and element 5 could even be another list). So, lists are helpful for housing more complex types of data. We will learn more about lists later. For now, we will work with atomic vectors. Remember that atomic vectors can only contain ONE TYPE OF DATA.

#CODING CHALLENGE #1. Your turn! Create an atomic vector called y that contains the integers from 1 to 100.
y <- 1:100

?class

###PART 3 - GETTING STARTED USING FUNCTIONS----

#Functions in R can be thought of as verbs. They are the main way that we get R to **DO** something with our data. Here, we will briefly introduce the syntax for using a function in R that was already created by someone else. (In a future lesson, we will go over the formal definition of a function, and we will be creating our own functions to perform custom actions!)

#class() is a function, which we used above to check our data type.

#"Arguments" are specified inside the parentheses. We will use an analogy to a sentence in English to understand the syntax for running a basic function in R.

#If we think of the function name as the verb, arguments can be thought of as the object (noun upon which the verb acts) of the sentence in R. A given function may also take arguments that are analogous to adverbs (adverbs are words that describe how a particular verb is done). An example sentence in English:

#I throw the ball underhand.

#"I" is the subject of the above sentence (the one doing the verb). In an R "sentence", it is understood that you want R to do something. So, you don't have to specify that R is the subject of the sentence. So, the above sentence, rewritten in R syntax, might look like this:

#throw(object = ball, how = underhand)

#In this analogy, our function is called throw. It is the verb describing what we want to do. It is helpful when R function names indicate clearly what they do.

#The first argument for the function throw is called "object" and is specified above as ball. What object will I throw? Answer: ball. For many R functions, one of the arguments will be the data the function will act upon. Above, we passed our vector called numbers as an argument to the function class.

#In our analogy here, the argument called "how" is specified as being underhand. Arguments often have a controlled vocabulary. For example, in the above hypothetical function, the "how" argument might be defined as either "underhand" or "overhand". In other cases, an argument may need to be defined as a number (such as as a parameter value).

#How do we know what arguments a particular function takes? And, what are the options that are available to us? 

#****Always check the documentation**** when using a function for the first time or if you want to refresh yourself about the options available for a given function.

#In the documentation, default settings are indicated in capital letters. We don't have to type out the defaults if we want to use those settings.

#Have a look at the documentation for the function class() by running the following line:
?class()

#In the case of class(), it is a simple function with just one argument, called x. The function class() needs to act upon an R object. So, we can input any R object into class() to find out what type of object it is. It is always important to know what type of data we are working with, in case data reformatting or reassignment of data to a different data type or object class is needed for downstream analysis. This is very common.

#When we use arguments in the default order, we don't have to type the name of the argument. However, I would recommend to type the name of the argument in the case of a more complex function that takes a lot of arguments to make sure there are no errors relating to the ordering of the functions.

#Check if the below are the same:
class(numbers)
class(x = numbers)

#In the first instance, we just gave numbers as the R object we want to check. There is no problem with this, as we are specifying the arguments in the default order. In this case, the function takes just one argument, and there is no ambiguity. The second line has the same meaning, and saves us typing.

#Tip: It is often very helpful to look at the examples at the bottom of the documentation page to see how that function may be used. Reading R documentation can take some getting used to if you are new to R. We will go through a few documentation pages together over the coming weeks.

#Now take a look at the documentation for a function with more optional arguments.
?sample()

#NOTE: If you just type the function name, you will get the script for the function, rather than the help. So, you need the question mark if you want to access the documentation.

#It can be helpful to look at the source script, such as if you want to see how the function works or if you want to adapt the code in a different way.
sample

#We will take a random sample from among the integers from 1 to 20 and assign the results to a vector called sampleA.
sampleA <- sample(1:20)

#Let's look at our sample.
samplea

#Ooops. What happened? Note that R is case sensitive!
sampleA

#That's better.

#Try creating a new sample by running the above line again. Note you can use the up arrow to scroll through lines of code you recently ran. Do you get the same result or different? 

#What does the below do? Is this the same or different as when leaving the replace argument set to the default?
sampleB <- sample(1:20, replace = TRUE)
sampleB

#Sampling with and without replacement of the elements might be done in different contexts. We will return to that issue in a later lesson (such as when using random resampling to estimate the diversity in sample).

#Try making a simple plot. What does this do? We will return to the topic of plots many times, but it is fun to make our first plot in R!
plot(sampleA, sampleB)

#Practise having a look at the documentation to see what arguments the function plot() takes.
?plot

#CODING CHALLENGE #2: Read the documentation from the function sample(). Through using the "size" argument, take a sample of 20 integers between 1 and 100 (without replacement, which is the default) from the vector called y that you created above.

###PART 4 - ARITHMETIC----

#Now moving on... we can use R to perform basic calculations. Run the following:
3 + 2
3 - 2
3^2
3 * 3
log(10)
log10(10)

#We can use parentheses to govern the order of arithmetic operations.
(3 + 2) * 10
3 + (2 * 10)

#We can then scale up to repeat operations across many values. Here, the number 2 is recycled across the elements of x. We are using our vector of integers created above, called numbers.
numbers + 2

#Here, we are adding the elements in two vectors together. Check. What did that do?
sampleA + sampleB

#What happens when you run these?
numbers + 1:5

#The integers 1:5 get recycled through the length of the vector numbers.

#That reminds me that length() is another extremely commonly used function available through base R. Try this:
length(numbers)
?length

#Tip 1: Try playing with variations of the above lines of code in your own way to get used to the different functions.

#Tip 2: Try using R for various tasks you have to do, to get coding practice and have fun!
#Can't decide what to do this Saturday night? We can create a vector of choices. (Note the usage of quotation marks to designate character data.)
choices <- c("swirl.tutorial", "get.take.out", "watch.movie", "rest")
choices
sample(choices, 1)

#CODING CHALLENGE #3: Create a line of R code for a simple decision you need to make. e.g. Will you have coffee or tea? You need to choose 2 of 10 books to take on a trip... write a line of code to do it for you, and show us next class!

###PART 5 - CHARACTER DATA AND INTRO TO MATRIX----

#In this section, we will introduce another data type, further practise using functions, and introduce the matrix.

#Now we are going to create a different type of data object, also an atomic vector but with a different type of element.
#Note the quotation marks. What happens if you leave those out?
class_names <- c("Jane", "Xin", "Joe", "Bob", "Yves", "Emma")

#Checking class of the object class_names
class(class_names)

#What kind of element does the atomic vector class_names contain?

#Above, we used the concatenate or combine function c() to combine the character elements into a vector.

#c() "is a generic function which combines its arguments". (from the R help)
#"The default method combines its arguments to form a vector." (from the help)
?c()

#This hypothetical class will have a group project. Let's use R to make groups.

#We can use the sample() function again.

#This time, we are adding an optional argument, to take a sample of size 3 from the vector class_names. Note we can just specify our argument settings, if we use the arguments in default order.
sample(class_names, 3)

#Note we can also specify the argument names. That is best practice to be explicit. This is especially important when there are a lot of arguments.
sample(x = class_names, size = 3)

#Do you predict the the following line will work? Why or why not? Make your prediction BEFORE running this line.
sample(3, class_names)

#Do you predict the following line will work as expected? Why or why not?
sample(size = 3, x = class_names)

#You can try different numbers for size. Play with that now.

#What happens if you take out that argument altogether to specify a sample size?
sample(class_names)

#So, what have we learned? What is the default for the argument size? Cross check your observations against the documentation.
?sample()

#Random sampling is very commonly used in diverse research settings, such as study design.

#Imagine you want to take 6 samples randomly along a transect.
#You have a central starting point in this case, and you are interested in sampling flowers up to 10 m in each direction.
#You could select random starting points along the transect.
#The below code would generate 6 random starting points between -10 and +10.
sample(-10:10, 6)

#So, notice that further above we used a named object, class_names, as an argument for the sample() function.
#Here, instead we specified the data directly. -10:10 gives us a vector of integers from -10 to 10.
#R often provides you with multiple ways to do things. This is just one example.
#Sometimes, choosing one over another is a matter of personal preference.
#In other cases, there are better and worse ways to do things, considering clarity,
#reproducibility, and efficiency.

#Run the above line again. Do you get the same or different numbers?
sample(-10:10, 6)

#If, for our purposes, we want to create a perfectly reproducible script, we can rectify the above issue by using the function set.seed().

#Run this:
set.seed(203)

#Now run the sampling function multiple times, using the same setting for set.seed() each time. 
sample(-10:10, 6)

set.seed(203)
sample(-10:10, 6)

set.seed(203)
sample(-10:10, 6)

#What do you observe?

#The sample() function can be useful in many contexts.

#Imagine we have a class with 20 students (as in my Arctic Ecology field course).

#We could randomly create 5 groups, each containing 4 students, very easily.

#Here, we will enter our hypothetical students:
students <- c("Ann", "Xin", "Jo", "Bev", "Adam", "Jose", "Maria", "David", "Sal", "Lea", "Bo", "Aria", "Evan", "Tia", "Eric", "Cris", "Deb", "Pat", "Prya", "Rose")
students
class(students)
length(students)

#Using the function matrix(), we will create a new data object, of the class matrix.
#We will do this by sampling from the vector called students. Our matrix will have 5 rows.
#By reading across the rows, we will see which students are in which groups.
#A matrix is a two-dimensional, rectangular arrangement of elements.
groups1 <- matrix(sample(students), nrow = 5)
groups1

#Question: In this case, would we want to sample with or without replacement? You can try both. What happens?

#Reminder about viewing the help.
?matrix

#Note we didn't use all the available arguments. Often, the defaults are used.
#However, it is wise always to look up the defaults.
#If you are happy with the default setting for a given argument for what you want to do, then you don't have to type it out.
#For example, I was fine with random sampling being used to fill the matrix by columns.
#Next time we will add another argument to give names to the groups.

#Also, try out the function is.matrix() on some of the data objects we have created so far.
is.matrix(groups1)
is.matrix(numbers)
is.matrix(class_names)

#Note we can use a period or underscore in object names. Don't use other characters.
#Also, R is case-sensitive. For example, try this:
is.matrix(Groups1)

#Try this. What is it telling us about our matrix called groups1? This is a very helpful function.
dim(groups1)
?dim

#Play around with this example code on your own, and bring questions with you to class.

#See you next time!

