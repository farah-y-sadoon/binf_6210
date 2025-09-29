#### IN CLASS ACTIVITY BINF 6210 - SOFTWARE TOOLS

### 1. Open a blank script in RStudio. Then, so that we will get consistent results across runs and computers, run:
set.seed(123)

## Next, write a line of code that creates a matrix called y. This matrix includes 100 values sampled from the random normal distribution, with mean of 10 and standard deviation of 5. The matrix has 20 rows and 5 columns
y <- matrix((rnorm(100, mean = 10, sd = 5)), nrow = 20, ncol = 5)

### 2. Extract columns 3 and 4 using indexing by position
(cols_3_4 <- y[, 3:4])

### 3. Calculate the mean of the entire matrix y
(mean_y <- mean(y))


### 4. Calculate the mean for row 2 (in one line of code)
(mean_col_2 <- mean(y[2]))

### 5. Run the same code as for question #1 again. Set the seed again and generate a matrix with the same parameters, this time called z. Then, check if the two matrices are equal.

set.seed(49) # different seed from the first question to generate different values
z <- matrix((rnorm(100, mean = 10, sd = 5)), nrow = 20, ncol = 5)

set.seed(123) # same seed as the first question to generate the same values
z <- matrix((rnorm(100, mean = 10, sd = 5)), nrow = 20, ncol = 5)

identical(y, z)

### 6. Open example Script #5 and run up to line 132 (i.e. loading "Cnidaria_BOLD_data.tsv" file and assigning to data frame named dfBOLD). Check what class of data is the column called “lat”, containing latitude data, using indexing by name.
class(df_bold$lat)

### 7. Create a new data frame, called dfBOLDnuc, containing columns 1, 2, 14, 47, and 72 from dfBOLD, using indexing by position. (Hint: remember to check this worked. You can look in the viewer to look at your new data frame. You should have the same number of rows as in the original and the correct five columns.)
(df_BOLDnuc <- df_bold[, c(1, 2, 14, 47, 72)])

### 8. Using dfBOLDnuc, what are the unique orders present in this dataset? (You can use indexing by position or name)
(unique_orders <- unique(df_BOLDnuc$order_name))

### 9. How many records does each order have?
(count_order_records <- table(df_BOLDnuc$order_name))

### 10. Find the name of the order with the most records.
(most_records <- sort(table(df_BOLDnuc$order_name), decreasing = TRUE))

### 11. What is the mean latitude for records in this dataset? Hints:
## Remember that a single column of a data frame is an atomic vector. So, you can pass numerical atomic vectors from data frames to mathematical functions.
## Note this dataset contains some missing data. So, set the argument called na.rm to TRUE for the function mean().
(mean_lat <- mean(df_BOLDnuc$lat, na.rm = TRUE))

### 12. A little trickier. Using dfBOLDnuc, what is the mean latitude among records for the order with the most data? Hints:
## You can use indexing by condition
## Exactly equals to is specified as ==
## Remember you need quotation marks for character data
## Recall that if you get stuck, look at the V2 sheet (with example answers) and then work through it to ensure you understand it. Then, try the next question on your own before looking at the answer.
(mean_lat_most_records <- mean(df_BOLDnuc$lat[df_BOLDnuc$order_name == "Alcyonacea"], na.rm = TRUE))

### 13. Let’s keep practising using additional examples.
## A) What is the maximum latitude in the entire dataset?
(max_lat <- max(df_BOLDnuc$lat, na.rm = TRUE))

## B) And, what is the standard deviation of latitude for records in the order Scleractinia? (Hint: Keep setting na.rm to TRUE.)
(sd_lat <- sd(df_BOLDnuc$lat, na.rm = TRUE))
