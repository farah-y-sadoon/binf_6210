#### QUIZ 1 

soil.nitrogen <- c(45, 30, 20, 16, 18, 32, 48, 33)
plant.height <- c(50, 42, 33, 5, 8, 35, 48, 31)

class(soil.nitrogen)

plant.model <- lm(plant.height ~ soil.nitrogen)
summary(plant.model)

Plant.Model


x <- c(1, 3, 5, 7, 9) * c(10, 10, 10, 10, 10)
x

x <- c(1, 3, 5, 7, 9) * 10
x

my.vector <- c(1, 3, 5, 7, 9)

x <- my.vector * 10
x
