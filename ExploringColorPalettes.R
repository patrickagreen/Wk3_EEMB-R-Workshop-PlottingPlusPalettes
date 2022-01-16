library(tidyverse)

# create a random number table 
randomTable <- sample(500,50)
dim(randomTable) = c(10,5)
randomTable
class(randomTable)

colnames(randomTable) <- c("")