#library(pryr)

# Display all generics that have a method for class 'CMBDataFrame' (none at the time of writing)
methods(class="CMBDataFrame")

# Display all methods for the generic 'print' (many, corresponding to different packages)
methods("print")

# Define CMBDF methods for some generics
print.CMBDataFrame <- function(cmbdf) cat("This is a CMBDataFrame")
mean.CMBDataFrame <- function(cmbdf) cat("This is the mean of a CMBDF")
summary.CMBDataFrame <- function(cmbdf) cat("This is the summary of a CMBDF")
plot.CMBDataFrame <- function(cmbdf) cat("This is a plot of a CMBDF")

### Note:
# The var function is provided by package {stats} and is not a generic.
# Hence we should avoid overriding var in rcosmo.
###

# Construct a CMBDF object
cmbdf <- CMBDataFrame(CMBData = '../CMB_map_smica1024.fits', sampleSize=100)

# Call our generics on the CMBDF object
print(cmbdf)
mean(cmbdf)
summary(cmbdf)
plot(cmbdf)
cmbdf #this just calls print(cmbdf)
print.data.frame(cmbdf) #calls the data.frame method of the plot generic

# Again, Display all generics that have a method for class 'CMBDataFrame' (now we have two: mean & print)
methods(class = "CMBDataFrame")


# Using replacement functions ---------------------------------------------

coords <- function(x) {
  attr(cmbdf, "coords")
}

`coords<-` <- function(cmbdf,...,value) {
  cat("The coords were: ", attr(cmbdf, "coords"), "\n", sep = "")
  attr(cmbdf, "coords") <- value
  cat("The coords are now: ", attr(cmbdf, "coords"), "\n", sep = "")
  cat("This code could actually do the transformation")
  cmbdf
}





# Understanding the role of lists in R ------------------------------------

library(listviewer)

# Using str versus using listviewer
str(cmbdf)
listviewer::jsonedit(cmbdf, mode = "view")


# Usage of formulas -------------------------------------------

cars$qspeed <- cut(cars$speed, breaks = quantile(cars$speed),
                   include.lowest = TRUE)
is.factor(cars$qspeed)

plot(dist ~ speed, data = cars)
plot(dist ~ speed, data = cars)


# New style (S4) classes -------------------------------------------------------

setClass("Person",
         slots = list(name = "character", age = "numeric"))
setClass("Employee",
         slots = list(boss = "Person"),
         contains = "Person")

getClass("Person")

alice <- new("Person", name = "Alice", age = 40)
bob <- new("Employee", name = "Bob", age = 20, boss = alice)

alice@age
bob@boss

# Create a class which inherits from a base types and S3 objects
setClass("RangedNumeric",
         contains = "numeric",
         slots = list(min = "numeric", max = "numeric"))

rn <- new("RangedNumeric", 1:10, min = 1, max = 10)
rn@min
rn@.Data # The .Data slot contains the base type or S3 object

# Define some generic methods
getGeneric("mean") # The standard S4 generic 'mean' already exists

setMethod("mean",
          c(x = "Person"),
          function(x) {
            cat("This is the mean of a Person: ", x@name)
          }
)

setMethod("mean",
          c(x = "Employee"),
          function(x) {
            cat("This is the mean of an Employee: ", x@name)
          }
)

# Print is an S3 generic
print.Person <- function(p) cat("This is the S3 print of a Person: ",
                                p@name)

# Print is also an s4 generic
getGeneric("print")
setMethod("print",
          c(x = "Person"),
          function(x) {
            cat("This is the S4 print of a Person: ", x@name)
          }
)

# Show is also a standard S4 generic
getGeneric("show")
setMethod("show",
          "Person",
          function(object = "Person") {
            cat("Showing ", object@name, "\n")
            print(object)
          }
)

# Var is not a standard s4 generic so let's define it
getGeneric("var")
setGeneric("var", function(x) {
  standardGeneric("var")
})
setMethod("var",
          "Person",
          function(x = "Person") {
            cat("This is the variance of a Person")
          }
)
var(alice)
var(bob)

# We can make methods that use two arguments
getGeneric("summary")
setMethod("Summary",
          c("Person","Employee"),
          function(object = "Person", y = "Employee") {
            cat("Summary of Person and Employee")
          }
)

