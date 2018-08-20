#'@export
coords <- function(x, ...) UseMethod("coords", x)

#'@export
`coords<-` <- function(x, ..., value) UseMethod("coords<-", x)

#'@export
geoArea <- function(x) UseMethod("geoArea", x)

#'@export
maxDist <- function(x) UseMethod("maxDist", x)

#'@export
pix <- function(x, ...) UseMethod("pix", x)

#'@export
`pix<-` <- function(x, ..., value) UseMethod("pix<-", x)

#'@export
nside <- function(x) UseMethod("nside", x)

#'@export
ordering <- function(x, ...) UseMethod("ordering", x)

#'@export
`ordering<-` <- function(x, ..., value) UseMethod("ordering<-", x)

#'@export
window <- function(x, ...) UseMethod("window", x)

#'@export
`window<-` <- function(x, ..., value) UseMethod("window<-", x)
