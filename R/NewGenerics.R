#'@export
coords <- function(x, new.coords) UseMethod("coords", x)

#'@export
`coords<-` <- function(x, ..., value) UseMethod("coords<-", x)

#'@export
geoArea <- function(x) UseMethod("geoArea", x)


#'@export
maxDist <- function(x) UseMethod("maxDist", x)
