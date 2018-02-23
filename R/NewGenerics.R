#'@export
coords <- function(x, new.coords) UseMethod("coords", x)

#'@export
`coords<-` <- function(x, ..., value) UseMethod("coords<-", x)

#'@export
area <- function(x) UseMethod("area", x)
