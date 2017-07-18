#' Plot CMB Data
#'
#' This function produces a plot from a CMB Data Frame.
#'
#'@param CMBDataFrame a CMB Data Frame with either spherical or cartesian coordinates.
#'
#'@return
#'A plot of the CMB data.
#'
#'@examples
#' ##This is a place holder.
#'
#'@export
plotCMB <- function(CMBDataFrame) {
  N <- nrow(CMBDataFrame)
  coords <- attr(CMBDataFrame,"coords")
  
  try(if(coords != "spherical" && coords != "cartesian") stop("Coordinates must be spherical or cartesian"))
  
  if (coords == "spherical") {
    sm <- matrix(c(CMBDataFrame$long, CMBDataFrame$lat, rep(1,N)), nrow = N)
    smx <- sph2car(sm, deg = FALSE)
  } else {
  # Else coords are cartesian
    smx <- matrix(c(CMBDataFrame$x, CMBDataFrame$y, CMBDataFrame$z))
  }
  
  cols <- colmap[cut(CMBDataFrame$I,length(colmap))]
  
}