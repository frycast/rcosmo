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
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sampleSize = 800000)
#' plotCMB(df)
#'
#'@export
plotCMB <- function(CMBDataFrame) {
  N <- nrow(CMBDataFrame)
  coords <- attr(CMBDataFrame,"coords")

  try(if(coords != "spherical" && coords != "cartesian") stop("Coordinates must be spherical or cartesian"))

  if (coords == "spherical") {
    sm <- matrix(c(CMBDataFrame$long, CMBDataFrame$lat, rep(1,N)), nrow = N)
    smx <- sphereplot::sph2car(sm, deg = FALSE)
  } else {
  # Else coords are cartesian
    smx <- matrix(c(CMBDataFrame$x, CMBDataFrame$y, CMBDataFrame$z))
  }

  cols <- colmap[cut(CMBDataFrame$I,length(colmap))]
  rgl::open3d()
  rgl::bg3d("black")
  rgl::plot3d(smx, col = cols, type = "p", cex = 5, pch = 3, box = FALSE, axes = FALSE)
}
