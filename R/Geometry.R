#' Geodesic distance on the unit sphere
#'
#' Get geodesic distance between points on the unit sphere
#'
#'@param p1 A \code{\link{data.frame}} with rows
#'specifying numeric points located on the unit sphere.
#'It should have columns labelled x,y,z
#'for Cartesian or theta, phi for spherical colatitude and
#'longitude respectively.
#'@param p2 Same as p1.
#'@param include.names Boolean. If TRUE then the row and
#'column names of the returned matrix will be taken from
#'the points in \code{p1} and \code{p2} (see examples
#'below).
#'
#'@return Let \eqn{n} denote the number of rows of \code{p1}
#'and let \eqn{m} denote the number of rows of \code{p2}.
#'Then the returned object is an \eqn{n} by \eqn{m} matrix
#'whose entry in position \eqn{ij} is the geodesic distance
#'from the \eqn{i}th row of \code{p1} to the
#'\eqn{j}th row of \code{p2}.
#'
#'@examples
#'
#' p1 <- data.frame(diag(3))
#' colnames(p1) <- c("x", "y", "z")
#' p1
#' p2 <- data.frame(x=c(1,0), y=c(0,3/5), z=c(0,4/5))
#' p2
#' geoDist(p1, p2, include.names = FALSE)
#'
#'@export
geoDist <- function(p1,p2, include.names = FALSE) {

  p1 <- as.data.frame(p1)
  p2 <- as.data.frame(p2)

  if ( include.names )
  {
    n1 <- paste(round(as.data.frame(t(p1)), digits = 2))
    n1 <- substring(n1, 2)
    n2 <- paste(round(as.data.frame(t(p2)), digits = 2))
    n2 <- substring(n2, 2)
  }

  p1 <- rcosmo::coords(p1, new.coords = "cartesian")
  p2 <- rcosmo::coords(p2, new.coords = "cartesian")

  # Dot product matrix using helper function pf1
  dp <- t(apply(p1, MARGIN = 1, p1f, p2 = p2))
  dp <- CleanFPErrors(dp)


  if (include.names)
  {
    rownames(dp) <- n1
    colnames(dp) <- n2
  }

  return(acos(dp))
}

## Helper function 1 for geoDist
CleanFPErrors <- function(x) {pmin(pmax(x,-1.0),1.0)}

## Helper function 2 for geoDist
p1f <- function(p1, p2) {
  apply(p2, MARGIN = 1, function(p2) {
    sum(p1*p2)
  })
}


#'Get the minimum geodesic distance between points
#'
#'
#'Get the minimum geodesic distance either between all points in a data.frame
#'pairwise, or between all points in a data.frame and one target
#'point.
#'
#'@param df A \code{data.frame} with columns x,y,z for cartesian
#'or theta, phi for spherical colatitude and longitude respectively.
#'The rows must correspond to points on the unit sphere.
#'If this is a \code{\link{HPDataFrame}} or \code{\link{CMBDataFrame}}
#'and coordinate columns are missing, then coordinates will be
#'assigned based on HEALPix pixel indices.
#'@param point An optional target point on the unit sphere
#'in cartesian coordinates, in which case all distances are
#'calculated between \code{point} and the points in \code{df}.
#'
#'@return If \code{point} is specified: the shortest
#'distance from \code{point} to the
#'points specified by the rows of \code{df}.
#'If \code{point} is not specified: the shortest
#'distance pairwise between points in \code{df}.
#'
#'@name minDist
#'
#'@examples
#'
#' ## Using a CMBDataFrame with HEALPix coordinates only
#' cmbdf <- CMBDataFrame(nside = 1, spix = c(1,5,12), ordering = "ring")
#' plot(cmbdf, hp.boundaries = 1, col = "blue", size = 5)
#' p <- c(0,0,1)
#' minDist(cmbdf, p) # no need to have coordinates
#'
#' ## Using a HPDataFrame with HEALPix coordinates only
#' hp <- HPDataFrame(nside = 1, I = rep(0,3), spix = c(1,5,12) )
#' minDist(hp, p) # notice no need to have coordinates
#'
#' ## Using a data.frame with cartesian coordinates
#' coords(hp) <- "cartesian"
#' df <- data.frame(x = hp$x, y = hp$y, z = hp$z)
#' minDist(df, p)
#'
#' ## Using a data.frame with spherical coordinates
#' coords(hp) <- "spherical"
#' df <- data.frame(theta = hp$theta, phi = hp$phi)
#' minDist(df, p)
#'
#' ## min distance between points in cmdf
#' minDist(cmbdf)
#'
#'
#'@export
minDist <- function(df, point)
{
  df <- coords(df, new.coords = "cartesian")

  if (!missing(point)) {

    minDist_internal1(df[,c("x","y","z")], point)
  } else {

    minDist_internal2(df[,c("x","y","z")])
  }
}




#'Get the maximum geodesic distance between points
#'
#'Get the maximum geodesic distance either between all points in a data.frame
#'pairwise, or between all points in a data.frame and one target
#'point.
#'
#'@param df A \code{data.frame} with columns x,y,z for cartesian
#'or theta, phi for spherical colatitude and longitude respectively.
#'The rows must correspond to points on the unit sphere.
#'If this is a \code{\link{HPDataFrame}} or \code{\link{CMBDataFrame}}
#'and coordinate columns are missing, then coordinates will be
#'assigned based on HEALPix pixel indices.
#'@param point An optional target point on the unit sphere
#'in cartesian coordinates, in which case all distances are
#'calculated between \code{point} and the points in \code{df}.
#'
#'@return If \code{point} is specified: the longest geodesic
#'distance from \code{point} to the
#'points specified by the rows of \code{df}.
#'If \code{point} is not specified: the longest geodesic
#'distance pairwise between points in \code{df}.
#'
#'@name maxDist
#'
#'@examples
#'
#' ## Using a CMBDataFrame with HEALPix coordinates only
#' cmbdf <- CMBDataFrame(nside = 1, spix = c(1,5,12), ordering = "ring")
#' plot(cmbdf, hp.boundaries = 1, col = "blue", size = 5)
#' p <- c(0,0,1)
#' maxDist(cmbdf, p) # no need to have coordinates
#'
#' ## Using a HPDataFrame with HEALPix coordinates only
#' hp <- HPDataFrame(nside = 1, I = rep(0,3), spix = c(1,5,12) )
#' maxDist(hp, p) # notice no need to have coordinates
#'
#' ## Using a data.frame with cartesian coordinates
#' coords(hp) <- "cartesian"
#' df <- data.frame(x = hp$x, y = hp$y, z = hp$z)
#' maxDist(df, p)
#'
#' ## Using a data.frame with spherical coordinates
#' coords(hp) <- "spherical"
#' df <- data.frame(theta = hp$theta, phi = hp$phi)
#' maxDist(df, p)
#'
#' ## max distance between points in cmdf
#' maxDist(cmbdf)
#'
#'@export
maxDist <- function(df, point)
{
  df <- coords(df, new.coords = "cartesian")

  if (!missing(point)) {

    maxDist_internal1(df[,c("x","y","z")], point)
  } else {

    maxDist_internal2(df[,c("x","y","z")])
  }
}


# NEXT: geoAngle

