#' Geodesic distance on the unit sphere
#'
#'@param p1 a 3 element vector on the unit sphere
#'given in Cartesian coordinates (x,y,z),
#'or a 2 element vector (theta, phi) giving spherical coordinates,
#'can also be a \code{\link{data.frame}} or matrix with rows specifying vectors
#'@param p2 a vector on the unit sphere given in Cartesian coordinates (x,y,z)
#'or a named vector (theta, phi) giving spherical coordinates,
#'can also be a \code{\link{data.frame}} or matrix with rows specifying vectors
#'
#'@return The geodesic distance between \code{p1} and \code{p2}
#'
#'@export
geoDist <- function(p1,p2) {

  if (is.matrix(p1)) p1 <- as.data.frame(p1)
  if (is.matrix(p2)) p2 <- as.data.frame(p2)


  p1 <- convertToXYZDataFrame(p1)
  p2 <- convertToXYZDataFrame(p2)

  if ( nrow(p1) != nrow(p2) )
  {
    stop("p1 and p2 must have the same number of rows")
  }


  return(acos(sapply(c(p1$x*p2$x + p1$y*p2$y + p1$z*p2$z),
                     function(x) {max(-1,min(x,1))})))
}


# Helper function for geoDist
convertToXYZDataFrame <- function(p)
{
  if ( (is.null(ncol(p)) && length(p) == 2)
       || ncol(p) == 2 )
  {
    names(p) <- c("theta","phi")
    return(rcosmo::sph2car(data.frame(theta = p["theta"], phi = p["phi"])))
  }
  else if ( (is.null(ncol(p)) && length(p) == 3)
             || ncol(p) == 3 )
  {
    names(p) <- c("x","y","z")
    return(data.frame(x = p["x"], x = p["y"], x = p["z"]))
  }
  else
  {
    stop("points have have length 2 (theta, phi) or 3 (x,y,z)")
  }

}




# NEXT: geoAngle


a <- c(1,2,3)
p <- data.frame(theta = c(0,1,1), phi = c(0,0,1))
p1 <- c(theta = 1, phi = 1)
