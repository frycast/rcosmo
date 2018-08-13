#' Geodesic distance on the unit sphere
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





# NEXT: geoAngle

