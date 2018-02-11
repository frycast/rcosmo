#' geodesic distance on the unit sphere
#'
#'@param p1 a point on the unit sphere given in cartesian coordinates (x,y,z)
#'@param p2 a point on the unit sphere given in cartesian coordinates (x,y,z)
#'
#'@return the geodesic distance between \code{p1} and \code{p2}
#'
#'@export
geoDist <- function(p1,p2) {
  return(acos(sum(p1*p2)))
}




#' Use Haversine formula
#'
#' Uses the Haversine formula to give the
#' geodesic distance between two points on the unit sphere given
#' in latitude and longitude. The Haversine formula is favoured
#' for its numerical stability
#'
#'
#'@param p1 a 2 element vector (lat, long) specifying a point on the
#'unit sphere
#'@param p2 a 2 element vector (lat, long) specifying a point on the
#'unit sphere
#'
#'@return the geodesic distance between \code{p1} and \code{p2}
#'
#'@export
haversineDist <- function(p1,p2) {

  dlat <- abs((p1[1] - p2[1])/2)
  dlon <- abs((p1[2] - p2[2])/2)
  return(2*asin(sqrt( (sin(dlat))^2 + cos(p1[1])*cos(p2[1])*(sin(dlon))^2 )))
}





#' Convert from spherical to cartesian
#'
#' @param p a point in latitude and longitude \code{c(lat, long)} or
#' a data.frame containing columns lat and long
#'
#' @return the cartesian coordinates of p
#'
#' @export
sph2car <- function(p)
{
  if ( is.data.frame(p) ) {

    df <- data.frame(lat = p$lat, long = p$long)
    car <- apply(df, 1, sph2car_helper)
    car <- as.data.frame(t(car))
    names(car) <- c("x", "y", "z")

    return(car)

  } else {
    return(sph2car_helper(p))
  }
}

sph2car_helper <- function(p)
{
  theta <- pi/2 - p[1]
  return(c(sin(theta)*cos(p[2]),
           sin(theta)*sin(p[2]),
           cos(theta)))
}



## HELPER FUNCTION, CROSS PRODUCT
vector_cross <- function(a, b) {
  if(length(a)!=3 || length(b)!=3){
    stop("Cross product is only defined for 3D vectors.");
  }
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (a[i1]*b[i2] - a[i2]*b[i1])
}







# HELPER FUNCTION: RODRIGUES ROTATION FORMULA SEE WIKIPEDIA PAGE:
# Rotation axis k is defined by being away from a and towards b.
# This function rotates a directly to b (rotating whole sphere points p_xyz).
rodrigues <- function(a,b,p.xyz)
{
  norm.a <- sqrt(sum(a * a))
  norm.b <- sqrt(sum(b * b))
  if ( !isTRUE( all.equal(a/norm.a, b/norm.b, check.attributes = FALSE, use.names = FALSE) ) )
  {
    k <- vector_cross(a,b)
    k <- k/sqrt(sum(k^2)) # normalised k
    K <- matrix( c(  0   , -k[3],  k[2],
                     k[3], 0    , -k[1],
                    -k[2], k[1] ,  0  ), nrow = 3, byrow = TRUE)
    theta <- acos( sum(a*b) / ( norm.a*norm.b ) ) # The angle between a and b
    I <- diag(c(1,1,1))
    R <- I + sin(theta)*K + (1-cos(theta))*K%*%K #Rodrigues' Formula.
    p.xyz <- t(R%*%t(p.xyz))
  }

  p.xyz
}
