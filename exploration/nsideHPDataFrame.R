

#'minDist2nside
#'
#'Get an nside such that all points belong
#'to different pixels, when the minimum
#'distance between the points is \code{dist}.
#'
#'We treat pixels as circles. From the total surface
#'area of a unit sphere and the radius of a circle with
#'area 4*pi/(12*nside^2), we derive a sufficiently
#'large nside to achieve the desired separation.
#'
#'@param dist The minimum distance between any two points in a data.frame
#'of points that lie on S^2
#'
#'@keywords internal
#'
#'@export
minDist2nside <- function(dist)
{
  n <- 2/(dist*sqrt(3))

  # Round to the next power of 2
  return(2^ceiling(log2(n)))
}


library(rcosmo)
sky <- CMBDataFrame(nside = 16, ordering = "nested")


