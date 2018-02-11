#' Nest Search
#'
#' Finds the closest HEALPix pixel center to a given \code{target} point,
#' specified in cartesian coordinates, using an efficient nested search
#' algorithm. HEALPix indices are all assumed to be in the "nested"
#' ordering scheme.
#'
#' @param target is a vector of Cartesian coordinates for the target
#' point on S^2
#' @param Nside is the Nside for which the HEALPix points are searched for
#'
#' @return if \code{index_only} TRUE then the output will be a HEALPix index.
#' If \code{index_only} FALSE then the output is the list containing the HEALPix index
#' and Cartesian coordinate vector of the HEALPix point closest to \code{tp}.
#'
#' @examples
#' # Find the pix index and Cartesian coordinates of the HEALPix point
#' # at Nside closest to the target point c(0,0,1)
#' h <- nestSearch(c(0,0,1),Nside=1024,index_only=FALSE, plot_points=TRUE )
#' cat("Closest HEALPix point to (0,0,1) at Nside = 1024 is (",h$xyz,")")
#'
#' @export
nestSearch <- function(target = c(0,0,1), Nside = 16,
                       index.only = TRUE,
                       j = c(min(3,log2(Nside)),
                             min(7,log2(Nside)),
                             log2(Nside)) )
{
  j <- c(0,j)
  for ( i in 2:length(j) )
  {
    h <- nest_search( target, j2=j[i], j1=j[i-1], pix.j1 = h$pix)
  }

  if (index.only==TRUE) {
    return(h$pix)
  }
  else {
    return(list( xyz=h$xyz, pix=h$pix ))
  }
}
