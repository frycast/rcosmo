# Compute the empirical covariance between observations in a CMBDataFrame
# assuming homogeneity and isotropy.
# INPUTS: CMBDataFrame, vector r of radii
# OUTPUTS: C(r), for each r


#' Covariance for CMB
#'
#' This function provides an empirical covariance estimate for data
#' in a CMBDataFrame or data.frame. It places data into bins.
#'
#' @param cmbdf is a CMBDataFrame or data.frame
#' @param num.bins specifies the number of bins
#' @param sample.size optionally specify the size of a simple random
#' sample to take before calculating covariance. This may be useful if
#' the full covariance computation is too slow.
#' @param max.dist an optional number between 0 and pi specifying the
#' maximum geodesic distance between any two points in cmbdf.
#' For example, if cmbdf represents a full sky map
#' or a random sample of a full sky map then \code{max.dist = pi}.
#' If max.dist is known then specifying it may reduce
#' computation time.
#'
#' @return
#'
#' @examples
#'
#' @export
covCMB <- function(cmbdf,
                   num.bins = 10,
                   sample.size,
                   max.dist)
{

  if ( class(cmbdf) != "CMBDataFrame" ) {

    stop("cmbdf must be a CMBDataFrame")

  }

  if ( !missing(sample.size) ) {

    cmbdf <- sampleCMB(cmbdf, sample.size = sample.size)

  }

  if ( coords(cmbdf) != "cartesian" ) {

    coords(cmbdf) <- "cartesian"

  }

  if ( !all.equal(names(cmbdf), c("x", "y", "z", "I")) )
  {
    stop("cmbdf must have columns named 'x', 'y', 'z', 'I' in that order")
  }

  if (missing(max.dist)) {

    return(covCMB_internal2(cmbdf, num.bins))

  } else {

    breaks <- seq(0, max.dist, length.out = num.bins+1)
    return(covCMB_internal1(cmbdf, breaks))

  }

}




## OLD VERSION (DISUSED)
# covCMB <- function(cmbdf,
#                    num.bins = 10,
#                    sample.size,
#                    large.sample,
#                    max.dist)
# {
#
#   if ( !missing(sample.size) ) {
#
#     cmbdf <- sampleCMB(cmbdf, sample.size)
#
#   }
#
#   if ( coords(cmbdf) != "cartesian" ) {
#
#     coords(cmbdf) <- "cartesian"
#
#   }
#
#   if ( missing(large.sample) ) {
#
#     large.sample <- ifelse( missing(sample.size),
#                             nrow(cmbdf) > 10000,
#                             sample.size > 10000)
#
#   }
#
#   if ( large.sample == TRUE ) {
#
#     return( covCMB_big(cmbdf) )
#
#   } else if ( large.sample == FALSE ) {
#
#     return( covCMB_small(cmbdf) )
#
#   } else {
#
#     stop("large.sample must be TRUE or FALSE")
#
#   }
# }
#
# covCMB_big <- function( cmbdf )
# {
#   cat("Using covCMB_big...\n")
#   print(cmbdf)
# }
#
# covCMB_small <- function( cmbdf )
# {
#   cat("Using covCMB_small...\n")
#
#   if (missing(max.dist))
#   {
#     return(covCMB_internal2(cmbdf, num.bins))
#   }
#   else
#   {
#     cut.points <- seq(0, max.dist, length.out = num.bins+1)
#     return(covCMB_internal1(cmbdf, cut.points))
#   }
#
#
#   return(dists)
# }


