# Compute the empirical covariance between observations in a CMBDataFrame
# assuming homogeneity and isotropy.
# INPUTS: CMBDataFrame, vector r of radii
# OUTPUTS: C(r), for each r


#' Covariance for CMB
#'
#' This function provides an empirical covariance estimate for data
#' in a CMBDataFrame or data.frame. It places data into bins.
#'
#' @param cmbdf is a \code{\link{CMBDataFrame}} or \code{data.frame}
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
#' @param breaks optionally specify the breaks manually using a
#' vector giving the break points between cells, i.e., the vector
#' has length \code{num.bins - 1}.
#'
#' @return
#' A data.frame containing sample covariance values, bin centers,
#' and number \code{n} of data point pairs whose distance falls in
#' the corresponding bin. The first
#' row of this data.frame corresponds to the sample variance.
#' If \code{breaks} are specified manually then locations of bin
#' centers are not returned.
#'
#' @examples
#'
#' @export
covCMB <- function(cmbdf,
                   num.bins = 10,
                   sample.size,
                   max.dist,
                   breaks)
{

  if ( !is.CMBDataFrame(cmbdf) ) {

    stop("cmbdf must be a CMBDataFrame")

  }

  if ( !missing(sample.size) ) {

    cmbdf <- sampleCMB(cmbdf, sample.size = sample.size)

  }

  if ( is.null(coords(cmbdf)) || coords(cmbdf) != "cartesian" ) {

    coords(cmbdf) <- "cartesian"

  }

  if ( !all.equal(names(cmbdf), c("x", "y", "z", "I")) )
  {
    stop("cmbdf must have columns named 'x', 'y', 'z', 'I' in that order")
  }

  if ( missing(breaks) )
  {

    if (missing(max.dist)) {

      max.dist <- maxDist_internal(cmbdf)

    }

    marks <- seq(0, max.dist, length.out = num.bins+1)[-1]
    bin.width <- marks[1]
    centers <- c(0, marks - bin.width/2)
    breaks <- marks[-num.bins]

    covs <- covCMB_internal1(cmbdf, breaks)
    result <- data.frame(dist = centers, cov = covs[,1], n = as.integer(c(covs[,2][1],covs[,2][-1]/2)) )
  }
  else
  {
    if ( !missing(max.dist) || num.bins != 10 )
    {
      warning(paste0("The max.dist and num.bins parameters have no effect",
                     " when breaks is specified"))
    }

    covs <- covCMB_internal1(cmbdf, breaks)
    result <- data.frame(cov = covs[,1], n = as.integer(c(covs[,2][1],covs[,2][-1]/2)) )
  }

  return(result)
}



