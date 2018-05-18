# Compute the empirical covariance between observations in a CMBDataFrame
# assuming homogeneity and isotropy.
# INPUTS: CMBDataFrame, vector r of radii
# OUTPUTS: C(r), for each r


#### WARNING: THE LAST BIN IS NOT THE RIGHT SIZE AS IT CONTAINS ALL
#### DISTANCES GREATER THAN max.dist SO IT SHOULD PERHAPS BE DISCARDED
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
#' maximum geodesic distance to use for calculating covariance. Only
#' used if \code{breaks} is unspecified.
#' @param breaks optionally specify the breaks manually using a
#' vector giving the break points between cells. This vector
#' has length \code{num.bins} since the last break point is taken
#' as \code{max.dist}. If \code{equiareal = TRUE} then
#' these breaks should be \eqn{cos(r_i)} where \eqn{r_i} are radii.
#' If \code{equiareal = FALSE} then these breaks should be \eqn{r_i}.
#' @param equiareal if TRUE then the bins have equal spherical
#' area. If false then the bins have equal annular widths.
#' Default is TRUE.
#' @param calc.max.dist if TRUE then the \code{max.dist} will be
#' calculated from the locations in cmbdf. Otherwise
#' either \code{max.dist} must be specified or max.dist will
#' default to pi.
#'
#' @return
#' An object of class CMBCovariance consisting of a \code{data.frame}
#' containing sample covariance
#' values, bin centers, and number \code{n} of data point pairs
#' whose distance falls in
#' the corresponding bin.
#' The first
#' row of this data.frame corresponds to the sample variance.
#' The attribute "breaks" contains the break points used.
#' The returned \code{data.frame} has
#' \code{num.bins + 1} rows since the first row, the sample
#' variance, is not counted as a bin.
#'
#' @examples
#'
#' @export
covCMB <- function(cmbdf,
                   num.bins = 10,
                   sample.size,
                   max.dist = pi,
                   breaks,
                   equiareal = TRUE,
                   calc.max.dist = FALSE)
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


    if (calc.max.dist)
    {
      max.dist <- maxDist_internal(cmbdf)
    }


    if (equiareal)
    {
      breaks <- seq(cos(max.dist), 1, length.out = num.bins+1)[-(num.bins+1)]
    }
    else
    {
      breaks <- cos(rev(seq(0, max.dist, length.out = num.bins+1)[-1]))
    }

  }

  covs <- covCMB_internal2(cmbdf, breaks)

  # Reverse order since cosine is decreasing
  covs[2:nrow(covs),1] <- rev(covs[2:nrow(covs),1])
  covs[2:nrow(covs),2] <- rev(covs[2:nrow(covs),2])
  v <- c(0,rev(acos(breaks)))
  # Drop the throw-away bin (distances greater than max.dist)
  covs <- covs[-nrow(covs),]

  # Find the bin centers (break_{i+1} - break_i)/2 with break_0 = 0.
  centers <- c(0,v[1:(length(v)-1)] + (v[2:length(v)] - v[1:(length(v) - 1)])/2)

  result <- data.frame(dist = centers, cov = covs[,1], n = as.integer(c(covs[,2][1], covs[,2][-1]/2)) )
  class(result) <- c("CMBCovariance", "data.frame")
  attr(result, "breaks") <- breaks

  return(result)
}



