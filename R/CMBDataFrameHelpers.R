#' Area of a HEALPix pixel
#'
#' Get the area of a single HEALPix pixel
#'
#'@param cmbdf a \code{\link{CMBDataFrame}}
#'
#'@return the area of a single HEALPix pixel
#' at the \code{nside} resolution of \code{cmbdf}
#'
#' @examples
#'
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#' pixelArea(df)
#'
#' df1 <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' pixelArea(df1)
#'
#'@export
pixelArea <- function(cmbdf)
{
  nside <- nside(cmbdf)
  if ( !is.numeric(nside) ) stop("problem with cmbdf nside parameter")
  return(pi/(3*nside^2))
}










#' Get the FITS headers from a \code{\link{CMBDataFrame}}
#'
#'
#'
#'@param cmbdf a CMBDataFrame.
#'
#'@return
#' The FITS headers belonging to the FITS file from which cmbdf
#' data was imported
#'
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#' df.sample <- CMBDataFrame(df, sample.size = 10000)
#' header(df.sample)
#'
#'@export
header <- function( cmbdf )
{
  # Check that argument is a CMBDF
  if ( !rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  return( c(attr( cmbdf, "header1" ),attr( cmbdf, "header2" )) )
}




#' Get the arcmin resolution from a \code{\link{CMBDataFrame}}
#'
#'@param cmbdf a CMBDataFrame.
#'
#'@return
#' The arcmin resolution as specified by the FITS file where the
#' data was sourced
#'
#' @examples
#'
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#' resolution(df)
#'
#'@export
resolution <- function( cmbdf )
{
  # Check that argument is a CMBDF
  if ( !rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  return( attr( cmbdf, "resolution" ) )
}





#' Take a simple random sample from a CMBDataFrame
#'
#' This function returns a CMBDataFrame with size sample.size,
#' whose rows comprise a simple random sample of the rows
#' from the input CMBDataFrame.
#'
#'@param cmbdf a CMB Data Frame.
#'@param sample.size the desired sample size.
#'
#'@return
#' A CMBDataFrame with size sample.size,
#' whose rows comprise a simple random sample of the rows
#' from the input CMBDataFrame.
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#' plot(sampleCMB(df, sample.size = 800000))
#'
#'@export
sampleCMB <- function(cmbdf, sample.size)
{
  if ( !rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  srows <- sample(1:nrow(cmbdf), sample.size)
  cmbdf[srows, ]
}




#' First Minkowski functional
#'
#' This function returns an area of the spherical region where measured values
#' are above of the specified threshold level \eqn{alpha}.
#'
#'@param cmbdf a CMB Data Frame.
#'@param \eqn{alpha} a threshold level
#'@param varindex an index of CMBDataFrame column with measured values
#'@return
#' The area of the exceedance region
#'
#' @references  Leonenko N., Olenko A. (2014) Sojourn measures of Student
#' and Fisher-Snedecor random fields.  Bernoulli, 20:1454-1483.
#'
#'@examples
#'
#' n <- 64
#' cmbdf <- CMBDataFrame(nside=n, I = rnorm(12*n^2), coords = "cartesian", ordering = "nested")
#' fmf(cmbdf, 0, 4)
#' fmf(cmbdf, 2, 4)
#'
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' cmbdf.win <- window(cmbdf, new.window = win)
#' fmf(cmbdf.win, 0, 4)
#'
#'@export
fmf <- function(cmbdf, alpha, varindex)
{
  if ( !rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  pixelArea(cmbdf)*sum(cmbdf[,varindex] > alpha)
}

