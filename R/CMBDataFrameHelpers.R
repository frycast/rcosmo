#' Area of a HEALPix pixel
#'
#' Get the area of a single HEALPix pixel
#'
#'@param nsideObject  \code{\link{CMBDataFrame}}, a
#'\code{\link{HPDataFrame}}, or an integer
#'giving the nside parameter.
#'
#'@return the area of a single HEALPix pixel
#' at the \code{nside} resolution of \code{nsideObject}
#'
#' @examples
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # pixelArea(df)
#'
#' df1 <- CMBDataFrame(nside = 64,
#'                     coords = "cartesian",
#'                     ordering = "nested")
#' pixelArea(df1)
#'
#'@export
pixelArea <- function(nsideObject)
{
  if ( is.CMBDataFrame(nsideObject) || is.HPDataFrame(nsideObject) ) {
    nside <- nside(nsideObject)

  } else if ( is.numeric(nsideObject) ) {

    nside <- as.integer(nsideObject)
  } else {

    stop("nsideObject must be a CMBDataFrame, HPDataFrame, or integer")
  }

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
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # df.sample <- CMBDataFrame(df, sample.size = 10000)
#' # header(df.sample)
#'
#'@export
header <- function( cmbdf )
{
  # Check that argument is a CMBDF
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
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
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # resolution(df)
#'
#'@export
resolution <- function( cmbdf )
{
  # Check that argument is a CMBDF
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  return( attr( cmbdf, "resolution" ) )
}




