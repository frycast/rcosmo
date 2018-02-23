







#'Restrict a \code{\link{CMBDataFrame}} to a \code{\link{CMBWindow}}
#'
#'@param cmbdf a CMBDataFrame
#'@param win a CMBWindow
#'
#'@return a CMBDataFrame which is restricted to the
#'region of the sky specified by \code{win}
#'
#'@examples
#'
#'@export
subWindow <- function(cmbdf, win)
{
  if ( !is.CMBDataFrame(cmbdf) ) {
    stop("'cmbdf' must be a CMBDataFrame")
  }

  if ( !is.CMBWindow(win) ) {

    if ( !is.list(win) )
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }

    if (!all(sapply(new.window, is.CMBWindow)))
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }
  }
  else
  {
    win <- list(win)
  }


  # subsequent operations will require cartesian coordinates
  win.xyz <- lapply(win, coords, new.coords = "cartesian")
  cmbdf.xyz <- coords(cmbdf, new.coords = "cartesian")



  # Triangulate all non-convex polygons into lists of convex polygons,
  # then make win.xyz into a new list of only convex polygons
  win.conv <- list()
  for ( w in win.xyz )
  {
    #Triangulate each w and add all triangles to win.conv
    if ( winType(w) == "polygon" || winType(w) == "minus.polygon" )
    {
      win.conv <- append(win.conv, ifelse(assumedConvex(w), w, triangulate(w)))
    }
    else
    {
      win.conv <- append(win.conv, w)
    }
  }



  # All the work is done by pointInConvexPolygon: src/CMBDataFrameHelpers.cpp

  keep <- rep(FALSE, nrow(cmbdf))
  for ( w in win.xyz )
  {
    # union all windows
    keep <- keep |
      switch(winType(w),
             polygon = pointInConvexPolygon(cmbdf[,c("x","y","z")], w),
             minus.polygon = !pointInConvexPolygon(cmbdf[,c("x","y","z")],w),
             disc = stop(paste("(development stage) subWindow not",
                               "implemented for discs in window list")),
             minus.disc = stop(paste("(development stage) subWindow not",
                                     "implemented for discs in window list")),
             stop("Failed to determine window type using rcosmo::winType"))
  }

  cmbdf.new <- cmbdf[keep,]

  attr(cmbdf.new, "window") <- win

  return(cmbdf.new)
}

## HELPER FUNCTION FOR subWindow
triangulate <- function(win)
{
  win.conv <- list()

  # LOOP --------- action happens here:
    tri <- win # dummy line

    # all triangles are convex and share the winType of win
    attr(tri, "assumedConvex") <- TRUE
    attr(tri, "winType") <- winType(win)
    win.conv <- append(win.conv, tri)
  # -----------------------------------

  return(win.conv)
}






#' Window attribute of \code{\link{CMBDataFrame}}
#'
#' When new.window is unspecified this function returns the
#' \code{\link{CMBWindow}} attribute of a
#' CMBDataFrame. The return value is NULL if the window is full sky.
#' When new.window is specified this function instead returns
#' a new CMBDataFrame whose CMBWindow attribute is new.window
#'
#'@param cmbdf a CMBDataFrame.
#'@param new.window optionally specify a new window, in which case a
#'new CMBDataFrame is returned whose CMBWindow is new.window
#'
#'@return
#' The window attribute of cmbdf or, if new.window is specified, a
#' new CMBDataFrame.
#'
#'@examples
#' cmbdf <- CMBDataFrame(nside = 16, coords = "cartesian", ordering = "nested")
#'
#' ## Create a new CMBDataFrame with a window
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' cmbdf.win <- window(cmbdf, new.window = win)
#' plot(cmbdf.win)
#' window(cmbdf.win)
#'
#' ## Change the window of an existing CMBDataFrame
#' window(cmbdf) <- CMBWindow(theta = rep(0.1, 10),
#'                            phi = seq(0, 9*2*pi/10, length.out = 10))
#' plot(cmbdf)
#'
#'@export
window <- function(cmbdf, new.window)
{
  if ( !missing(new.window) )
  {
    return(subWindow(cmbdf, new.window))
  }

  return(attr(cmbdf, "window"))
}




#' Assign a new \code{\link{CMBWindow}} to a \code{\link{CMBDataFrame}}
#'@export
`window<-` <- function(cmbdf,...,value)
{
  return(subWindow(cmbdf, value))
}






#' HEALPix pixel indices from \code{\link{CMBDataFrame}}
#'
#' If new.pix is unspecified then this function returns the vector of
#' HEALPix pixel indices from a CMBDataFrame. If new.pix is specified then
#' this function returns a new CMBDataFrame with pixel indices new.pix
#'
#'@param cmbdf a CMBDataFrame.
#'@param new.pix optional vector of pixel indices
#'
#'@return
#' The vector of HEALPix pixel indices or, if new.pix is specified,
#' a new CMBDataFrame.
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' pix(df)
#'
#'@export
pix <- function(cmbdf, new.pix)
{
  # Check that argument is a CMBDF
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  if ( !missing(new.pix) )
  {
    row.names(cmbdf) <- new.pix
    return(cmbdf)
  }

  return(as.numeric( row.names(cmbdf) ))
}





#' Assign new pixel indices to a CMBDataFrame
#' @export
`pix<-` <- function(cmbdf,...,value) {
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  row.names(cmbdf) <- value
  cmbdf
}





#' HEALPix ordering scheme from a CMBDataFrame
#'
#' This function returns the HEALPix ordering scheme from a CMBDataFrame.
#' The ordering scheme is either "ring" or "nested".
#'
#' If a new ordering is specified, using e.g. new.ordering = "ring", the
#' ordering scheme of the CMBDataFrame will be converted.
#'
#'@param cmbdf a CMB Data Frame.
#'@param new.ordering specifies the new ordering ("ring" or "nest") if a change of ordering
#'scheme is desired.
#'
#'@return
#' The name of the HEALPix ordering scheme that is used in the CMBDataFrame cmbdf
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' ordering(df)
#' ordering(df, new.ordering = "ring")
#'
#'@export
ordering <- function( cmbdf, new.ordering )
{
  # Check that argument is a CMBDF
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  if ( missing(new.ordering) )
  {
    return(attr( cmbdf, "ordering" ))

  } else {
    new.ordering <- as.character(tolower(new.ordering))

    if( !identical(new.ordering, "ring") & !identical(new.ordering,"nested") )
    {
      stop("new.ordering must be either 'ring' or 'nested'")
    }

    if ( identical(as.character(attr(cmbdf, "ordering")), new.ordering) )
    {
      # Nothing to do

    } else if ( identical(new.ordering, "nested") )
    {
      cat("Converting to nested ordering...\n")
      # Convert 'ring' to 'nested' ordering
      cat("Conversion not completed as this function is under development\n")

    } else if ( identical(new.ordering, "ring") )
    {
      cat("Converting to ring ordering...\n")
      # Convert 'nested' to 'ring' ordering
      cat("Conversion not completed as this function is under development\n")
    }
  }
}

#' Assign new ordering scheme to CMBDataFrame
#' @export
`ordering<-` <- function(cmbdf,...,value) {
  ordering(cmbdf, new.ordering = value)
  cmbdf
}







#' HEALPix Nside parameter from a CMBDataFrame
#'
#' This function returns the HEALPix Nside parameter of a CMBDataFrame
#'
#'@param cmbdf a CMB Data Frame.
#'
#'@return
#' The HEALPix Nside parameter
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' nside(df)
#'
#'@export
nside <- function( cmbdf )
{
  # Check that argument is a CMBDF
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  attr( cmbdf, "nside" )
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
  if ( !is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  spix <- sample(pix(cmbdf), sample.size)
  cmbdf[spix,]
}




#' Check if an object is of class CMBDataFrame
#'
#' @param cmbdf Any R object
#'
#' @return TRUE if \code{cmbdf} is a CMBDataFrame, otherwise FALSE
#'
#' @export
is.CMBDataFrame <- function(cmbdf)
{
  identical(as.numeric(sum(class(cmbdf) == "CMBDataFrame")), 1)
}





#' as.CMBDataFrame
#'
#' Safely converts a data.frame to a CMBDataFrame
#'
#' @param df Any data.frame with a column labelled "I" for intensities
#' @param coords specifies the coordinate system to be "spherical",
#' "cartesian" or unspecified (HEALPix only). If "spherical" then df
#' must have columns named "theta" and "phi" (colatitude and longitude
#' respectively). If "cartesian" then df
#' must have columns named "x", "y", and "z"
#' @param ordering specifies the ordering scheme ("ring" or "nested")
#' @param nside specifies the Nside parameter
#'
#' @return A CMBDataFrame
#'
#' @export
as.CMBDataFrame <- function(df, coords, ordering, nside)
{
  if ( !is.data.frame(df) ) {

    stop(gettextf("'%s' is not a data.frame", deparse(substitute(df))))

  }

  if ( !("I"  %in% names(df) ) ) {

    stop(gettextf("'%s' does not have a column named 'I' for intensities",
                  deparse(substitute(df))))

  }


  if ( !is.CMBDataFrame(df) ) {
  ################ df IS NOT A CMBDataFrame ####################

    attr(df, "ordering") <- ordering
    attr(df, "nside") <- nside

    if ( missing(coords) ) {

      if ( any(c("theta", "phi", "x", "y", "z") %in% names(df) ) ) {
        warning("coords was unspecified and so coordinates
                were set to HEALPix only")
      }
      attr(df, "coords") <- NULL

    } else {

      coords <- tolower(coords)

      if ( coords == "spherical"
        && !("theta" %in% names(df)
        &&   "phi" %in% names(df)) ) {
        stop(gettextf("Since coords = spherical, '%s' must have
                      column names %s and %s",
                      deparse(substitute(df)),
                      dQuote("theta"),
                      dQuote("phi")))
      }

      if ( coords == "cartesian"
          && !("x" %in% names(df)
          &&   "y" %in% names(df)
          &&   "z" %in% names(df) ) ) {
        stop(gettextf("Since coords = cartesian, '%s' must have
                      column names %s, %s and %s",
                      deparse(substitute(df)),
                      dQuote("x"),
                      dQuote("y"),
                      dQuote("z")))
      }

      if (coords != "spherical" && coords != "cartesian") {
        stop("coords must be unspecified, 'spherical' or 'cartesian'")
      }
      attr(df, "coords") <- coords

    }

  } else {
  ################ df IS  A CMBDataFrame #######################

    coords(df) <- ifelse(!missing(coords), coords, coords(df))
    ordering(df) <- ifelse(!missing(ordering), ordering, ordering(df))

    if (!missing(nside) && nside(df) != tolower(nside) ) {
      stop(gettextf("Since '%s' is a CMBDataFrame already,
                    nside should be unspecified or match the nside
                    attribute of '%s'",
                    deparse(substitute(df)),
                    deparse(substitute(df))))
    }

  }

  class(df) <- c("CMBDataFrame", "data.frame")
  return(df)
}
