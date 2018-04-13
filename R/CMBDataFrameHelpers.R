
#'Restrict a \code{\link{CMBDataFrame}} to a \code{\link{CMBWindow}}
#'
#'A single CMBWindow or a list of CMBWindows can be passed to the \code{win}
#'argument
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows: Let \eqn{A} be the union of the
#'interiors of all windows whose winType does not include "minus",
#'and let \eqn{B} be the intersection of the exteriors of all the windows whose
#'\code{winType} does include "minus". Then, provided that
#'\code{intersect = TRUE} (the default), the returned CMBDataFrame will
#'be the intersection of the points in \code{cmbdf} with \eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned CMBDataFrame will
#'be the intersection of the points in \code{cmbdf} with the union of
#'\eqn{A} and \eqn{B}.
#'Note that if \eqn{A} (resp. \eqn{B}) is empty
#'then the returned CMBDataFrame will be the intersection of \eqn{B}
#'(resp. \eqn{A}) with \code{cmbdf}.
#'
#'@param cmbdf a CMBDataFrame
#'@param win a CMBWindow or a list of CMBWindows
#'@param intersect a boolean that determines
#'the behaviour when \code{win} is a list (see details).
#'
#'@return a CMBDataFrame which is restricted to the
#'region of the sky specified by \code{win}
#'
#'@examples
#'
#'@export
subWindow <- function(cmbdf, win, intersect = TRUE)
{
  if ( !rcosmo::is.CMBDataFrame(cmbdf) ) {
    stop("'cmbdf' must be a CMBDataFrame")
  }

  if ( !rcosmo::is.CMBWindow(win) ) {

    if ( !is.list(win) )
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }

    if (!all(sapply(win, rcosmo::is.CMBWindow)))
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }
  }
  else
  {
    win <- list(win)
  }

  # subsequent operations will require cartesian coordinates
  win.xyz <- lapply(win, rcosmo::coords, new.coords = "cartesian")
  cmbdf.xyz <- rcosmo::coords(cmbdf, new.coords = "cartesian")

  # Triangulate all non-convex polygons into lists of convex polygons,
  # making win.xyz into a new list of only convex polygons
  win.conv <- list()
  for ( w in win.xyz )
  {
    #Triangulate each w and add all triangles to win.conv
    if ( rcosmo:::contains("polygon", winType(w)) && !rcosmo::assumedConvex(w))
    {
        win.conv <- append(win.conv, rcosmo::triangulate(w))
    }
    else
    {
      win.conv[[length(win.conv)+1]] <- w
    }
  }
  win.xyz <- win.conv


  ## pointInside happens here
  exist.m <- FALSE
  exist.p <- FALSE
  keep.p <- rep(FALSE, nrow(cmbdf.xyz))
  keep.m <- rep(TRUE, nrow(cmbdf.xyz))
  for ( w in win.xyz )
  {
    switch(rcosmo::winType(w),
           polygon = keep.p <- keep.p |
             rcosmo::pointInConvexPolygon(cmbdf.xyz[,c("x","y","z")], w),
           minus.polygon = keep.m <- keep.m &
             !rcosmo::pointInConvexPolygon(cmbdf.xyz[,c("x","y","z")], w),
           disc = keep.p <- keep.p |
             rcosmo::pointInDisc(cmbdf.xyz[,c("x","y","z")], w),
           minus.disc = keep.m <- keep.m &
             !rcosmo::pointInDisc(cmbdf.xyz[,c("x","y","z")], w),
           stop("Failed to determine window type using rcosmo::winType"))

    if ( rcosmo:::contains("minus", rcosmo::winType(w)) )
    {
      exist.m <- TRUE
    }
    else
    {
      exist.p <- TRUE
    }

  }

  if ( intersect && exist.p )
  {
    keep <- keep.m & keep.p
  }
  else if ( exist.m && exist.p )
  {
    keep <- keep.m | keep.p
  }
  else if ( exist.p )
  {
    keep <- keep.p
  }
  else
  {
    keep <- keep.m
  }


  cmbdf.new <- cmbdf[keep, , drop = FALSE]

  attr(cmbdf.new, "window") <- win

  return(cmbdf.new)
}





#' Window attribute of \code{\link{CMBDataFrame}}
#'
#' When new.window is unspecified this function returns the
#' \code{\link{CMBWindow}} attribute of a
#' CMBDataFrame. The return value is NULL if the window is full sky.
#' When new.window is specified this function instead returns
#' a new CMBDataFrame whose CMBWindow attribute is new.window
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows: Let \eqn{A} be the union of the
#'interiors of all windows in the list, \code{new.window}, whose winType
#'does not include "minus",
#'and let \eqn{B} be the intersection of the exteriors of all the windows whose
#'\code{winType} does include "minus". Then, provided that
#'\code{intersect = TRUE} (the default), the returned CMBDataFrame will
#'be the intersection of the points in \code{cmbdf} with \eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned CMBDataFrame will
#'be the intersection of the points in \code{cmbdf} with the union of
#'\eqn{A} and \eqn{B}.
#'Note that if \eqn{A} (resp. \eqn{B}) is empty
#'then the returned CMBDataFrame will be the intersection of \eqn{B}
#'(resp. \eqn{A}) with \code{cmbdf}.
#'
#'@param cmbdf a CMBDataFrame.
#'@param new.window optionally specify a new window
#'in which case a new CMBDataFrame is returned whose CMBWindow is new.window.
#'\code{new.window} may also be a list (see details section).
#'@param intersect a boolean that determines
#'the behaviour when \code{win} is a list (see details).
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
window <- function(cmbdf, new.window, intersect = TRUE)
{
  if ( !missing(new.window) )
  {
    return(rcosmo::subWindow(cmbdf = cmbdf, win = new.window,
                              intersect = intersect))
  }

  return(attr(cmbdf, "window"))
}




#' Assign a new \code{\link{CMBWindow}} to a \code{\link{CMBDataFrame}}
#'@export
`window<-` <- function(cmbdf,...,value)
{
  return(rcosmo::window(cmbdf, new.window = value))
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
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  if ( !missing(new.pix) )
  {
    row.names(cmbdf) <- new.pix
    return(cmbdf)
  }

  return(as.integer( row.names(cmbdf) ))
}





#' Assign new pixel indices to a CMBDataFrame
#' @export
`pix<-` <- function(cmbdf,...,value) {
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
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
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  if ( missing(new.ordering) )
  {
    return(attr( cmbdf, "ordering" ))

  } else {
    new.ordering <- as.character(tolower(new.ordering))

    if ( identical(as.character(attr(cmbdf, "ordering")), new.ordering) )
    {
      # Nothing to do

    } else if ( identical(new.ordering, "nested") ) {

      message("Converting to nested ordering...\n")
      pix(cmbdf) <- ring2nest(nside = nside(cmbdf), pix = pix(cmbdf))
      cmbdf <- cmbdf[order(pix(cmbdf)),]
      attr(cmbdf, "ordering") <- "nested"

    } else if ( identical(new.ordering, "ring") ) {

      message("Converting to ring ordering...\n")
      pix(cmbdf) <- nest2ring(nside = nside(cmbdf), pix = pix(cmbdf))
      cmbdf <- cmbdf[order(pix(cmbdf)),]
      attr(cmbdf, "ordering") <- "ring"

    } else {

      stop("new.ordering must be either 'ring' or 'nested'")

    }

    return(cmbdf)
  }
}

#' Assign new ordering scheme to CMBDataFrame
#' @export
`ordering<-` <- function(cmbdf,...,value) {
  rcosmo::ordering(cmbdf, new.ordering = value)
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
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }

  return( as.integer(attr( cmbdf, "nside" )) )
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
  if ( !rcosmo::is.CMBDataFrame(cmbdf) )
  {
    stop("Argument must be a CMBDataFrame")
  }
  spix <- sample(pix(cmbdf), sample.size)
  cmbdf[spix, ]
}




#' First Minkowski functional
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

