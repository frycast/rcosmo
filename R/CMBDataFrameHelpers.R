#' pixelArea
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









#'subWindow
#'
#'Restricts a \code{\link{CMBDataFrame}}, \code{CMBDat} object,
#'or \code{\link{data.frame}} to a \code{\link{CMBWindow}} region.
#'A single CMBWindow or a list of CMBWindows can be passed to the \code{win}
#'argument.
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows.
#'
#'If the argument is a list of CMBWindows, then interious of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned CMBDataFrame will
#'be the points of \code{cmbdf} in the the intersection of \eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned CMBDataFrame
#'consists of the points of \code{cmbdf} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned CMBDataFrame
#'will be the points of \code{cmbdf} in \eqn{B} (resp. \eqn{A}).
#'
#'
#'@param cmbdf a \code{\link{CMBDataFrame}}, a \code{data.frame},
#'or CMBDat object. If this is a data.frame then it must have
#'columns labelled x,y,z specifying cartesian coordinates, or
#'columns labelled theta, phi specifying colatitude and longitude
#'respectively.
#'@param win a \code{\link{CMBWindow}} or a list of CMBWindows
#'@param intersect a boolean that determines
#'the behaviour when \code{win} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param in.pixels a vector of pixels at resolution
#'\code{in.pixels.res} whose union contains the
#'window(s) \code{win} entirely. This will only be used
#'if \code{cmbdf} is a \code{CMBDataFrame}
#'@param in.pixels.res a resolution
#'(i.e., \eqn{j} such that nside = \code{2^j|}) at
#'which the \code{in.pixels} parameter is specified
#'
#'@return a CMBDataFrame, or just a data.frame,
#'which is restricted to the
#'region of the sky specified by \code{win}
#'
#'@examples
#'
#'@export
subWindow <- function(cmbdf, win, intersect = TRUE, in.pixels,
                      in.pixels.res = 0)
{
  if ( !is.CMBDataFrame(cmbdf) && !is.CMBDat(cmbdf) )
  {
    if ( is.data.frame(cmbdf) )
    {
      cmbdf <- coords(cmbdf, new.coords = "cartesian")
    }
    else
    {
      stop("cmbdf must be a data.frame, CMBDataFrame or CMBDat object")
    }
  }

  if ( !rcosmo:::is.CMBWindow(win) ) {

    if ( !is.list(win) )
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }

    if (!all(sapply(win, rcosmo:::is.CMBWindow)))
    {
      stop("'win' must be a CMBWindow or list of CMBWindows")
    }
  }
  else
  {
    win <- list(win)
  }

  if ( !missing(in.pixels) && rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    if ( rcosmo:::ordering(cmbdf) != "nested" )
    {
      stop("in.pixel can only be used with nested ordering")
    }
    pixelWin <- rcosmo:::pixelWindow(in.pixels.res,
                                    log2(nside(cmbdf)),
                                    in.pixels)
    cmbdf <- cmbdf[pixelWin,]
  }

  # subsequent operations will require win to have cartesian coordinates
  win.xyz <- lapply(win, rcosmo:::coords, new.coords = "cartesian")

  # Triangulate all non-convex polygons into lists of convex polygons,
  # making win.xyz into a new list of only convex polygons
  win.conv <- list()
  for ( w in win.xyz )
  {
    #Triangulate each w and add all triangles to win.conv
    if ( rcosmo:::contains("polygon", winType(w)) && !rcosmo:::assumedConvex(w))
    {
        win.conv <- append(win.conv, rcosmo:::triangulate(w))
    }
    else
    {
      win.conv[[length(win.conv)+1]] <- w
    }
  }
  win.xyz <- win.conv


  exist.m <- FALSE
  exist.p <- FALSE
  if ( rcosmo:::is.CMBDataFrame(cmbdf) || is.data.frame(cmbdf) )
  {
    keep.p <- rep(FALSE, nrow(cmbdf))
    keep.m <- rep(TRUE, nrow(cmbdf))
  }
  else
  {
    if ( !rcosmo:::is.CMBDat(cmbdf) )
    {
      stop("cmbdf must be a CMBDataFrame, data.frame or CMBDat object")
    }
    keep.p <- rep(FALSE, 12*cmbdf$nside^2)
    keep.m <- rep(TRUE, 12*cmbdf$nside^2)
  }
  ## pointInside happens here
  if ( (rcosmo:::is.CMBDataFrame(cmbdf)
       && (!is.null(coords(cmbdf)) && coords(cmbdf) == "cartesian"))
       || (!rcosmo:::is.CMBDataFrame(cmbdf) && is.data.frame(cmbdf)))
  { ## This might be a little quicker if coords are already cartesian
    #################################################################
    if ( !all( c("x","y","z") %in% names(cmbdf) ) )
    {
     stop(paste0("If cmbdf is a data.frame it must have columns named x,y,z ",
                 "that specify Cartesian coordinates"))
    }

    for ( w in win.xyz )
    {
      switch(rcosmo:::winType(w),
             polygon = keep.p <- keep.p |
               rcosmo:::pointInConvexPolygon(cmbdf[,c("x","y","z")], w),
             minus.polygon = keep.m <- keep.m &
               !rcosmo:::pointInConvexPolygon(cmbdf[,c("x","y","z")], w),
             disc = keep.p <- keep.p |
               rcosmo:::pointInDisc(cmbdf[,c("x","y","z")], w),
             minus.disc = keep.m <- keep.m &
               !rcosmo:::pointInDisc(cmbdf[,c("x","y","z")], w),
             stop("Failed to determine window type using rcosmo::winType"))

      if ( rcosmo:::contains("minus", rcosmo:::winType(w)) )
      {
        exist.m <- TRUE
      }
      else
      {
        exist.p <- TRUE
      }
    }
  }
  else
  { ## This is quicker and more robust than assigning new coords in R
    #################################################################
    if ( rcosmo:::is.CMBDataFrame(cmbdf) )
    {
      nes <- (rcosmo:::ordering(cmbdf) == "nested")
      ns <- nside(cmbdf)
      spx <- rcosmo:::pix(cmbdf)
    }
    else if ( rcosmo:::is.CMBDat(cmbdf) )
    {
      nes <- (cmbdf$ordering == "nested")
      ns <- cmbdf$nside
      spx <- 1:(12*ns^2)
    }


    for ( w in win.xyz )
    {
      switch(rcosmo:::winType(w),
             polygon = keep.p <- keep.p |
               rcosmo:::pointInConvexPolygonHP(nside = ns,
                                              nested = nes,
                                              win = w, spix = spx),
             minus.polygon = keep.m <- keep.m &
               !rcosmo:::pointInConvexPolygonHP(nside = ns,
                                               nested = nes,
                                               win = w, spix = spx),
             disc = keep.p <- keep.p |
               rcosmo:::pointInDiscHP(nside = ns,
                                     nested = nes,
                                     win = w, spix = spx),
             minus.disc = keep.m <- keep.m &
               !rcosmo:::pointInDiscHP(nside = ns,
                                      nested = nes,
                                      win = w, spix = spx),
             stop("Failed to determine window type using rcosmo::winType"))

      if ( rcosmo:::contains("minus", rcosmo:::winType(w)) )
      {
        exist.m <- TRUE
      }
      else
      {
        exist.p <- TRUE
      }
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

  if ( is.CMBDataFrame(cmbdf) || is.data.frame(cmbdf) )
  {
    cmbdf.new <- cmbdf[keep, , drop = FALSE]
  }
  else if ( is.CMBDat(cmbdf) )
  {
    keep <- which(keep)
    cmbdf.new <- cmbdf$data[keep,]

    attr(cmbdf.new, "row.names") <- keep
    attr(cmbdf.new, "nside") <- cmbdf$nside
    class(cmbdf.new) <- c("CMBDataFrame","data.frame")
    attr(cmbdf.new, "ordering") <- cmbdf$ordering
    attr(cmbdf.new, "coords") <- NULL
    attr(cmbdf.new, "resolution") <- cmbdf$resoln
    attr(cmbdf.new, "header1") <- cmbdf$header1
    attr(cmbdf.new, "header2") <- cmbdf$header2
  }

  attr(cmbdf.new, "window") <- win

  return(cmbdf.new)
}








#' Window attribute of \code{\link{CMBDataFrame}}
#'
#' When new.window or in.pixels is unspecified this function returns the
#' \code{\link{CMBWindow}} attribute of a
#' CMBDataFrame. The return value is NULL if the window is full sky.
#' When new.window is specified this function instead returns
#' a new CMBDataFrame whose CMBWindow attribute is new.window
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows.
#'
#'If the argument is a list of CMBWindows, then interious of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned CMBDataFrame will
#'be the points of \code{cmbdf} in the the intersection of \eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned CMBDataFrame
#'consists of the points of \code{cmbdf} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned CMBDataFrame
#'will be the points of \code{cmbdf} in \eqn{B} (resp. \eqn{A}).
#'
#'@param cmbdf a CMBDataFrame.
#'@param new.window optionally specify a new window
#'in which case a new CMBDataFrame is returned whose CMBWindow is new.window.
#'\code{new.window} may also be a list (see details section and examples).
#'@param intersect a boolean that determines
#'the behaviour when \code{win} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param in.pixels a vector of pixels at resolution
#'\code{in.pixels.res} whose union contains the
#'window(s) \code{win} entirely, or if \code{new.window} is
#'unspecified then this whole pixel is returned
#'@param in.pixels.res a resolution
#'(i.e., \eqn{j} such that nside = \code{2^j|}) at
#'which the \code{in.pixels} parameter is specified
#'
#'@return
#' The window attribute of cmbdf or, if new.window/in.pixels is specified,
#' a new CMBDataFrame.
#'
#'@examples
#'
#'
#' ## Example 1: Create a new CMBDataFrame with a window
#'
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' cmbdf.win <- window(cmbdf, new.window = win)
#' plot(cmbdf.win)
#' window(cmbdf.win)
#'
#' ## Example 2: Change the window of an existing CMBDataFrame
#'
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' window(cmbdf) <- win2 <- CMBWindow(theta = c(pi/6,pi/3,pi/3, pi/6), phi = c(0,0,pi/6,pi/6))
#' plot(cmbdf)
#'
#' ## Example 3: union of windows
#'
#' ## Create 2 windows
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(2*pi/3,3*pi/4,3*pi/4, 2*pi/3), phi = c(pi/4,pi/4,pi/3,pi/3))
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the union of win1 and win2
#'
#'cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#'cmbdf.win <- window(cmbdf, new.window = list(win1, win2), intersect = TRUE)
#'plot(cmbdf.win)
#'
#'#' ## Example 4: intersection of windows
#'
#' ## Create 2 windows
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4), phi = c(pi/4,pi/4,pi/3,pi/3))
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the intersection of win1 and win2
#'
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' cmbdf.win1 <- window(cmbdf, new.window = win1)
#' cmbdf.win12 <- window(cmbdf.win1, new.window = win2)
#' plot(cmbdf.win12)
#' plot(win1)
#' plot(win2)
#'
#'
#'## Example 5: intersection of windows with "minus" type
#'
#' ## Create 2 windows with "minus" type
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2), set.minus =TRUE)
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4), phi = c(pi/4,pi/4,pi/3,pi/3), set.minus =TRUE)
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the intersection of win1 and win2
#'
#'cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#'cmbdf.win <- window(cmbdf, new.window = list(win1, win2))
#'plot(cmbdf.win)
#'
#'
#'## Example 6: intersection of windows with different types
#'
#' ##Create 2 windows, one with "minus" type
#'
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4), phi = c(pi/4,pi/4,pi/3,pi/3), set.minus =TRUE)
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the intersection of win1 and win2
#'
#'cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#'cmbdf.win <- window(cmbdf, new.window = list(win1, win2), intersect = TRUE)
#'plot(cmbdf.win)
#'
#' ## Example 7: union of windows with different types
#'
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2), set.minus =TRUE)
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4), phi = c(pi/4,pi/4,pi/3,pi/3))
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the union of win1 and win2
#'
#'cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#'cmbdf.win <- window(cmbdf, new.window = list(win1, win2), intersect = FALSE)
#'plot(cmbdf.win)
#'
#'
#'@export
window <- function(cmbdf, new.window, intersect = TRUE,
                   in.pixels, in.pixels.res = 0)
{
  if ( !missing(new.window) )
  {
    if ( !missing(in.pixels) )
    {
      if (max(in.pixels) > 12*(2^in.pixels.res)^2)
      {
        stop("in.pixels out of range specified by in.pixels.res")
      }

      if (ordering(cmbdf) != "nested")
      {
        stop("in.pixel can only be used with nested ordering")
      }

      return(rcosmo:::subWindow(cmbdf = cmbdf, win = new.window,
                               intersect = intersect,
                               in.pixels = in.pixels,
                               in.pixels.res = in.pixels.res))
    }
    else
    {
      return(rcosmo:::subWindow(cmbdf = cmbdf, win = new.window,
                               intersect = intersect))
    }
  }
  # We reach this point only if new.window is missing
  if ( !rcosmo:::is.CMBDataFrame(cmbdf) )
  {
    stop("If new.window is unspecified then cmbdf must be a CMBDataFrame")
  }
  if ( !missing(in.pixels) )
  {
    if (max(in.pixels) > 12*(2^in.pixels.res)^2)
    {
      stop("in.pixels out of range specified by in.pixels.res")
    }

    if (rcosmo:::ordering(cmbdf) != "nested")
    {
      stop("in.pixel can only be used with nested ordering")
    }

    pixelWin <- rcosmo:::pixelWindow(in.pixels.res,
                                    log2(nside(cmbdf)),
                                    in.pixels)

    return(cmbdf[pixelWin,])
  }

  return(attr(cmbdf, "window"))
}




#' Assign a new \code{\link{CMBWindow}} to a \code{\link{CMBDataFrame}}
#'
#' @keywords internal
#'
#'
#'@export
#'
`window<-` <- function(cmbdf,...,value)
{
  return(rcosmo:::window(cmbdf, new.window = value))
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

