
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
#'a \code{\link{HPDataFrame}}
#'or \code{\link{CMBDat}} object. If this is a data.frame then it must have
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
#'(i.e., \eqn{j} such that nside = \code{2^j}) at
#'which the \code{in.pixels} parameter is specified
#'
#'@keywords internal
#'
#'@return a CMBDataFrame, HPDataFrame, or just a data.frame,
#'which is restricted to the
#'region of the sky specified by \code{win}
#'
#'@export
subWindow <- function(cmbdf, win, intersect = TRUE, in.pixels,
                      in.pixels.res = 0)
{
  if ( !rcosmo::is.CMBDataFrame(cmbdf)
       && !rcosmo::is.CMBDat(cmbdf) )
  {
    if ( is.data.frame(cmbdf) )
    {
      # Note that HPDataFrames end up in here too
      cmbdf <- rcosmo::coords(cmbdf, new.coords = "cartesian")

      if ( is.HPDataFrame(cmbdf) )
      {
        spx <- rcosmo::pix(cmbdf)
      }
    }
    else
    {
      stop("cmbdf must be a data.frame, HPDataFrame, CMBDataFrame or CMBDat object")
    }
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

  if ( !missing(in.pixels) && rcosmo::is.CMBDataFrame(cmbdf) )
  {
    if ( rcosmo::ordering(cmbdf) != "nested" )
    {
      stop("in.pixel can only be used with nested ordering")
    }
    pixelWin <- rcosmo::pixelWindow(in.pixels.res,
                                     log2(nside(cmbdf)),
                                     in.pixels)
    cmbdf <- cmbdf[pixelWin,]
  }

  # subsequent operations will require win to have cartesian coordinates
  win.xyz <- lapply(win, rcosmo::coords, new.coords = "cartesian")

  # Triangulate all non-convex polygons into lists of convex polygons,
  # making win.xyz into a new list of only convex polygons
  win.conv <- list()
  for ( w in win.xyz )
  {
    #Triangulate each w and add all triangles to win.conv
    if ( contains("polygon", winType(w)) && !rcosmo::assumedConvex(w))
    {
      win.conv <- append(win.conv, triangulate(w))
    }
    else
    {
      win.conv[[length(win.conv)+1]] <- w
    }
  }
  win.xyz <- win.conv


  exist.m <- FALSE
  exist.p <- FALSE
  if ( rcosmo::is.CMBDataFrame(cmbdf)
       || is.data.frame(cmbdf)
       || is.HPDataFrame(cmbdf) )
  {
    keep.p <- rep(FALSE, nrow(cmbdf))
    keep.m <- rep(TRUE, nrow(cmbdf))
  }
  else
  {
    if ( !rcosmo::is.CMBDat(cmbdf) )
    {
      stop("cmbdf must be a CMBDataFrame, data.frame or CMBDat object")
    }
    keep.p <- rep(FALSE, 12*cmbdf$nside^2)
    keep.m <- rep(TRUE, 12*cmbdf$nside^2)
  }
  ## pointInside happens here
  if ( (rcosmo::is.CMBDataFrame(cmbdf)
        && (!is.null(coords(cmbdf)) && coords(cmbdf) == "cartesian"))
       || (!rcosmo::is.CMBDataFrame(cmbdf) && is.data.frame(cmbdf)))
  { ## This might be a little quicker if coords are already cartesian
    ## Note that HPDataFrames end up in here
    #################################################################
    if ( !all( c("x","y","z") %in% names(cmbdf) ) )
    {
      stop(paste0("If cmbdf is a data.frame it must have columns named x,y,z ",
                  "that specify Cartesian coordinates"))
    }

    for ( w in win.xyz )
    {
      switch(rcosmo::winType(w),
             polygon = keep.p <- keep.p |
               pointInConvexPolygon(cmbdf[,c("x","y","z")], w),
             minus.polygon = keep.m <- keep.m &
               !pointInConvexPolygon(cmbdf[,c("x","y","z")], w),
             disc = keep.p <- keep.p |
               pointInDisc(cmbdf[,c("x","y","z")], w),
             minus.disc = keep.m <- keep.m &
               !pointInDisc(cmbdf[,c("x","y","z")], w),
             stop("Failed to determine window type using rcosmo::winType"))

      if ( contains("minus", rcosmo::winType(w)) )
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
    if ( rcosmo::is.CMBDataFrame(cmbdf) || rcosmo::is.HPDataFrame(cmbdf) )
    {
      nes <- (rcosmo::ordering(cmbdf) == "nested")
      ns <- nside(cmbdf)
      spx <- rcosmo::pix(cmbdf)
    }
    else if ( rcosmo::is.CMBDat(cmbdf) )
    {
      nes <- (cmbdf$ordering == "nested")
      ns <- cmbdf$nside
      spx <- 1:(12*ns^2)
    }


    for ( w in win.xyz )
    {
      switch(rcosmo::winType(w),
             polygon = keep.p <- keep.p |
               pointInConvexPolygonHP(nside = ns,
                                               nested = nes,
                                               win = w, spix = spx),
             minus.polygon = keep.m <- keep.m &
               !pointInConvexPolygonHP(nside = ns,
                                                nested = nes,
                                                win = w, spix = spx),
             disc = keep.p <- keep.p |
               pointInDiscHP(nside = ns,
                                      nested = nes,
                                      win = w, spix = spx),
             minus.disc = keep.m <- keep.m &
               !pointInDiscHP(nside = ns,
                                       nested = nes,
                                       win = w, spix = spx),
             stop("Failed to determine window type using rcosmo::winType"))

      if ( contains("minus", rcosmo::winType(w)) )
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

  # Note that HPdataFrame and CMBDataFrame are both data.frames
  if ( is.data.frame(cmbdf) )
  {
    cmbdf.new <- cmbdf[keep, , drop = FALSE]

    # We just need to set the pix attribute if it is HPDataFrame
    if ( is.HPDataFrame(cmbdf) )
    {
      attr(cmbdf.new, "pix") <- spx[keep]
      row.names(cmbdf.new) <- 1:nrow(cmbdf.new)
    }
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





#' Get a sub window from a data.frame
#'
#' This function returns a
#' data.frame containing the data in \code{x} restricted to the
#' CMBWindow \code{new.window}
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows.
#'
#'If the argument is a list of CMBWindows, then interiors of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned data.frame will
#'be the points of \code{x} in the the intersection of
#'\eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned data.frame
#'consists of the points of \code{x} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned data.frame
#'will be the points of \code{x} in \eqn{B} (resp. \eqn{A}).
#'
#'@param x A data.frame. Must have
#'columns labelled x,y,z specifying cartesian coordinates, or
#'columns labelled theta, phi specifying colatitude and longitude
#'respectively.
#'@param new.window A single \code{\link{CMBWindow}} object or a list of them.
#'@param intersect A boolean that determines
#'the behaviour when \code{new.window} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param ... Unused arguments.
#'
#'@return
#' A data.frame containing the data in \code{x} restricted to the
#' CMBWindow \code{new.window}
#'
#'@examples
#'
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#'
#' cmbdf <- CMBDataFrame(nside = 4)
#' df2 <- coords(cmbdf, new.coords = "cartesian")
#' df <- as.data.frame(df2[,1:3])
#' df
#' df.win <- window(df, new.window = win1)
#' df.win
#'
#'@export
window.data.frame <- function(x, new.window, intersect = TRUE, ...)
{
  return(subWindow(x, win = new.window, intersect = intersect))
}
