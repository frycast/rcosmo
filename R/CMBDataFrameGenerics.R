#' Get a sub window from \code{\link{CMBDataFrame}}
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
#'If the argument \code{new.window} is a list of CMBWindows,
#'then interious of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned CMBDataFrame will
#'be the points of \code{cmbdf} in the the intersection of \eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned CMBDataFrame
#'consists of the points of \code{x} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned CMBDataFrame
#'will be the points of \code{x} in \eqn{B} (resp. \eqn{A}).
#'
#'@param x A \code{\link{CMBDataFrame}}.
#'@param new.window Optionally specify a new window
#'in which case a new CMBDataFrame is returned whose CMBWindow is new.window.
#'\code{new.window} may also be a list (see details section and examples).
#'@param intersect A boolean that determines
#'the behaviour when \code{new.window} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param in.pixels A vector of pixels at resolution
#'\code{in.pixels.res} whose union contains the
#'window(s) \code{new.window} entirely, or if \code{new.window} is
#'unspecified then this whole pixel is returned.
#'@param in.pixels.res An integer. Resolution
#'(i.e., \eqn{j} such that nside = \code{2^j}) at
#'which the \code{in.pixels} parameter is specified
#'@param ... Unused arguments.
#'
#'@return
#' The window attribute of \code{x} or, if new.window/in.pixels is specified,
#' a new CMBDataFrame.
#'
#'@examples
#'
#'
#' ## Example 1: Create a new CMBDataFrame with a window
#'
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian",
#'                       ordering = "nested")
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' cmbdf.win <- window(cmbdf, new.window = win)
#' plot(cmbdf.win)
#' window(cmbdf.win)
#'
#' ## Example 2: Change the window of an existing CMBDataFrame
#'
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' window(cmbdf) <- win2 <- CMBWindow(theta = c(pi/6,pi/3,pi/3, pi/6),
#'                                    phi = c(0,0,pi/6,pi/6))
#' plot(cmbdf)
#'
#' ## Example 3: union of windows
#'
#' ## Create 2 windows
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(2*pi/3,3*pi/4,3*pi/4, 2*pi/3),
#'                             phi = c(pi/4,pi/4,pi/3,pi/3))
#' plot(win1)
#' plot(win2)
#'
#'## Create CMBDataFrame with points in the union of win1 and win2
#'
#'cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#'cmbdf.win <- window(cmbdf, new.window = list(win1, win2), intersect = FALSE)
#'plot(cmbdf.win)
#'
#'#' ## Example 4: intersection of windows
#'
#' ## Create 2 windows
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4),
#'                   phi = c(pi/4,pi/4,pi/3,pi/3))
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
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2),
#'                   set.minus =TRUE)
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4),
#'                   phi = c(pi/4,pi/4,pi/3,pi/3),
#'                   set.minus =TRUE)
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
#' win2 <- CMBWindow(theta = c(pi/4,pi/3,pi/3, pi/4),
#'                   phi = c(pi/4,pi/4,pi/3,pi/3),
#'                   set.minus =TRUE)
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
window.CMBDataFrame <- function(x, new.window, intersect = TRUE,
                                in.pixels, in.pixels.res = 0, ...)
{
  cmbdf <- x

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

      return(subWindow(cmbdf = cmbdf, win = new.window,
                                intersect = intersect,
                                in.pixels = in.pixels,
                                in.pixels.res = in.pixels.res))
    }
    else
    {
      return(subWindow(cmbdf = cmbdf, win = new.window,
                                intersect = intersect))
    }
  }
  # We reach this point only if new.window is missing
  if ( !missing(in.pixels) )
  {
    if (max(in.pixels) > 12*(2^in.pixels.res)^2)
    {
      stop("in.pixels out of range specified by in.pixels.res")
    }

    if (rcosmo::ordering(cmbdf) != "nested")
    {
      stop("in.pixel can only be used with nested ordering")
    }

    pixelWin <- rcosmo::pixelWindow(in.pixels.res,
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
`window<-.CMBDataFrame` <- function(x, ..., value)
{
  return(rcosmo::window(x, new.window = value))
}




#' HEALPix ordering scheme from a CMBDataFrame
#'
#' This function returns the HEALPix ordering scheme from a CMBDataFrame.
#' The ordering scheme is either "ring" or "nested".
#'
#' If a new ordering is specified, using e.g. new.ordering = "ring", the
#' ordering scheme of the CMBDataFrame will be converted.
#'
#'@param x A \code{\link{CMBDataFrame}}.
#'@param new.ordering Specifies the new ordering ("ring" or "nest")
#'if a change of ordering
#'scheme is desired.
#'@param ... Unused arguments.
#'
#'@return
#' The name of the HEALPix ordering scheme that is
#' used in the CMBDataFrame x.
#'
#'@examples
#' df <- CMBDataFrame(nside = 1, ordering = "nested")
#' ordering(df)
#' df1 <- ordering(df, new.ordering = "ring")
#' ordering(df1)
#'
#'@export
ordering.CMBDataFrame <- function( x, new.ordering, ... )
{
  cmbdf <- x

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
      pix(cmbdf) <- rcosmo::ring2nest(nside = nside(cmbdf),
                                       pix = pix(cmbdf))
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



#' Assign a new ordering scheme to a \code{\link{CMBDataFrame}}
#'
#' @keywords internal
#'
#' @seealso \code{\link{ordering.CMBDataFrame}}
#'
#' @examples
#' cmbdf <- CMBDataFrame(n = 1, ordering = "ring")
#' ordering(cmbdf) <- "nested"
#' ordering(cmbdf)
#'
#' @export
`ordering<-.CMBDataFrame` <- function(x,...,value) {
  rcosmo::ordering(x, new.ordering = value)
}





#' HEALPix Nside parameter from a CMBDataFrame
#'
#' This function returns the HEALPix Nside parameter of a CMBDataFrame
#'
#'@param x A \code{\link{CMBDataFrame}}.
#'
#'@return
#' The HEALPix Nside parameter.
#'
#'@examples
#' df <- CMBDataFrame(nside = 16)
#' nside(df)
#'
#'@export
nside.CMBDataFrame <- function( x )
{
  return( as.integer(attr( x, "nside" )) )
}




#' HEALPix pixel indices from \code{\link{CMBDataFrame}}
#'
#' If new.pix is unspecified then this function returns the vector of
#' HEALPix pixel indices from a CMBDataFrame. If new.pix is specified then
#' this function returns a new CMBDataFrame with the same number of rows
#' as \code{cmbdf}, but with pix attribute \code{new.pix}. Thus,
#' \code{new.pix} must have length equal to \code{nrow(cmbdf)}.
#'
#'@param x A CMBDataFrame.
#'@param new.pix Optional vector of pixel indices with
#'length equal to \code{nrow(x)}.
#'@param ... Unused arguments.
#'
#'@return
#' The vector of HEALPix pixel indices or, if new.pix is specified,
#' a new CMBDataFrame.
#'
#'@examples
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' # pix(df)
#'
#' df <- CMBDataFrame(nside = 16, sample.size = 10, ordering = "nested")
#' pix(df)
#'
#'@export
pix.CMBDataFrame <- function(x, new.pix, ...)
{

  if ( !missing(new.pix) )
  {
    row.names(x) <- new.pix
    return(x)
  }

  return(as.integer( row.names(x) ))
}



#' Assign new pixel indices to a CMBDataFrame
#'
#' @keywords internal
#'
#' @export
`pix<-.CMBDataFrame` <- function(x,...,value) {
  row.names(x) <- value
  x
}



# This prevents cmbdf[,i] from dropping attributes of cmbdf
#'@export
`[.CMBDataFrame` <- function(x,i,j, ..., drop) {

  mdrop <- missing(drop)

  # nargs counts empty arguments so a[1,] has nargs() = 3
  Narg <- nargs() - (!mdrop)

  if ( Narg > 2 & mdrop ) drop <- FALSE
  if ( Narg > 2 ) r <- NextMethod("[", drop = drop)
  if ( Narg <= 2 ) r <- NextMethod("[")

  attr(r, "ordering") <- attr(x, "ordering")
  attr(r, "nside") <- attr(x, "nside")
  attr(r, "coords") <- attr(x, "coords")
  attr(r, "window") <- attr(x, "window")
  attr(r, "header1") <- attr(x, "header1")
  attr(r, "header2") <- attr(x, "header2")
  attr(r, "resolution") <- attr(x, "resolution")
  r
}





#'\code{\link{cbind}} for CMBDataFrames
#'
#'Add a new column or columns (vector, matrix or data.frame)
#'to a \code{\link{CMBDataFrame}}. Note that method dispatch
#'occurs on the first argument. So, the CMBDataFrame must
#'be the first argument
#'
#'See the documentation for \code{\link{cbind}}
#'
#'@param ... (generalized) vectors or matrices. Columns to bind.
#'@param deparse.level Integer controlling the construction of labels
#'in the case of non-matrix-like arguments.
#'
#'
#'@examples
#'cmbdf <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
#'cmbdf2 <- cbind(cmbdf, myData = rep(1, 12))
#'cmbdf2
#'
#'
#'@export
cbind.CMBDataFrame <- function(..., deparse.level = 1)
{
  args <- list(...)

  is.cmbdf <- sapply(args, is.CMBDataFrame)
  if ( sum(is.cmbdf) != 1 )
  {
    stop("Just 1 CMBDataFrame must be passed to '...'")
  }

  # Keep this CMBDataFrame so we can get its attributes
  cmbdf <- args[is.cmbdf][[1]]
  args[is.cmbdf][[1]] <- as.data.frame(args[is.cmbdf][[1]])

  df <- do.call(cbind, c(args, deparse.level = deparse.level))

  class(df) <- class(cmbdf)
  nam <- names(df)
  attr(df, "ordering") <- attr(cmbdf, "ordering")
  attr(df, "nside") <- attr(cmbdf, "nside")
  attr(df, "coords") <- attr(cmbdf, "coords")
  attr(df, "window") <- attr(cmbdf, "window")
  attr(df, "row.names") <- attr(cmbdf, "row.names")
  attr(df, "header1") <- attr(cmbdf, "header1")
  attr(df, "header2") <- attr(cmbdf, "header2")
  attr(df, "resolution") <- attr(cmbdf, "resolution")
  names(df) <- nam

  return(df)
}






#'\code{\link{rbind}} for CMBDataFrames
#'
#'Add a new row or rows to a \code{\link{CMBDataFrame}}.
#'All arguments passed to \code{...} must be CMBDataFrames.
#'If the CMBDataFrame arguments have overlapping pixel
#'indices then all but one of the non-unique rows will be
#'deleted unless \code{unsafe = TRUE}. If \code{unsafe = TRUE}
#'then a \code{\link{HPDataFrame}} will be returned instead
#'of a \code{\link{CMBDataFrame}}.
#'
#'@param ... A number of CMBDataFrames
#'@param unsafe A boolean. If the CMBDataFrame arguments have
#'overlapping pixel
#'indices then all but one of the non-unique rows will be
#'deleted unless \code{unsafe = TRUE}. If \code{unsafe = TRUE}
#'then a \code{\link{HPDataFrame}} will be returned instead
#'of a \code{\link{CMBDataFrame}}.
#'@param deparse.level See documentation
#'for \code{\link{rbind.data.frame}}.
#'
#'@seealso See the documentation for \code{\link{rbind}}
#'
#'@examples
#' df <- CMBDataFrame(nside = 1, I = 1:12)
#'
#' df.123 <- CMBDataFrame(df, spix = c(1,2,3))
#' df.123
#' df.234 <- CMBDataFrame(df, spix = c(2,3,4))
#' df.234
#'
#' df.1234 <- rbind(df.123, df.234)
#' df.1234
#' class(df.1234) # A CMBDataFrame
#' pix(df.1234)
#'
#' df.123234 <- rbind(df.123, df.234, unsafe = TRUE)
#' df.123234
#' class(df.123234) # A HPDataFrame
#' pix(df.123234)
#'
#'@export
rbind.CMBDataFrame <- function(..., deparse.level = 1, unsafe = FALSE)
{
  args <- list(...)

  if ( !all(sapply(args, rcosmo::is.CMBDataFrame)) )
  {
    stop("rbind.CMBDataFrame requires all arguments to be CMBDataFrames")
  }

  if ( !all(sapply(args, areCompatibleCMBDFs, cmbdf2 = args[[1]])) )
  {
    stop(paste0("The CMBDataFrames are not compatible for rbind.\n",
                "You may need to convert coordinates or ordering"))
  }

  if ( !unsafe )
  {
    # Make sure that no pixels are repeated, using pairwise intersections
    for (i in 1:(length(args)-1))
    {
      for (j in (i+1):length(args))
      {
        cap.ij <- intersect(pix(args[[i]]), pix(args[[j]]))

        if ( length(cap.ij) != 0 )
        {
          to.drop <- which(pix(args[[i]]) %in% cap.ij)
          args[[i]] <- args[[i]][-to.drop,]
        }
      }
    }
  }

  # Keep this CMBDataFrame and the windows so we can get attributes
  cmbdf <- args[[1]]
  wins <- sapply(args, rcosmo::window)

  if ( unsafe )
  {
    # Get pix before we convert to data.frames
    pix <- unlist(lapply(args, pix))
  }

  args <- lapply(args, as.data.frame)
  df <- do.call(rbind, c(args, deparse.level = deparse.level))

  if ( !unsafe )
  {
    class(df) <- class(cmbdf)
    attr(df, "nside") <- nside(cmbdf)
    attr(df, "ordering") <- ordering(cmbdf)
    attr(df, "coords") <- coords(cmbdf)
    attr(df, "window") <- wins
    attr(df, "header1") <- attr(cmbdf, "header1")
    attr(df, "header2") <- attr(cmbdf, "header2")
    attr(df, "resolution") <- attr(cmbdf, "resolution")
  }
  else
  {
    len <- nrow(df)
    if (len != length(pix)) stop("(development stage) unexpected error")

    row.names(df) <- 1:len
    attr(df, "pix") <- pix
    attr(df, "nside") <- nside(cmbdf)
    attr(df, "ordering") <- ordering(cmbdf)
    class(df) <- unique(c("HPDataFrame", class(df)))
  }

  return(df)
}



#Reduce(intersect, list(pix(a.w1),pix(a.w2),pix(a.w3)))

#' Check compatibleness of CMBDataFrames
#'
#' Compare attributes to decide if two CMBDataFrames are compatible
#'
#' If the CMBDataFrames do not have compatible attributes then
#' a message is printed indicating the attributes that do not match.
#' To suppress this use the \code{\link{suppressMessages}} function
#'
#' @param cmbdf1 a \code{\link{CMBDataFrame}}
#' @param cmbdf2 a \code{\link{CMBDataFrame}}
#' @param compare.pix A boolean. If TRUE then cmbdf1 and
#' cmbdf2 must share the same pixel indices to be considered
#' compatible
#'
#' @examples
#' a <- CMBDataFrame(nside = 2, ordering = "ring", coords = "cartesian")
#' b <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
#' areCompatibleCMBDFs(a,b)
#'
#' suppressMessages(areCompatibleCMBDFs(a,b))
#'
#' @export
areCompatibleCMBDFs <- function(cmbdf1, cmbdf2, compare.pix = FALSE)
{
  if ( !is.CMBDataFrame(cmbdf1) || !is.CMBDataFrame(cmbdf1) )
  {
    stop("cmbdf1 and cmbdf2 must be CMBDataFrames")
  }

  ns <- identical(nside(cmbdf1), nside(cmbdf2))
  ord <- identical(ordering(cmbdf1), ordering(cmbdf2))
  pix <- identical(pix(cmbdf1), pix(cmbdf2))

  reasons <- ""
  if (!ns)
  {
    reasons <- paste0(reasons, "nside mismatch (nside1 = ", nside(cmbdf1),
                     ", nside2 = ", nside(cmbdf2), ")\n")
  }
  if (!ord)
  {
    reasons <- paste0(reasons, "ordering mismatch (ordering1 = ",
                      ordering(cmbdf1),
                      ", ordering2 = ", ordering(cmbdf2), ")\n")
  }
  if (!pix && compare.pix)
  {
    reasons <- paste0(reasons, "pixels mismatch pix(cmbdf1) != pix(cmbdf2)\n")
  }

  if ( !(ns && ord && (pix || !compare.pix) ) )
  {
    message(reasons)
    return(FALSE)
  }
  else
  {
    return(TRUE)
  }
}







#' Convert dataframes to CMBDataFrames
#'
#'
#' Safely converts a \code{\link{data.frame}} to a CMBDataFrame. The
#' rows of the data.frame are assumed to be in the HEALPix order
#' given by \code{ordering}, and at the HEALPix resolution given
#' by \code{nside}. Coordinates, if present,  are assumed to correspond to
#' HEALPix pixel centers. The coordinates must be named either x,y,z
#' (cartesian) or theta, phi (spherical colatitude and longitude respectively).
#'
#' @param df Any \code{data.frame} whose rows are in HEALPix order
#' @param ordering character string that specifies the ordering scheme
#' ("ring" or "nested")
#' @param nside an integer \eqn{2^k} that specifies the Nside (resolution)
#' HEALPix parameter
#' @param spix an integer vector that specifies the HEALPix pixel index
#' corresponding to each row of \code{df}. If \code{spix} is left blank and
#' \code{df} is a \code{data.frame}, then \code{df} is assumed to contain data
#' for every pixel at resolution parameter \code{nside} (the full sky).
#' In other words,
#' in this case, the number of rows of \code{df} must be equal to 12*nside^2.
#' However, if \code{spix} is left blank and \code{df}
#' is a \code{CMBDataFrame},
#' then \code{spix} is set equal to \code{pix(df)}
#'
#' @return A CMBDataFrame
#'
#' @examples
#'
#' ## Example 1: Create df with no coords, then create CMBDataFrames cmbdf and
#' ## df2 with spherical coords
#'
#' df <- data.frame(I=rnorm(12))
#' df
#'
#' cmbdf <- as.CMBDataFrame(df,ordering= "ring", nside=1)
#' summary(cmbdf)
#' pix(cmbdf)
#' coords(cmbdf)
#'
#' df2 <- coords(cmbdf, new.coords = "spherical")
#' df2
#'
#' ## Example 2: Create CMBDataFrames for first 10 Healpix centers
#'
#' df <- data.frame(I=rnorm(10))
#' df
#' cmbdf <- as.CMBDataFrame(df,ordering= "ring", nside=2, spix=1:10)
#' summary(cmbdf)
#' pix(cmbdf)
#'
#' @export
as.CMBDataFrame <- function(df, ordering, nside, spix)
{
  if ( !is.data.frame(df) ) {

    stop(gettextf("'%s' is not a data.frame", deparse(substitute(df))))

  }

  if ( !is.CMBDataFrame(df) ) {
    ################ df IS NOT A CMBDataFrame ####################

    if(missing(nside) || missing(ordering))
    {
      stop(paste0("nside and ordering must both be specified if df is",
                  " not already a CMBDataFrame"))
    }

    if(missing(spix))
    {
      spix <- 1:(12*nside^2)
    }

    if ( nrow(df) != length(spix) )
    {
      stop(paste0("The number of rows of df must equal the length of spix, ",
                  "or if spix is unspecified then df must have number of ",
                  "rows equal to 12*nside^2"))
    }

    attr(df, "ordering") <- ordering
    attr(df, "nside") <- nside
    attr(df, "row.names") <- spix


    if ( ("theta" %in% names(df) && "phi" %in% names(df)) ) {
      attr(df, "coords") <- "spherical"
    } else if ( ("x" %in% names(df) && "y" %in% names(df)
          && "z" %in% names(df) ) ) {
      attr(df, "coords") <- "cartesian"
    } else {
      attr(df, "coords") <- NULL
    }

  } else {
    ################ df IS  A CMBDataFrame #######################

    ## We do some checks and change nothing

    if(!missing(spix) && spix != pix(df))
    {
      stop(paste0("df is a CMBDataFrame so spix must match",
                  " pix(df) or be unspecified"))
    }

    if (!missing(ordering) && rcosmo::ordering(df) != tolower(ordering) ) {
      stop(gettextf("Since '%s' is a CMBDataFrame already,
                    ordering should be unspecified or match the ordering
                    attribute of '%s'",
                    deparse(substitute(df)),
                    deparse(substitute(df))))
    }

    if (!missing(nside) && rcosmo::nside(df) != tolower(nside) ) {
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

















#' Check if an object is of class CMBDataFrame
#'
#' @param cmbdf Any R object
#'
#' @return TRUE if \code{cmbdf} is a CMBDataFrame, otherwise FALSE
#'
#'@examples
#'
#' df <- CMBDataFrame(nside = 16)
#' is.CMBDataFrame(df)
#' df2 <- coords(df, new.coords = "cartesian")
#' is.CMBDataFrame(df2)
#'
#'
#' @export
is.CMBDataFrame <- function(cmbdf)
{
  identical(as.numeric(sum(class(cmbdf) == "CMBDataFrame")), 1)
}











#' Geodesic area covered by a \code{\link{CMBDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{cmbdf}
#'
#'@param x a CMBDataFrame.
#'
#'@return The sum of the areas of all pixels (rows) in x.
#'
#'@examples
#'
#' ## At low resolution, a few data points can
#' ## occupy a large pixel area, e.g.:
#'
#' cmbdf <- CMBDataFrame(nside = 1, spix = c(1,2,3))
#' pix(cmbdf)
#'
#' ## The total number of Healpix points at nside=1 equals 12. As cmbdf has 3 Helpix
#' ## it occupies pi = 1/4*(surface area of unit sphere)
#'
#' geoArea(cmbdf)
#' plot(cmbdf, size = 5, hp.boundaries = 1)
#'
#'@export
geoArea.CMBDataFrame <- function(x)
{
  nside <- rcosmo::nside(x)
  if ( !is.numeric(nside) ) stop("problem with cmbdf nside attribute")
  return(pi/(3*nside^2)*nrow(x))
}




#' Coordinate system from a \code{\link{CMBDataFrame}}
#'
#' If \code{new.coords} is unspecified then
#' this function returns the coordinate system used in the CMBDataFrame
#' \code{cmbdf}.
#' The coordinate system is either "cartesian" or "spherical".
#' If a new coordinate system is specified, using e.g.
#' \code{new.coords = "spherical"}, then this function instead
#' returns a new CMBDataFrame whose coordinates are of the specified
#' type. The original CMBDataFrame, \code{cmbdf}, is unaffected.
#' If you would like to change \code{cmbdf} without creating a new
#' variable, then use \code{\link{coords<-.CMBDataFrame}} (see
#' examples below).
#'
#'
#'@param x A CMBDataFrame, \code{cmbdf}.
#'@param new.coords Specifies the new coordinate
#'system ("spherical" or "cartesian")
#'if a change of coordinate system is desired.
#'@param ... Unused arguments.
#'
#'@return
#' If new.coords is unspecified, then the name of the coordinate system
#' of \code{cmbdf} is returned. Otherwise a new CMBDataFrame is returned
#' equivalent to \code{cmbdf} but having the desired change of coordinates
#'
#'@examples
#'
#' ## Create df with no coords, then create df2 with cartesian coords
#' df <- CMBDataFrame(nside = 16)
#' df
#' coords(df)
#' df2 <- coords(df, new.coords = "cartesian")
#' coords(df2)
#'
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' coords(df)
#'
#'
#'
#'@export
coords.CMBDataFrame <- function( x, new.coords, ... )
{
  cmbdf <- x

  # If new.coords argument is missing then return the coordinate type
  if ( missing(new.coords) )
  {
    return(attr(cmbdf, "coords"))
  }
  else
  {
    new.coords <- as.character(tolower(new.coords))

    if ( new.coords != "spherical" && new.coords != "cartesian" )
    {
      stop("new.coords must be 'spherical' or 'cartesian'")
    }

    if ( is.null(attr(cmbdf, "coords")) )
    {
      cart <- (new.coords == "cartesian")
      nest <- (ordering(cmbdf) == "nested")
      ns <- rcosmo::nside(cmbdf)
      sp <- rcosmo::pix(cmbdf)

      if (cart)
      {
        nc <- 3
        nam <- c("x","y","z")
      }
      else
      {
        nc <- 2
        nam <- c("theta","phi")
      }

      crds <- pix2coords_internal(nside = ns, nested = nest,
                                  cartesian = cart, spix = sp)[,1:nc]
      crds <- as.data.frame(crds)
      names(crds) <- nam

      cmbdf <- cbind.CMBDataFrame(crds, cmbdf)
    }
    else if ( attr(cmbdf, "coords") == new.coords )
    {
      # Nothing to do
    }
    else if ( new.coords == "spherical" )
    {
      # Convert to spherical

      x.i <- which(names(cmbdf) == "x")
      y.i <- which(names(cmbdf) == "y")
      z.i <- which(names(cmbdf) == "z")

      crds <- cmbdf[,c(x.i, y.i, z.i)]
      crds <- car2sph(crds)
      other <- cmbdf[,-c(x.i, y.i, z.i), drop = FALSE]
      cmbdf <- cbind.CMBDataFrame(crds, other)
      attr(cmbdf, "coords") <- "spherical"
    }
    else if ( new.coords == "cartesian" )
    {
      # convert to cartesian

      theta.i <- which(names(cmbdf) == "theta")
      phi.i <- which(names(cmbdf) == "phi")

      crds <- cmbdf[,c(theta.i, phi.i)]
      crds <- sph2car(crds)
      other <- cmbdf[,-c(theta.i, phi.i), drop = FALSE]
      cmbdf <- cbind.CMBDataFrame(crds, other)
      attr(cmbdf, "coords") <- "cartesian"
    }

    attr(cmbdf, "coords") <- new.coords
    return(cmbdf)
  }
}



#' Assign new coordinate system to a \code{\link{CMBDataFrame}}
#'
#'@keywords internal
#'
#' @seealso \code{\link{coords.CMBDataFrame}}
#'
#' @examples
#'
#' ## Create df with no coords, then create df2 with cartesian coords
#' df <- CMBDataFrame(nside = 16)
#' df
#' coords(df)
#' df2 <- coords(df, new.coords = "cartesian")
#' coords(df2)
#' coords(df)
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' df
#'
#' @export
`coords<-.CMBDataFrame` <- function(x, ..., value) {
  return(coords(x, new.coords = value))
}







#' Plot CMB Data
#'
#' This function produces a plot from a \code{\link{CMBDataFrame}}.
#'
#'@param x A \code{\link{CMBDataFrame}}.
#'@param intensities The name of a column that specifies CMB intensities.
#'This is only used if \code{col} is unspecified.
#'@param add If TRUE then this plot will be added to any existing plot.
#'Note that if \code{back.col} (see below) is specified then a new plot
#'window will be opened and \code{add = TRUE} will have no effect.
#'@param sample.size Optionally specifies the size of a simple random
#'sample to take before plotting. This can make the plot less
#'computationally intensive.
#'@param type A single character indicating the type of item to plot.
#'Supported types are: 'p' for points, 's' for spheres, 'l' for lines,
#''h' for line segments from z = 0, and 'n' for nothing.
#'@param size The size of plotted points.
#'@param box Whether to draw a box.
#'@param axes Whether to draw axes.
#'@param aspect Either a logical indicating whether to adjust the
#'aspect ratio, or a new ratio.
#'@param col Specify the colour(s) of the plotted points.
#'@param back.col Optionally specifies the background colour of
#'the plot. This argument is passed to rgl::bg3d.
#'@param labels Optionally specify a vector of labels to plot,
#'such as words or vertex indices. If this is specified then
#'\code{rgl::text3d} is used instead of \code{rgl::plot3d}. Then
#'\code{length(labels)} must equal \code{nrow(x)}.
#'@param hp.boundaries Integer. If greater than 0 then HEALPix
#'pixel boundaries at \code{nside = hp.boundaries} will be
#'added to the plot.
#'@param hpb.col Colour for the \code{hp.boundaries}.
#'@param ... Arguments passed to rgl::plot3d.
#'
#'@return
#'A plot of the CMB data
#'
#'@examples
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # sky <- CMBDataFrame("CMB_map_smica1024.fits")
#' # plot(sky, sample.size = 800000)
#'
#'@export
plot.CMBDataFrame <- function(x, intensities = "I",
                              add = FALSE, sample.size,
                              type = "p", size = 1, box = FALSE,
                              axes = FALSE, aspect = FALSE,
                              col, back.col = "black", labels,
                              hp.boundaries = 0, hpb.col = "gray",
                              ...) {

  cmbdf <- x
  if ( !missing(sample.size) ) {

    spix <- sample(rcosmo::pix(cmbdf), sample.size)
    cmbdf <- cmbdf[rcosmo::pix(cmbdf) %in% spix,]
  }

  if (is.null(coords(cmbdf))) {

    rcosmo::coords(cmbdf) <- "cartesian"
  }

  if (missing(col)) {

     col <- tryCatch(colscheme(cmbdf[,intensities, drop = TRUE],
                           breaks1024, colmap),
                     error = function(e) {
                      stop(paste0("Problem producing CMB colour scheme. ",
                                  "Make sure that the intensities ",
                                  "parameter matches a column name of cmbdf. ",
                                  "Or, try specifying colours manually ",
                                  "with the col parameter."),
                           call. = FALSE)
                     })
  }

  ## Change coordinates if necessary
  cmbdf.xyz <- rcosmo::coords(cmbdf, new.coords = "cartesian")

  ## Do the plotting
  if ( !add ) {

    rgl::open3d()
    rgl::bg3d(back.col)
  }

  if ( missing(labels) ) {

    rgl::plot3d(cmbdf.xyz$x, cmbdf.xyz$y, cmbdf.xyz$z,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  } else {

    rgl::text3d(cmbdf.xyz$x, cmbdf.xyz$y, cmbdf.xyz$z, labels,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }

  if ( hp.boundaries > 0 ) {

    rcosmo::displayPixelBoundaries(nside = hp.boundaries,
                                    col = hpb.col)
  }
}

# Helper function for plot.CMBDataFrame
colscheme <- function(I, breaks, colmap) {

  intervals<- findInterval(I, breaks[2:256], rightmost.closed = FALSE,
                           all.inside = FALSE, left.open = TRUE)
  cols <- colmap[intervals+1]
  return(cols)
}







#' Summarise a \code{\link{CMBDataFrame}}
#'
#' This function produces a summary from a CMBDataFrame.
#'
#'@param object A \code{\link{CMBDataFrame}}.
#'@param intensities the name of a column specifying CMB intensities
#'(or potentially another numeric quantity of interest)
#'@param ... Unused arguments.
#'
#'@return
#'A summary includes window's type and area,
#' total area covered by observations,
#' and main statistcs for intensity values
#'
#'
#'@examples
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # df <- CMBDataFrame("CMB_map_smica1024.fits")
#' # df.sample <- CMBDataFrame(df, sample.size = 800000)
#' # summary(df.sample)
#'
#' ns <- 16
#' df <- CMBDataFrame(I = rnorm(12*ns^2), nside = ns,
#'                    ordering = "nested")
#'
#' win1 <- CMBWindow(x=0,y=3/5,z=4/5,r=0.8)
#' df.sample1 <- window(df, new.window = win1)
#' summary(df)
#'
#'@export
summary.CMBDataFrame <- function(object, intensities = "I", ...)
{
  cmbdf <- object
  ans <- list(intensities = summary(cmbdf[,intensities, drop = TRUE]))

  if ( is.null(rcosmo::coords(cmbdf)) )
  {
    ans$coords <- "HEALPix only"
  }
  else
  {
    ans$coords <- rcosmo::coords(cmbdf)
  }

  if ( is.null(rcosmo::window(cmbdf)) )
  {
    ans$window <- "full sky"
  }
  else
  {
    ans$window <- rcosmo::window(cmbdf)
  }

  if ( is.null(resolution(cmbdf)) )
  {
    ans$resolution <- "unknown"
  }
  else
  {
    ans$resolution <- rcosmo::resolution(cmbdf)
  }

  ans$ordering <- rcosmo::ordering(cmbdf)
  ans$nside <- rcosmo::nside(cmbdf)
  ans$pix <- rcosmo::pix(cmbdf)
  ans$n <- nrow(cmbdf)
  ans$area <- rcosmo::geoArea(cmbdf)
  ans$method <- rcosmo::header(cmbdf)[grepl("METHOD  =",
                                            rcosmo::header(cmbdf))]

  class(ans) <- "summary.CMBDataFrame"
  return(ans)
}


#'Print a summary of a CMBDataFrame
#'
#'@keywords internal
#'
#'
#'@param x a \code{summary.CMBDataFrame} object, i.e.,
#'the output of \code{\link{summary.CMBDataFrame}}
#'
#'@export
#'
#'
print.summary.CMBDataFrame <- function(x, ...)
{
  cat(
    cli::rule(center = " CMBDataFrame Object ", line = "bar4"), "\n",
    sep = ""
  )

  # Window loop details here in boxes
  if ( any(x$window != "full sky") )
  {
    cat("Number of CMBWindows: ", length(x$window), "\n" )
    if ( length(x$window) <= 5 )
    {
      lapply(x$window, function(x) { print(summary(x)); cat("\n\n") } )
    }
    else
    {
      cat("Too many windows to print them all here", "\n\n")
    }
  }
  else
  {
    cat("Full sky map\n")
  }
  cat(x$method, "\n")

  cat("Total area covered by all pixels: ", x$area, "\n")



  # Summary of intensities
  cat(cli::rule(line = "~"), "\n", sep = "")
  cli::cat_line("Intensity quartiles")
  print(x$intensities)

  # Finishing line
  cat(
    cli::rule(line = "="),
    sep = ""
  )
}


#' Print CMBDataFrame
#'
#' This function neatly prints the contents of a CMBDataFrame.
#'
#'@param x A \code{\link{CMBDataFrame}}.
#'@param ... arguments passed to \code{\link{print.tbl_df}}
#'
#'@return
#'Prints contents of the CMB data frame to the console.
#'
#'
#'@export
print.CMBDataFrame <- function(x, ...)
{
  cat("A CMBDataFrame\n")
  print(tibble::as.tibble(x), ...)
}





