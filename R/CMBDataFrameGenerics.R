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
#' df <- CMBDataFrame(nside = 1, ordering = "nested")
#' ordering(df)
#' ordering(df, new.ordering = "ring")
#'
#'@export
ordering.CMBDataFrame <- function( cmbdf, new.ordering )
{

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
      pix(cmbdf) <- rcosmo:::ring2nest(nside = nside(cmbdf),
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

#' Assign new ordering scheme to CMBDataFrame
#' @export
`ordering<-.CMBDataFrame` <- function(cmbdf,...,value) {
  rcosmo:::ordering(cmbdf, new.ordering = value)
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
#' df <- CMBDataFrame(nside = 16)
#' nside(df)
#'
#'@export
nside.CMBDataFrame <- function( cmbdf )
{
  return( as.integer(attr( cmbdf, "nside" )) )
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
pix.CMBDataFrame <- function(cmbdf, new.pix)
{

  if ( !missing(new.pix) )
  {
    row.names(cmbdf) <- new.pix
    return(cmbdf)
  }

  return(as.integer( row.names(cmbdf) ))
}



#' Assign new pixel indices to a CMBDataFrame
#' @export
`pix<-.CMBDataFrame` <- function(cmbdf,...,value) {
  row.names(cmbdf) <- value
  cmbdf
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






#'Like \code{\link{rbind}} for CMBDataFrames
#'
#'Add a new row or rows to a \code{\link{CMBDataFrame}}.
#'All arguments passed to \code{...} must be CMBDataFrames.
#'
#'@param unsafe defaults to FALSE. If \code{unsafe = TRUE} then
#'overlapping pixel coordinates will not throw an error (faster).
#'
#'See the documentation for \code{\link{rbind}}
#'
#'@export
rbind.CMBDataFrame <- function(..., deparse.level = 1, unsafe = FALSE)
{
  args <- list(...)

  if ( !all(sapply(args, rcosmo:::is.CMBDataFrame)) )
  {
    stop("rbind.CMBDataFrame requires all arguments to be CMBDataFrames")
  }

  if ( !all(sapply(args, rcosmo:::areCompatibleCMBDFs, cmbdf2 = args[[1]])) )
  {
    stop(paste0("The CMBDataFrames are not compatible for rbind.\n",
                "You may need to convert coordinates or ordering"))
  }

  if ( unsafe == FALSE )
  {
    # Make sure that no pixels are repeated, using pairwise intersections
    for (i in 1:(length(args)-1))
    {
      for (j in (i+1):length(args))
      {
        len <- length(intersect(pix(args[[i]]), pix(args[[j]])))

        if ( len != 0 )
        {
          stop("The CMBDataFrames passed to rbind overlap somewhere")
        }
      }
    }
  }

  # Keep this CMBDataFrame and the windows so we can get attributes
  cmbdf <- args[[1]]
  wins <- sapply(args, rcosmo:::window)

  args <- lapply(args, as.data.frame)
  df <- do.call(rbind, c(args, deparse.level = deparse.level))


  class(df) <- class(cmbdf)
  attr(df, "nside") <- nside(cmbdf)
  attr(df, "ordering") <- ordering(cmbdf)
  attr(df, "coords") <- coords(cmbdf)
  attr(df, "window") <- wins
  attr(df, "header1") <- attr(cmbdf, "header1")
  attr(df, "header2") <- attr(cmbdf, "header2")
  attr(df, "resolution") <- attr(cmbdf, "resolution")

  return(df)
}



#Reduce(intersect, list(pix(a.w1),pix(a.w2),pix(a.w3)))

#' areCompatibleCMBDFs
#'
#' Compare attributes to decide if two CMBDataFrames are compatible
#'
#' If the CMBDataFrames do not have compatible attributes then
#' a message is printed indicating the attributes that do not match.
#' To suppress this use the \code{\link{suppressMessages}} function
#'
#' @param cmbdf1 a \code{\link{CMBDataFrame}}
#' @param cmbdf2 a \code{\link{CMBDataFrame}}
#'
#' @examples
#' a <- CMBDataFrame(nside = 2, ordering = "ring", coords = "cartesian")
#' b <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
#' areCompatibleCMBDFs(a,b)
#'
#' suppressMessages(areCompatibleCMBDFs(a,b))
#'
#' @export
areCompatibleCMBDFs <- function(cmbdf1, cmbdf2)
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
    reasons <- paste0(reasons, "ordering mismatch (ordering1 = ", ordering(cmbdf1),
                      ", ordering2 = ", ordering(cmbdf2), ")\n")
  }
  if (!pix)
  {
    reasons <- paste0(reasons, "pixels mismatch pix(cmbdf1) != pix(cmbdf2)\n")
  }

  if ( !(ns && ord && pix) )
  {
    message(reasons)
    return(FALSE)
  }
  else
  {
    return(TRUE)
  }
}




#'Get the maximum distance between all points
#'in a \code{\link{CMBDataFrame}}
#'
#'@param cmbdf a CMBDataFrame object
#'
#'@export
maxDist.CMBDataFrame <- function(cmbdf)
{
  coords(cmbdf) <- "cartesian"
  return(maxDist_internal(cmbdf))
}






#' as.CMBDataFrame
#'
#' Safely converts a \code{\link{data.frame}} to a CMBDataFrame. The
#' rows of the data.frame are assumed to be in the HEALPix order
#' given by \code{ordering}, and at the HEALPix resolution given
#' by \code{nside}. Coordinates, if present, are checked to correspond
#' to HEALPix pixel centers. The coordinates must be named either x,y,z
#' (cartesian) or theta, phi (spherical colatitude and longitude respectively).
#'
#' @param df Any \code{data.frame} whose rows are in HEALPix order
#' @param ordering character string that specifies the ordering scheme
#' ("ring" or "nested")
#' @param nside an integer that specifies the Nside (resolution)
#' HEALPix parameter
#' @param spix a vector that specifies the HEALPix pixel index
#' corresponding to each row of \code{df}. If \code{spix} is left blank and
#' \code{df} is a \code{data.frame}, then \code{df} is assumed to contain data
#' for every pixel at resolution parameter \code{nside} (the full sky).
#' However, if \code{spix} is left blank and \code{df} is a \code{CMBDataFrame},
#' then \code{spix} is set equal to \code{pix(df)}
#'
#' @return A CMBDataFrame
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

    if (!missing(ordering) && ordering(df) != tolower(ordering) ) {
      stop(gettextf("Since '%s' is a CMBDataFrame already,
                    ordering should be unspecified or match the ordering
                    attribute of '%s'",
                    deparse(substitute(df)),
                    deparse(substitute(df))))
    }

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



#' Check if an object is of class CMBDat
#'
#' @param cmbdf Any R object
#'
#' @return TRUE if \code{cmbdf} is a CMBDat object, otherwise FALSE
#'
#' @export
is.CMBDat <- function(cmbdf)
{
  identical(as.numeric(sum(class(cmbdf) == "CMBDat")), 1)
}







#' Geodesic area covered by a \code{\link{CMBDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{cmbdf}
#'
#'@param cmbdf a CMBDataFrame
#'
#'@return the sum of the areas of all pixels (rows) in cmbdf
#'
#'@examples
#' ## At low resolution, a few data points can
#' ## occupy a large pixel area, e.g.:
#' cmbdf <- CMBDataFrame(nside = 1, spix = c(1,2,3))
#' pix(cmbdf)
#' geoArea(cmbdf) # pi = 1/4*(surface area of unit sphere)
#' plot(cmbdf, size = 5, hp.boundaries = 1)
#'
#'@export
geoArea.CMBDataFrame <- function(cmbdf)
{
  nside <- rcosmo:::nside(cmbdf)
  if ( !is.numeric(nside) ) stop("problem with cmbdf nside attribute")
  return(pi/(3*nside^2)*nrow(cmbdf))
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
#'@param cmbdf A CMBDataFrame.
#'@param new.coords Specifies the new coordinate system ("spherical" or "cartesian")
#'if a change of coordinate system is desired.
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
#' coords(df)
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' df
#'
#'
#'
#'@export
coords.CMBDataFrame <- function( cmbdf, new.coords )
{
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
      ns <- nside(cmbdf)
      sp <- pix(cmbdf)

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

      crds <- rcosmo:::pix2coords_internal(nside = ns, nested = nest,
                                           cartesian = cart, spix = sp)[,1:nc]
      crds <- as.data.frame(crds)
      names(crds) <- nam

      cmbdf <- rcosmo:::cbind.CMBDataFrame(crds, cmbdf)
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
      crds <- rcosmo:::car2sph(crds)
      other <- cmbdf[,-c(x.i, y.i, z.i), drop = FALSE]
      cmbdf <- rcosmo:::cbind.CMBDataFrame(crds, other)
      attr(cmbdf, "coords") <- "spherical"
    }
    else if ( new.coords == "cartesian" )
    {
      # convert to cartesian

      theta.i <- which(names(cmbdf) == "theta")
      phi.i <- which(names(cmbdf) == "phi")

      crds <- cmbdf[,c(theta.i, phi.i)]
      crds <- rcosmo:::sph2car(crds)
      other <- cmbdf[,-c(theta.i, phi.i), drop = FALSE]
      cmbdf <- rcosmo:::cbind.CMBDataFrame(crds, other)
      attr(cmbdf, "coords") <- "cartesian"
    }

    attr(cmbdf, "coords") <- new.coords
    return(cmbdf)
  }
}



#' Assign new coordinate system to a \code{\link{CMBDataFrame}}
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
`coords<-.CMBDataFrame` <- function(cmbdf,...,value) {
  return(coords(cmbdf, new.coords = value))
}







#' Plot CMB Data
#'
#' This function produces a plot from a \code{\link{CMBDataFrame}}.
#'
#'@param cmbdf A \code{\link{CMBDataFrame}}.
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
#'\code{length(labels)} must equal \code{nrow(cmbdf)}.
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
#' filename <- "CMB_map_smica1024.fits"
#' sky <- CMBDataFrame(filename)
#' plot(sky, sample.size = 800000)
#'
#'@export
plot.CMBDataFrame <- function(cmbdf, intensities = "I",
                              add = FALSE, sample.size,
                              type = "p", size = 1, box = FALSE,
                              axes = FALSE, aspect = FALSE,
                              col, back.col = "black", labels,
                              hp.boundaries = 0, hpb.col = "gray",
                              ...)
{
  if ( !missing(sample.size) )
  {
    spix <- sample(pix(cmbdf), sample.size)
    cmbdf <- cmbdf[pix(cmbdf) %in% spix,]
  }

  if (is.null(coords(cmbdf)))
  {
    coords(cmbdf) <- "cartesian"
  }

  if (missing(col))
  {
     col <- tryCatch(colscheme(cmbdf[,intensities, drop = TRUE],
                           rcosmo:::breaks1024, rcosmo:::colmap),
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
  cmbdf.xyz <- coords(cmbdf, new.coords = "cartesian")

  ## Do the plotting
  if ( !add )
  {
    rgl::open3d()
    rgl::bg3d(back.col)
  }

  if ( missing(labels) )
  {
    rgl::plot3d(cmbdf.xyz$x, cmbdf.xyz$y, cmbdf.xyz$z,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }
  else
  {
    rgl::text3d(cmbdf.xyz$x, cmbdf.xyz$y, cmbdf.xyz$z, labels,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }

  if ( hp.boundaries > 0 )
  {
    plotHPBoundaries(nside = hp.boundaries, col = hpb.col)
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
#'@param cmbdf a CMBDataFrame.
#'@param intensities the name of a column specifying CMB intensities
#'(or potentially another numeric quantity of interest)
#'
#'@return
#'A summary
#'
#'@export
summary.CMBDataFrame <- function(cmbdf, intensities = "I")
{
  ans <- list(intensities = summary(cmbdf[,intensities, drop = TRUE]))

  if ( is.null(coords(cmbdf)) )
  {
    ans$coords <- "HEALPix only"
  }
  else
  {
    ans$coords <- coords(cmbdf)
  }

  if ( is.null(window(cmbdf)) )
  {
    ans$window <- "full sky"
  }
  else
  {
    ans$window <- window(cmbdf)
  }

  if ( is.null(resolution(cmbdf)) )
  {
    ans$resolution <- "unknown"
  }
  else
  {
    ans$resolution <- resolution(cmbdf)
  }

  ans$ordering <- ordering(cmbdf)
  ans$nside <- nside(cmbdf)
  ans$pix <- pix(cmbdf)
  ans$n <- nrow(cmbdf)
  ans$area <- geoArea(cmbdf)
  ans$method <- header(cmbdf)[grepl("METHOD  =", header(cmbdf))]

  class(ans) <- "summary.CMBDataFrame"
  return(ans)
}


#'Print a summary of a CMBDataFrame
#'
#'@param x a \code{summary.CMBDataFrame} object, i.e.,
#'the output of \code{\link{summary.CMBDataFrame}}
#'
#'@export
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















#' Print CMB Data
#'
#' This function neatly prints the contents of a CMB Data Frame.
#'
#'@param cmbdf a CMB Data Frame.
#'@param ... arguments passed to \code{\link{print.tbl_df}}
#'
#'@return
#'Prints contents of the CMB data frame to the console.
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' print(df)
#' df
#'
#'@export
print.CMBDataFrame <- function(cmbdf,...)
{
  cat("A CMBDataFrame\n")
  print(tibble::as.tibble(cmbdf), ...)
}





