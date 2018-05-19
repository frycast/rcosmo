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

  if ( !all(sapply(args, rcosmo::is.CMBDataFrame)) )
  {
    stop("rbind.CMBDataFrame requires all arguments to be CMBDataFrames")
  }

  if ( !all(sapply(args, rcosmo::areCompatibleCMBDFs, cmbdf2 = args[[1]])) )
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
  wins <- sapply(args, rcosmo::window)

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
#' a <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
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
  crd <- identical(coords(cmbdf1), coords(cmbdf2))

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
  if (!crd)
  {
    reasons <- paste0(reasons, "coords mismatch (coords1 = ", coords(cmbdf1),
                      ", coords2 = ", coords(cmbdf2), ")")
  }

  if ( !(ns && ord && crd) )
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







#' Geodesic area covered by a \code{\link{CMBDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{cmbdf}
#'
#'@param cmbdf a CMBDataFrame
#'
#'@return the sum of the areas of all pixels (rows) in cmbdf
#'
#'@export
geoArea.CMBDataFrame <- function(cmbdf)
{
  nside <- nside(cmbdf)
  if ( !is.numeric(nside) ) stop("problem with cmbdf nside parameter")
  return(4*pi/(12*nside^2)*nrow(cmbdf))
}









#' Coordinate system from a CMBDataFrame
#'
#' This function returns the coordinate system used in a CMBDataFrame.
#' The coordinate system is either "cartesian" or "spherical"
#'
#' If a new coordinate system is specified, using e.g. new.coords = "spherical", the
#' coordinate system of the CMBDataFrame will be converted.
#'
#'@param cmbdf a CMBDataFrame.
#'@param new.coords specifies the new coordinate system ("spherical" or "cartesian")
#'if a change of coordinate system is desired.
#'
#'@return
#' If new.coords is unspecified, then the name of the coordinate system
#' of \code{cmbdf} is returned. Otherwise a new CMBDataFrame is returned
#' equivalent to \code{cmbdf} but having the desired change of coordinates
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' coords(df)
#' coords(df, new.coords = "cartesian")
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
      crds <- rcosmo::car2sph(crds)
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
      crds <- rcosmo::sph2car(crds)
      other <- cmbdf[,-c(theta.i, phi.i), drop = FALSE]
      cmbdf <- rcosmo:::cbind.CMBDataFrame(crds, other)
      attr(cmbdf, "coords") <- "cartesian"
    }

    attr(cmbdf, "coords") <- new.coords
    return(cmbdf)
  }
}



#' Assign new \code{\link{coords}} system to \code{\link{CMBDataFrame}}
#' @export
`coords<-.CMBDataFrame` <- function(cmbdf,...,value) {
  return(coords(cmbdf, new.coords = value))
}








#### CURRENTLY THE DATA USED IN THE PLOT FUNCTION IS TOO LARGE FOR CRAN ##
#' Plot CMB Data
#'
#' This function produces a plot from a CMB Data Frame.
#'
#'@param cmbdf a CMB Data Frame with either spherical or cartesian coordinates.
#'@param add if TRUE then this plot will be added to any existing plot.
#'Note that if \code{back.col} (see below) is specified then a new plot
#'window will be opened and \code{add = TRUE} will have no effect
#'@param sample.size optionally specifies the size of a simple random
#'sample to take before plotting. This can make the plot less
#'computationally intensive
#'@param type a single character indicating the type of item to plot.
#'Supported types are: 'p' for points, 's' for spheres, 'l' for lines,
#''h' for line segments from z = 0, and 'n' for nothing.
#'@param size the size of plotted points
#'@param box whether to draw a box
#'@param axes whether to draw axes
#'@param aspect either a logical indicating whether to adjust the
#'aspect ratio, or a new ratio.
#'@param col specify the colour(s) of the plotted points
#'@param back.col optionally specifies the background colour of
#'the plot. This argument is passed to rgl::bg3d.
#'@param labels optionally specify a vector of labels to plot,
#'such as words or vertex indices. If this is specified then
#'\code{rgl::text3d} is used instead of \code{rgl::plot3d}. Then
#'\code{length(labels)} must equal \code{nrow(cmbdf)}
#'@param ... arguments passed to rgl::plot3d
#'
#'@return
#'A plot of the CMB data
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#' plot(df, sample.size = 800000)
#'
#'@export
plot.CMBDataFrame <- function(cmbdf, add = FALSE, sample.size,
                              type = "p", size = 1, box = FALSE,
                              axes = FALSE, aspect = FALSE,
                              col, back.col = "black", labels, ...)
{

  if (is.null(coords(cmbdf)))
  {
    # When the coods are HEALPix only there is some issue
    # with sampling in the plot function then converting to cartesian.
    # The plot function crashes / is slow.
    stop("(development stage) cannot plot when coords are NULL")
  }

  if ( !missing(sample.size) )
  {
    spix <- sample(pix(cmbdf), sample.size)
    cmbdf <- cmbdf[pix(cmbdf) %in% spix,]
  }

  if (missing(col))
  {
    col <- rcosmo:::colmap[cut(cmbdf$I, length(colmap))]
  }


  # ## Stored data is used to make colours if col is missing
  # if ( missing(col) )
  # {
  #   if ( nside(cmbdf) == 1024 )
  #   {
  #     if (missing(sample.size))
  #     {
  #       col <- rcosmo:::CMBcols1024
  #     }
  #     else
  #     {
  #       col <-  CMBcols1024[spix]
  #       stop("(development stage) colours not assigned to sample")
  #     }
  #
  #     warning(paste("(development stage) the colour map for",
  #             "nside = 1024 may not be ideal"))
  #
  #   }
  #   else if ( nside(cmbdf) == 2048 )
  #   {
  #
  #     stop("(development stage) colours not assigned")
  #
  #     ## The following code must be replaced to work with nside = 2048
  #     if (missing(sample.size))
  #     {
  #       col <- rcosmo:::CMBcols1024
  #     }
  #     else
  #     {
  #       col <-  CMBcols1024[spix]
  #     }
  #
  #     warning(paste("(development stage) the colour map used was not",
  #             "generated for nside = 2048"))
  #   }
  #   else
  #   {
  #     col <- "blue"
  #   }
  # }


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
    rgl::plot3d(cmbdf.xyz, col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }
  else
  {
    rgl::text3d(cmbdf.xyz$x, cmbdf.xyz$y, cmbdf.xyz$z, labels,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }
}









#' Summarise a \code{\link{CMBDataFrame}}
#'
#' This function produces a summary from a CMBDataFrame.
#'
#'@param cmbdf a CMBDataFrame.
#'
#'@return
#'A summary
#'
#'@export
summary.CMBDataFrame <- function(cmbdf)
{
  ans <- list(intensities = summary(cmbdf$I))

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
  if ( x$window != "full sky" )
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
  print(tibble::as.tibble(cmbdf), ...)
}





