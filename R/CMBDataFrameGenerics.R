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







#' Area of a \code{\link{CMBDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{cmbdf}
#'
#'@param cmbdf a CMBDataFrame
#'
#'@return the sum of the areas of all pixels (rows) in cmbdf
#'
#'@export
area.CMBDataFrame <- function(cmbdf)
{
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

    if ( is.null(attr(cmbdf, "coords")) )
    {
      stop(paste("(development stage) not yet implemented",
                 "pix2coords in coords function"))
    }

    # Make sure that new.coords doesn't match current coords
    if ( attr(cmbdf, "coords") == new.coords )
    {
      # Nothing to do
    }
    else if ( new.coords == "spherical" )
    {
      # Convert to spherical
      n <- ncol(cmbdf)
      other.names <- names(cmbdf)[-c(which(names(cmbdf) == "x"),
                                     which(names(cmbdf) == "y"),
                                     which(names(cmbdf) == "z"))]

      xyz <- cmbdf[,c("x", "y", "z")]
      others <- cmbdf[, other.names]

      cmbdf[,1:2] <- car2sph(xyz)
      cmbdf[,3:(n-1)] <- others
      cmbdf[,n] <- NULL
      names(cmbdf) <- c("theta", "phi", other.names)
    }
    else if ( new.coords == "cartesian" )
    {
      # convert to cartesian

      n <- ncol(cmbdf)
      other.names <- names(cmbdf)[-c(which(names(cmbdf) == "theta"),
                                     which(names(cmbdf) == "phi"))]

      sph <- cmbdf[,c("theta", "phi")]
      others <- cmbdf[, other.names]

      cmbdf[,1:3] <- sph2car(sph)
      cmbdf[,4:(n+1)] <- others
      names(cmbdf) <- c("x","y","z", other.names)
    }

    attr(cmbdf, "coords") <- new.coords
    return(cmbdf)
  }
}



#' Assign new coordinate system to CMBDataFrame
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
                              col, back.col, ...)
{

  if ( !missing(sample.size) )
  {
    spix <- sample(pix(cmbdf), sample.size)
    cmbdf <- cmbdf[spix,]
  }

  ## Stored data is used to make colours if col is missing
  if ( missing(col) )
  {
    if ( nside(cmbdf) == 1024 )
    {
      col <- ifelse(missing(sample.size), CMBcols1024, CMBcols1024[spix])

      warning(paste("(development stage) the colour map for",
              "nside = 1024 may not be ideal"))

    }
    else if ( nside(cmbdf) == 2048 )
    {

      ## The following code must be replaced to work with nside = 2048
      col <- ifelse(missing(sample.size), CMBcols1024, CMBcols1024[spix])

      warning(paste("(development stage) the colour map used was not",
              "generated for nside = 2048"))
    }
    else
    {
      col <- "blue"
    }
  }

  ## Change coordinates if necessary
  coords <- coords(cmbdf)

  try(if(coords != "spherical" && coords != "cartesian")
    stop("Coordinates must be spherical or cartesian"))

  if (coords == "spherical") {
    cmbdf.xyz <- rcosmo::sph2car(cmbdf[,c("theta","phi")])
  } else {
    # Else coords are already cartesian
    cmbdf.xyz <- data.frame(x = cmbdf$x, y = cmbdf$y, z = cmbdf$z)
  }



  ## Do the plotting
  if ( !missing(back.col) )
  {
    rgl::open3d()
    rgl::bg3d(back.col)
  }
  rgl::plot3d(cmbdf.xyz, col = col, type = type, size = size,
              box = box, axes = axes, add = add, aspect = aspect, ...)
}









#' Summarise CMB Data
#'
#' This function produces a summary from a CMB Data Frame.
#'
#'@param cmbdf a CMB Data Frame.
#'
#'@return
#'A summary of the CMB data.
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' summary(df)
#'
#'@export
summary.CMBDataFrame <- function(cmbdf)
{
  if ( sum(names(df) == "I") %>% as.numeric() %>% identical(1) )
  {
    ans <- list(intensities = summary(df$I))
  }

  ans[["coords"]] <- ifelse(is.null(coords(cmbdf)),
                            "HEALPix only",coords(cmbdf))
  ans[["ordering"]] <- ordering(cmbdf)
  ans[["nside"]] <- nside(cmbdf)
  ans[["window"]] <- "The CMBWindow class is under development"

  ans
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
  print(tibble::as.tibble(cmbdf),...)
}





