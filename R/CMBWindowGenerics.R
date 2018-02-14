
#' Coordinate system from a \code{\link{CMBWindow}}
#'
#' This function returns the coordinate system used in a
#' \code{\link{CMBWindow}}. The coordinate system is either
#' "cartesian" or "spherical"
#'
#' If a new coordinate system is specified, using e.g.
#' new.coords = "spherical", the
#' coordinate system of the CMBWindow will be converted
#'
#'@param cmbdf a CMBWindow.
#'@param new.coords specifies the new coordinate system
#'("spherical" or "cartesian")
#'if a change of coordinate system is desired
#'
#'@return
#' If new.coords is unspecified, then the name of the coordinate system
#' of \code{win} is returned. Otherwise a new CMBWindow is returned
#' equivalent to \code{win} but having the desired change of coordinates
#'
#'@examples
#' df <- CMBDataFrame("CMB_map_smica1024.fits", sample.size = 800000)
#' coords(df)
#' coords(df, new.coords = "cartesian")
#'
#'@export
coords.CMBWindow <- function( win, new.coords )
{
  # If new.coords argument is missing then return the coordinate type
  if ( missing(new.coords) )
  {
    return(attr(win, "coords"))
  }
  else
  {
    new.coords <- as.character(tolower(new.coords))

    # Make sure that new.coords doesn't match current coords
    if ( attr(win, "coords") == new.coords )
    {
      # Nothing to do
    }
    else if ( new.coords == "spherical" )
    {
      # Convert to spherical
      win[,1:2] <- car2sph(win)
      names(win) <- c("theta", "phi")
      win[,3] <- NULL
      attr(win, "coords") <- "spherical"
    }
    else if ( new.coords == "cartesian" )
    {
      # Convert to cartesian
      win[,1:3] <- sph2car(win)
      names(win) <- c("x", "y", "z")
      attr(win, "coords") <- "cartesian"
    }

    return(win)
  }
}




#'Assign new coordinate system to CMBWindow
#'@export
`coords<-.CMBWindow` <- function(win,...,value) {
  value <- tolower(value)
  if (coords(win) == value)
  {
    return(win)

  } else {

    return(coords(win, new.coords = value))
  }
}
