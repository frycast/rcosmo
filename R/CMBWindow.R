



#' CMBWindow
#'
#' Create a CMBWindow: Either a polygon or a disc type
#'
#'If \code{r} is unspecified then the rows of \code{...} correspond to
#'counter-clockwise ordered vertices defining a spherical polygon
#'on the unit sphere.
#'In this case, there must be at least 3 rows (vertices).
#'On the other hand,
#'if \code{r} is specified then \code{...} must have just one row, and this
#'row is taken to be the center of a disc of radius \code{r}
#'
#'@param ... these arguments are compulsory and must be labelled either x, y, z
#'(cartesian) or theta, phi (spherical, colatitude and longitude respectively).
#'Alternatively, a single data.frame may be passed in with columns labelled
#'x, y, z or theta, phi.
#'@param r if a disc type window is required then this specifies the
#'radius of the disc
#'@param set.minus when \code{TRUE} the window will be the unit
#'sphere minus the window specified
#'@param assume.convex when \code{TRUE} the window is assumed to be convex
#'resulting in a faster computation time when the window is used with functions
#'such as \code{\link{subWindow}}. This argument is irrelevant when the window
#'is not a polygon
#'
#'@return
#'
#'@examples
#'win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#'
#'@export
CMBWindow <- function(..., r, set.minus = FALSE, assume.convex = FALSE) {

  args <- list(...)

  # theta AND phi WERE PASSED TO '...'
  if ( all( c("theta", "phi") %in% names(args) ) )
  {
    if (length(args) != 2) stop(paste("please only specify theta and phi,",
                                "extra arguments were given to '...'"))

    window <- data.frame(theta = args[["theta"]], phi = args[["phi"]])
    coords <- "spherical"

  # x, y, z WERE PASSED TO '...'
  } else if ( all( c("x", "y", "z") %in% names(args) ) ) {

    if (length(args) != 3) stop(paste("please only specify x, y and z,",
                                "extra arguments were given to '...'"))

    coords <- "cartesian"
    window <- data.frame(x = args[["x"]], y = args[["y"]], z = args[["z"]])

  # A data.frame WAS PASSED TO '...'
  } else if ( length(args) == 1 ) {

    df <- args[[1]]
    if ( !is.data.frame(df) ) stop(paste("only one argument was passed to",
                                        "'...' and it was not a data.frame"))

    if ( all( c("theta", "phi") %in% names(df) ) ) {

      coords <- "spherical"
      window <- data.frame(theta = args[["theta"]], phi = args[["phi"]])

    } else if ( all( c("x", "y", "z") %in% names(df) ) ) {

      coords <- "cartesian"
      window <- data.frame(x = args[["x"]], y = args[["y"]], z = args[["z"]])

    } else {
      stop(paste("the data.frame does not have columns labelled 'x','y','z'",
                  "or 'theta', 'phi'"))
    }

  # INCORRECT ARGUMENTS PASSED TO ...
  } else {
    stop(paste("must specify either x, y and z or theta and phi.",
               "\nOr else pass in a data.frame containing those."))
  }

  if ( missing(r) && nrow(window) < 3 )
  {
    stop(paste("'r' is unspecified so the window is a polygon,",
               "but a polygonal window must have",
                "at least 3 vertices"))
  }

  ## This is all we have to do for disc windows
  if ( !missing(r) )
  {
    if ( nrow(window) != 1 )
    {
      stop(paste("'r' is specified so the window is a disc,",
                 "but a disc window must have",
                 "just one center (one row)"))
    }
    else
    {
      window <- cbind(window, r = r)
    }
  }


  class(window) <- c("CMBWindow", "data.frame")
  attr(window, "coords") <- coords
  if ( missing(r) )
  {
    attr(window, "winType") <- ifelse(set.minus, "minus.polygon", "polygon")
    attr(window, "assumedConvex") <- assume.convex
  }
  else
  {
    attr(window, "winType") <- ifelse(set.minus, "minus.disc", "disc")
    attr(window, "assumedConvex") <- TRUE
  }

  return(window)
}









