#' CMBWindow class.
#'
#' The function \code{CMBWindow} creates objects of class \code{CMBWindow}.
#' It is either a polygon or a disc type.
#'
#'If \code{r} is unspecified then the rows of \code{...} correspond to
#'counter-clockwise ordered vertices defining a spherical polygon
#'lying entirely within one open hemisphere on the unit sphere.
#'Counter-clockwise is understood from the perspective outside the
#'sphere, facing the hemisphere that contains the polygon, looking
#'toward the origin. Note that there must be at least 3 rows (vertices)
#'to define a polygon (we exlude bygons).
#'On the other hand,
#'if \code{r} is specified then \code{...} must specify just one row, and this
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
#'@examples
#' win <- CMBWindow(theta = c(pi/2,pi/2,pi/3, pi/3), phi = c(0,pi/3,pi/3,0))
#' plot(win)
#'
#'
#' ## Create a disc type window
#' win1<- CMBWindow(x=0,y=3/5,z=4/5,r=0.8, set.minus =TRUE)
#' plot(win1)
#'
#'
#' ## Apply a disc type window to CMBDataFrame
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian", ordering = "nested")
#' window(cmbdf) <- CMBWindow(x=0,y=3/5,z=4/5,r=0.8, set.minus =TRUE)
#' plot(cmbdf)
#'@export
CMBWindow <- function(..., r, set.minus = FALSE, assume.convex = FALSE) {

  args <- list(...)

  # theta AND phi WERE PASSED TO '...'
  if ( all( c("theta", "phi") %in% names(args) ) )
  {
    if (length(args) != 2) stop(paste("please only specify theta and phi,",
                                "extra arguments were given to '...'.",
                                "The arguments given to '...' were labelled",
                                paste(dQuote(names(args)), collapse = " and ")))

    window <- data.frame(theta = args[["theta"]], phi = args[["phi"]])
    coords <- "spherical"

  # x, y, z WERE PASSED TO '...'
  } else if ( all( c("x", "y", "z") %in% names(args) ) ) {

    if (length(args) != 3) stop(paste("please only specify x, y and z,",
                                "extra arguments were given to '...'.",
                                "The arguments given to '...' were labelled",
                                paste(dQuote(names(args)), collapse = " and ")))

    coords <- "cartesian"
    window <- data.frame(x = args[["x"]], y = args[["y"]], z = args[["z"]])

    if ( !isTRUE(all.equal(apply(window, 1, function(x) {sum(as.numeric(x)^2)}),
                  rep(1,nrow(window)))) )
    {
      warning(paste("One or more vertices in the CMBWindow",
                      "do not lie on the unit sphere"))
    }

  # A data.frame WAS PASSED TO '...'
  } else if ( length(args) == 1 ) {

    df <- args[[1]]
    if ( !is.data.frame(df) ) stop(paste("only one argument was passed to",
                                        "'...' and it was not a data.frame"))

    if ( is.CMBWindow(df) )
    {

      window <- df

    } else if ( all( c("theta", "phi") %in% names(df) ) ) {

      coords <- "spherical"
      window <- df[,c("theta","phi")]

    } else if ( all( c("x", "y", "z") %in% names(df) ) ) {

      coords <- "cartesian"
      window <- df[,c("x","y","z")]

      if ( !isTRUE(all.equal(apply(window, 1, function(x) {sum(as.numeric(x)^2)}),
                  rep(1,nrow(window)))) )
      {
        warning(paste("One or more vertices in the CMBWindow",
                      "do not lie on the unit sphere"))
      }

    } else {
      stop(paste("the data.frame does not have columns labelled 'x','y','z'",
                  "or 'theta', 'phi'"))
    }

  # INCORRECT ARGUMENTS PASSED TO ...
  } else {
    stop(paste("must specify either x, y and z or theta and phi.",
               "\nOr else pass in a data.frame containing those.",
               "The arguments given to '...' were labelled",
               paste(dQuote(names(args)), collapse = " and ")))
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









