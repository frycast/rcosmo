



#' CMBWindow
#'
#' Add description here
#'
#' Add further details here
#'
#'@param ... arguments that must be labelled either x, y, z
#'(cartesian) or theta, phi (spherical, colatitude and latitude respectively).
#'Alternatively, a single data.frame may be passed in with columns labelled
#'x, y, z or theta, phi.
#'The rows correspond to clockwise ordered
#'vertices defining a spherical polygon on the surface of the unit sphere.
#'There must be at least 3 rows (i.e., spherical bigons are excluded).
#''Clockwise' is understood from a perspective outside the sphere looking
#'in towards the origin.
#'
#'@return
#'
#'@examples
#'win <- CMBWindow(theta = c(0,1,2), phi = c(0.5, 1.2))
#'
#'@export
CMBWindow <- function(...) {

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

  if ( nrow(window) < 3 ) stop(paste("the window must have",
                                     "at least 3 vertices"))

  # Create temporary window in cartesian coordinates for dist and area
  if ( coords == "cartesian" ) {

    window.xyz <- window

  } else {

    window.xyz <- sph2car(window)
  }

  # Calculate maximum distance within the spherical polygon
  max.dist <- 0
  for ( i in 1:(nrow(window.xyz) - 1) )
  {
    for ( j in (i+1):nrow(window.xyz) )
    {
      dist <- geoDist(window.xyz[i,], window.xyz[j,])
      if ( dist > max.dist ) max.dist <- dist
    }
  }

  # Calculate the area of the spherical polygon
  area <- sphericalArea(window.xyz)


  class(window) <- c("CMBWindow", "data.frame")
  attr(window, "coords") <- coords
  attr(window, "maxDist") <- as.numeric(max.dist)
  attr(window, "area") <- area

  return(window)
}








## HELPER FUNCTION to calculate area of spherical polygon
# win must be a data.frame in cartesian coordinates:
sphericalArea <- function(win)
{
  n <- nrow(win)

  # Each iteration finds the angle at vertex i + 1 (A2)
  angles <- vector(mode = "numeric", length = n)
  for ( i in 1:n )
  {
    # Get vertex coordinates
    A1 <- as.numeric(win[ i, ])
    A2 <- as.numeric(win[ i %% n + 1, ])     # i + 1 (cyclic)
    A3 <- as.numeric(win[ (i+1) %% n + 1, ]) # i + 2 (cyclic)

    # Calculate normal vectors to define planes through origin
    n1 <- vector_cross(A1, A2)
    n2 <- vector_cross(A2, A3)

    # Moduli of normal vectors
    mod.n1 <- sqrt(sum(n1*n1))
    mod.n2 <- sqrt(sum(n2*n2))

    # Angle between normal vectors is equal to angle at A2
    theta <- acos( -sum(n1*n2)/(mod.n1*mod.n2) )

    # Decide if angleat A2 is obtuse or acute using scalar triple product
    if ( det(matrix(c(A1, A2, A3), nrow = 3)) <= 0 )
    {
      angles[i] <- theta
    } else {
      angles[i] <- 2*pi - theta
    }
  }

  # Use Gauss-Bonnet Theorem for area of spherical polygon
  return(sum(angles) - (n-2)*pi)
}
