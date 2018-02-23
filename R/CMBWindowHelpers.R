#'Get the maximum distance between all points
#'in a \code{\link{CMBWindow}}
#'
#'@param win a CMBWindow object
#'
#'@export
maxDist <- function(win)
{
  if ( !is.CMBWindow(win) ) stop("Argument 'win' must be a CMBWindow")

  # Create temporary window in cartesian coordinates for dist and area
  if ( coords(win) == "cartesian" ) {

    win.xyz <- win

  } else {

    win.xyz <- sph2car(win)
  }

  # Calculate maximum distance
  max.dist <- switch(winType(win),
                     polygon = polygonMaxDist(win.xyz),
                     minus.polygon = pi,
                     disc = 2*as.numeric(win$r),
                     minus.disc = pi,
                     stop(paste("Could not determine window type",
                                "using rcosmo::winType")))

  return(max.dist)
}

## HELPER FUNCTION FOR maxDist
polygonMaxDist <- function(win)
{
  max.dist <- 0
  for ( i in 1:(nrow(win.xyz) - 1) )
  {
    for ( j in (i+1):nrow(win.xyz) )
    {
      dist <- geoDist(win.xyz[i,], win.xyz[j,])
      if ( dist > max.dist ) max.dist <- dist
    }
  }

  return(max.dist)
}







#'Get the type (polygon or disk) of a \code{\link{CMBWindow}}
#'
#'@param win a CMBWindow object
#'
#'@export
winType <- function(win)
{
  if ( !is.CMBWindow(win) )
  {
    stop("'win' must be a CMBWindow")
  }

  return(attr(win, "winType"))
}



#' Check if an object is a CMBWindow
#'
#' @param win any object
#'
#' @return TRUE or FALSE depending if win is a CMBWindow
#'
#'@export
is.CMBWindow <- function(win)
{
  return(identical(as.numeric(sum(class(win) == "CMBWindow")), 1))
}
