#'Get the maximum distance between any two points in a CMBWindow
#'
#'@param win a CMBWindow object
#'
#'@export
maxDist <- function(win)
{
  if ( !is.CMBWindow(win) ) stop("Argument 'win' must be a CMBWindow")

  return(attr(win, "maxDist"))
}


#' Get the spherical area inside a CMBWindow
#'
#' @param win a CMBWindow
#'
#' @return Tthe spherical area inside win
#'
#'@export
area <- function(win)
{
  if ( !is.CMBWindow(win) ) stop("Argument 'win' must be a CMBWindow")

  return(attr(win, "area"))
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
