
#'Triangulate a polygonal \code{\link{CMBWindow}}
#'
#'@param win a CMBWindow object
#'
#'@return a list of CMBWindow polygons or minus.polygons,
#'each having 3 vertices and representing a triangle.
#' If winType of \code{win} does not include
#' "minus" then these triangles have
#' pairwise disjoint interiors and their union
#' is equal to the original polygon,
#' \code{win}.
#' Otherwise, if winType of \code{win}
#' does include "minus" the triangles are
#' the same as for the non-minus type above, but have "minus" types.
#'
#'@examples
#'
#'## Example 1
#'
#' win <- CMBWindow(theta = c(2*pi/3,3*pi/4,3*pi/4, 2*pi/3),
#'                  phi = c(pi/4,pi/4,pi/3,pi/3))
#' win
#' plot(win)
#' win1 <- triangulate(win)
#' win1
#' summary(win1[[1]])
#' plot(win1[[1]], add= FALSE, col="green")
#' plot(win1[[2]], col="blue")
#'
#' ## Example 2: triangilation minus-type polygon
#'
#' win <- CMBWindow(theta = c(pi/5,pi/3,pi/4, pi/3, pi/5),
#'                  phi = c(pi/5,pi/5, pi/4 ,pi/3,pi/3), set.minus =TRUE)
#' win
#' plot(win)
#' summary(win)
#' win1 <- triangulate(win)
#' win1
#' plot(win1[[1]], add= FALSE, col="green")
#' plot(win1[[2]], col="blue")
#' plot(win1[[3]], col="yellow")
#' summary(win1[[1]])
#' summary(win1[[2]])
#' summary(win1[[3]])
#'
#'@export
triangulate <- function(win)
{
  if ( !is.CMBWindow(win) ) stop("'win' must be a CMBWindow")

  if ( winType(win) != "polygon" && winType(win) != "minus.polygon" )
  {
    stop("The winType of 'win' must be 'polygon' or 'minus.polygon'")
  }

  coords(win) <- "cartesian"

  triangles <- list()
  i <- 1
  reset.i <- TRUE

  while( nrow(win) > 3 )
  {
    n <- nrow(win)
    i <- i + 1

    if ( reset.i )
    {
      i <- 1
      reset.i <- FALSE
    }

    if (i == nrow(win))
    {
      break;
      stop(paste("Triangulation failed. Ensure the polygon",
                 "is oriented counter-clockwise.",
                 "If you are sure the polygon is convex then",
                 "perhaps try assume.convex = TRUE as an",
                 "argument to CMBWindow"))
    }

    i.1 <- i
    i.2 <- i %% n + 1       # i + 1 cyclic
    i.3 <- (i+1) %% n + 1   # i + 2 cyclic

    V1 <- as.numeric(win[i.1,])
    V2 <- as.numeric(win[i.2,])
    V3 <- as.numeric(win[i.3,])

    # Check if the vertex at V2 is a concave-up corner.
    tri <- matrix(c(V1,V2,V3), nrow = 3, byrow = TRUE)
    if ( det(tri) > 0 ) # V1 cross V2 dot V3
    {
      # Check if there are any vertices inside the triangle (v1,v2,v3)
      # other than V1,V2,V3 themselves
      vertex.inside <- FALSE
      for ( j in seq(1,n)[-c(i.1, i.2, i.3)] )
      {
        Vj <- as.numeric(win[j,])
        if ( det(matrix(c(V1,V2,Vj), nrow = 3)) <= 0
             || det(matrix(c(V2,V3,Vj), nrow = 3)) <= 0
             || det(matrix(c(V3,V1,Vj), nrow = 3)) <= 0)
        {
          ## Vj is outside, do nothing
        }
        else
        {
          vertex.inside <- TRUE
          break;
        }
      }


      ## If there are no vertices inside then (V1,V2,V3) is an ear
      if ( !vertex.inside )
      {
        tri <- data.frame(tri)
        names(tri) <- c("x","y","z")
        tri <- CMBWindow(tri, assume.convex = TRUE)
        attr(tri, "winType") <- winType(win) # in case set.minus = TRUE
        triangles[[length(triangles)+1]] <- tri

        win <- win[-i.2,]
        reset.i <- TRUE
      }
    }

  }

  ## while loop has exited so the remaining vertices form a triangle
  attr(win, "assumedConvex") <- TRUE
  triangles[[length(triangles)+1]] <- win

  return(triangles)
}










#'Check if a \code{\link{CMBWindow}} is assumed convex.
#'
#'Initially any \code{\link{CMBWindow}} is not assumed convex.
#'The assumedConvex
#' attribute can be change for any \code{\link{CMBWindow}}.
#'
#'@param win a CMBWindow object
#'@param assume.convex optionally change the assumedConvex
#'attribute to TRUE or FALSE
#'
#'@examples
#'
#'win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#'assumedConvex(win1)
#'win2 <- assumedConvex(win1, assume.convex = TRUE)
#'assumedConvex(win2)
#'assumedConvex(win1) <- TRUE
#'assumedConvex(win1)
#'
#'@export
assumedConvex <- function(win, assume.convex)
{
  if ( !is.CMBWindow(win) ) stop("Argument 'win' must be a CMBWindow")

  if ( !missing(assume.convex) )
  {
    if (!is.logical(assume.convex))
    {
      stop("assume.convex must be logical TRUE or FALSE")
    }

    attr(win, "assumedConvex") <- assume.convex

    if ( winType(win) == "disc" || winType(win) == "minus.disc" )
    {
      warning(paste("Changing the assumedConvex attribute of a disc",
                    "is strange and may lead to undefined behaviour"))
    }

    return(win)
  }
  else
  {
    return(attr(win, "assumedConvex"))
  }
}


#' Change the \code{\link{assumedConvex}} boolean of a \code{\link{CMBWindow}}
#'
#' @keywords internal
#'
#' @export
`assumedConvex<-` <- function(win, ..., value)
{
  return(assumedConvex(win, assume.convex = value))
}










#'Get/change winType
#'
#'Get/change the winType (polygon or disk) of a \code{\link{CMBWindow}}.
#'If \code{new.type} is missing then the \code{winType} of win
#'is returned. Otherwise, a new window is returned with \code{winType}
#'equal to \code{new.type}. If you want to change the
#'winType of \code{win} directly, then use \code{\link{winType<-}}, see
#'the examples below.
#'
#'@param win a \code{CMBWindow} object or a list of such
#'@param new.type optionally specify a new type. Use this to
#'change between "polygon" and "minus.polygon" or to change
#'between "disc" and "minus.disc"
#'
#'@return If \code{new.type} is missing then the \code{winType} of win
#'is returned. Otherwise a new window is returned with \code{winType}
#'equal to \code{new.type}
#'
#'
#'@examples
#'
#' win <- CMBWindow(theta = c(pi/2,pi/2,pi/3, pi/3), phi = c(0,pi/3,pi/3,0))
#' winType(win)
#'
#' win1 <- CMBWindow(x=0,y=3/5,z=4/5,r=0.8)
#' winType(win1)
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian",
#'                       ordering = "nested")
#' cmbdf.win1 <- window(cmbdf, new.window = win1)
#' plot(cmbdf.win1)
#'
#'
#' winType(win1) <- "minus.disc"
#' winType(win1)
#' cmbdf <- CMBDataFrame(nside = 64, coords = "cartesian",
#'                       ordering = "nested")
#' cmbdf.win1 <- window(cmbdf, new.window = win1)
#' plot(cmbdf.win1)
#'
#'
#'@export
winType <- function(win, new.type)
{
  if ( !is.CMBWindow(win) )
  {
    if ( !is.list(win) || !all(sapply(win, is.CMBWindow)) )
    {
      stop("'win' must be a CMBWindow of list of CMBWindows")
    }

    return(sapply(win, winType))
  }

  if ( missing(new.type) )
  {
    return(attr(win, "winType"))
  }
  else
  {
    # Only change the winType if new.type is compatible
    if ( new.type == "disc" || new.type == "minus.disc" )
    {
      if ( contains("disc", winType(win)) )
      {
        attr(win, "winType") <- new.type
        return(win)
      }
      else
      {
        stop("cannot change from the specified type to a disc type")
      }
    } else if ( new.type == "polygon" || new.type == "minus.polygon" ) {

      if ( contains("polygon", winType(win)) )
      {
        attr(win, "winType") <- new.type
        return(win)
      }
      else
      {
        stop("cannot change from the specified type to a polygon type")
      }
    }
    else
    {
      stop("the specified new.type is invalid")
    }
  }
}


#' Assign new \code{\link{winType}} to a \code{\link{CMBWindow}}
#'
#' @keywords internal
#'
#' @seealso \code{\link{winType}}
#'
#' @examples
#'
#' win <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5, set.minus = TRUE)
#' winType(win)
#' winType(win) <- "disc"
#' winType(win)
#'
#' @export
`winType<-` <- function(win, ..., value)
{
  return(winType(win, new.type = value))
}





#'Get the maximum distance between all points
#'in a \code{\link{CMBWindow}}
#'
#'@param x A CMBWindow object.
#'
#'@return The maximum distance between window's points.
#'
#'@examples
#'
#' ## win is a equilateral spherical triangle which sides are pi/2
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' maxWindowDist(win)
#'
#'@export
maxWindowDist <- function(x)
{

  # Create temporary window in cartesian coordinates for dist and area
  if ( coords(x) == "cartesian" ) {

    win.xyz <- x

  } else {

    win.xyz <- sph2car(x)
  }

  # Calculate maximum distance
  max.dist <- switch(winType(x),
                     polygon = polygonMaxDist(win.xyz),
                     minus.polygon = pi,
                     disc = 2*as.numeric(x$r),
                     minus.disc = pi,
                     stop(paste("Could not determine window type",
                                "using rcosmo::winType")))

  return(max.dist)
}

## HELPER FUNCTION FOR maxWindowDist
polygonMaxDist <- function(win)
{
  max.dist <- 0
  for ( i in 1:(nrow(win) - 1) )
  {
    for ( j in (i+1):nrow(win) )
    {
      dist <- geoDist(win[i,], win[j,])
      if ( dist > max.dist ) max.dist <- dist
    }
  }

  return(max.dist)
}

