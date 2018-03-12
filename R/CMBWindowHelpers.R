
#'Triangulate a polygonal \code{\link{CMBWindow}}
#'
#'@param win a CMBWindow object
#'
#'@return a list of CMBWindow polygons or minus.polygons,
#'each having 3 vertices and representing a triangle.
#'These triangles have pairwise disjoint interiors and their
#'union is equal to the original polygon, \code{win}.
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










#'Check if a \code{\link{CMBWindow}} is assumed convex
#'
#'@param win a CMBWindow object
#'@param assume.convex optionally change the assumedConvex
#'attribute to TRUE or FALSE
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
#' @export
`assumedConvex<-` <- function(win, ..., value)
{
  return(assumedConvex(win, assume.convex = value))
}



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
#'@param win a \code{CMBWindow} object or a list of such
#'@param new.type optionally specify a new type. Use this to
#'change between "polygon" and "minus.polygon" or to change
#'between "disc" and "minus.disc"
#'
#'@return If \code{new.type} is missing then the \code{winType} of win
#'is returned. Otherwise a new window is returned with \code{winType}
#'equal to \code{new.type}
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
        stop("cannot change from disc types to the specified new type")
      }
    } else if ( new.type == "polygon" || new.type == "minus.polygon" ) {

      if ( contains("polygon", winType(win)) )
      {
        attr(win, "winType") <- new.type
        return(win)
      }
      else
      {
        stop("cannot change from polygon types to the specified new type")
      }
    }
    else
    {
      stop("the specified new.type is invalid")
    }
  }
}


#' Assign new \code{\link{winType}} to a \code{\link{CMBWindow}}
#' @export
`winType<-` <- function(win, ..., value)
{
  return(winType(win, new.type = value))
}
