#' Summarise a \code{\link{CMBWindow}}
#'
#' This function produces a summary from a CMBWindow
#'
#'@param object A CMBWindow.
#'@param ... Unused arguments.
#'
#'@return
#'A summary includes window's type and area
#'
#'@examples
#'
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' summary(win)
#'
#' win1<- CMBWindow(x=0,y=3/5,z=4/5,r=0.8, set.minus = TRUE)
#' summary(win1)
#'
#'
#'@export
summary.CMBWindow <- function(object, ...)
{
  win <- object
  ans <- list()
  ans$area <- geoArea(win)
  ans$coords <- coords(win)
  ans$assumedConvex <- assumedConvex(win)
  ans$winType <- winType(win)

  if ( coords(win) == "cartesian" )
  {
    crds <- car2sph(win)
  }
  else
  {
    crds <- win[,c("theta","phi")]
  }

  ans$coords <- crds

  if ( winType(win) == "polygon" || winType(win) == "minus.polygon" )
  {
    ans$rangeTheta <- c(min(crds[,"theta"]), max(crds[,"theta"]))
    ans$rangePhi <- c(min(crds[,"phi"]), max(crds[,"phi"]))
  }

  if ( winType(win) == "disc" || winType(win) == "minus.disc" )
  {
    ans$center <- crds[, c("theta","phi")]
    ans$radius <- win[,"r"]
  }

  class(ans) <- "summary.CMBWindow"

  return(ans)
}



#'Print a summary of a \code{\link{CMBWindow}}
#'
#'@param x a \code{summary.CMBWindow} object, i.e.,
#'the output of \code{\link{summary.CMBWindow}}
#'
#'@keywords internal
#'
#'@export
print.summary.CMBWindow <- function(x, ...)
{
  cat(
    cli::boxx(c(paste0("Window type: ", x$winType),
                paste0("Window area: ", round(x$area,4)))
    ),
    sep = ""
  )
}




#' Check if an object is a CMBWindow
#'
#' @param win any object
#'
#' @return TRUE or FALSE depending if win is a CMBWindow
#'
#' @examples
#'
#' win <- CMBWindow(x=0,y=3/5,z=4/5,r=0.8, set.minus = TRUE)
#' is.CMBWindow(win)
#'
#'@export
is.CMBWindow <- function(win)
{
  return(identical(as.numeric(sum(class(win) == "CMBWindow")), 1))
}





#' Visualise a \code{\link{CMBWindow}}
#'
#'@param x A CMBWindow.
#'@param add if TRUE then this plot will be added to any existing plot.
#'Note that if \code{back.col} (see below) is specified then a new plot
#'window will be opened and \code{add = TRUE} will have no effect
#'@param type a single character indicating the type of item to plot.
#'Supported types are: 'p' for points, 's' for spheres, 'l' for lines,
#''h' for line segments from z = 0, and 'n' for nothing.
#'@param col specify the colour(s) of the plotted points
#'@param size the size of plotted points
#'@param box whether to draw a box
#'@param axes whether to draw axes
#'@param aspect either a logical indicating whether to adjust the
#'aspect ratio, or a new ratio.
#'@param back.col specifies the background colour of
#'the plot. This argument is passed to rgl::bg3d.
#'@param ... arguments passed to rgl::plot3d
#'
#'@examples
#'
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' win2 <- CMBWindow(theta = c(2*pi/3,3*pi/4,3*pi/4, 2*pi/3), phi = c(pi/4,pi/4,pi/3,pi/3))
#' plot(win1)
#' plot(win2)
#'
#'@export
plot.CMBWindow <- function(x, add = TRUE, type = "l",
                           col = "red",
                           size = 2, box = FALSE,
                           axes = FALSE, aspect = FALSE,
                           back.col, ...)
{
  win <- x

  if ( coords(win) == "spherical" )
  {
    rcosmo::coords(win) <- "cartesian"
  } else if ( rcosmo::coords(win) == "cartesian" ) {
    # do nothing
  } else {
   stop("'win' must be in either spherical or cartesian coordinates")
  }

  if ( contains("polygon", rcosmo::winType(win)) )
  {
    boundary <- polygonBoundary(win)
  }
  else
  {
    boundary <- discBoundary(win)
  }

  if ( !missing(back.col) )
  {
    rgl::open3d()
    rgl::bg3d(back.col)
  }
  rgl::plot3d( boundary, add = add, type = type, col = col, size = size,
               box = box, axes = axes, aspect = aspect, ... )

}

## HELPER FUNCTION FOR plot.CMBWindow
discBoundary <- function( disc, eps = 0.01 )
{
  rows <- nrow(disc)
  if (rows != 1) {
    stop("Disc must have just 1 center point and 1 radius.");
  }

  center <- as.numeric(disc[1,c("x","y","z")])

  phi <- seq(0,2*pi,by = 0.01)
  r <- as.numeric(disc[1,"r"])
  theta <- rep(r, length(phi))

  boundary <- sph2car(data.frame(theta = theta, phi = phi))
  boundary <- rodrigues(c(0,0,1), center, boundary)

  return(boundary)
}

## HELPER FUNCTION FOR plot.CMBWindow
polygonBoundary <- function( vertices.xyz, eps = 0.01 )
{
  rows <- nrow(vertices.xyz)
  if(rows < 3){
    stop("Polygon must have at least 3 vertices.")
  }

  boundary <- data.frame()
  for ( row in 1:rows )
  {
    #Geodesic is the shortest distance between two_points
    V1 <- vertices.xyz[row,]
    V2 <- vertices.xyz[1 + (row %% rows),]

    normal <- vector_cross(V1, V2)
    ## Rotate so that V1XV2 is moved to (0,0,1)
    rotated <- rodrigues(as.matrix(normal)[1,], c(0,0,1),
                                  rbind(V1,V2))

    ## This atan2 function returns negative thetas which are then corrected
    phi <- atan2(rotated[,2], rotated[,1])
    phi[phi < -1e-14] <- 2*pi + phi[phi < -1e-14]

    ## Take the shortest path
    phi <- sort(phi)
    if ( phi[2] - phi[1] > pi )
    {
      phi[1] <- phi[1] + 2*pi
      phi <- sort(phi)
    }
    line.phi <- seq(phi[1], phi[2], by = eps)

    ## Rotate so that all values are in [0,2*pi)
    line.phi[line.phi >= 2*pi] <- line.phi[line.phi >= 2*pi] - 2*pi

    line <- data.frame(phi = line.phi,
                       theta = rep(pi/2,length(line.phi)))
    line <- sph2car( line )

    line.rotated <-  as.data.frame(rodrigues(c(0,0,1),
                                                      as.matrix(normal)[1,],
                                                      line))
    names(line.rotated) <- c("x","y","z")
    boundary <- rbind(boundary, line.rotated)
  }
  boundary
}









#' Geodesic area of a \code{\link{CMBWindow}}
#'
#' @param x A CMBWindow.
#'
#' @return The spherical area inside the CMBWindow \code{x}.
#'
#'@examples
#'
#' ## A window that covers 1/8 of the unit sphere is constructed and its area is
#' ## pi/2 = 1/8*(surface area of unit sphere)
#'
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' geoArea(win)
#'
#'@export
geoArea.CMBWindow <- function(x)
{
  win <- x
  # Calculate the area of the spherical polygon
  win.xyz <- rcosmo::coords(win, new.coords = "cartesian")

  a <- switch(rcosmo::winType(win),
              polygon = polygonArea(win.xyz),
              minus.polygon = 4*pi - polygonArea(win.xyz),
              disc = 2*pi*(1 - cos(as.numeric(win$r))),
              minus.disc = 2*pi*(1 + cos(as.numeric(win$r))),
              stop("Could not determine window type using rcosmo::winType"))

  return(a)
}

## HELPER FUNCTION to calculate area of spherical polygon
# win must be a data.frame in cartesian coordinates.
# Works with counter-clockwise oriented polygons
polygonArea <- function(win)
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
    # swap these for clockwise oriented polygon if desired
    if ( det(matrix(c(A1, A2, A3), nrow = 3)) <= 0 )
    {
      angles[i] <- 2*pi - theta
    } else {
      angles[i] <- theta
    }
  }

  # Use Gauss-Bonnet Theorem for area of spherical polygon
  return(sum(angles) - (n-2)*pi)
}












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
#'@param x a CMBWindow, \code{win}.
#'@param new.coords specifies the new coordinate system
#'("spherical" or "cartesian")
#'if a change of coordinate system is desired
#'@param ... Unused arguments.
#'
#'@return
#' If new.coords is unspecified, then the name of the coordinate system
#' of \code{win} is returned. Otherwise a new CMBWindow is returned
#' equivalent to \code{win} but having the desired change of coordinates
#'
#'@examples
#' ## Create win with sperical coords, then change it to win1 with cartesian coords
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' coords(win)
#' win1 <- coords(win, new.coords = "cartesian")
#' coords(win1)
#'
#' ## Change back to spherical coordinates
#'
#' coords(win1) <- "spherical"
#' coords(win1)
#'
#'@export
coords.CMBWindow <- function( x, new.coords, ... )
{
  win <- x

  # If new.coords argument is missing then return the coordinate type
  if ( missing(new.coords) )
  {
    return(attr(win, "coords"))
  }
  else
  {
    new.coords <- as.character(tolower(new.coords))

    if ( contains("disc", winType(win)) )
    {
      r <- win[,"r"]
    }

    # Make sure that new.coords doesn't match current coords
    if ( attr(win, "coords") == new.coords )
    {
      # Nothing to do
    }
    else if ( new.coords == "spherical" )
    {
      # Convert to spherical
      win[,1:2] <- car2sph(win[,c("x", "y", "z")])

      # Add the radius if needed
      if ( contains("disc", winType(win))  )
      {
        win[,3] <- r
        names(win) <- c("theta", "phi", "r")
        win[,4] <- NULL
      }
      else
      {
        names(win) <- c("theta", "phi")
        win[,3] <- NULL
      }

      attr(win, "coords") <- "spherical"
    }
    else if ( new.coords == "cartesian" )
    {
      # Convert to cartesian
      win[,1:3] <- sph2car(win[,c("theta", "phi")])

      # Add the radius if needed
      if ( contains("disc", winType(win))  )
      {
        win[,4] <- r
        names(win) <- c("x", "y", "z", "r")
      }
      else
      {
        names(win) <- c("x", "y", "z")
      }

      attr(win, "coords") <- "cartesian"
    }

    return(win)
  }
}




#'Assign new coordinate system to CMBWindow
#'
#'@keywords internal
#'
#'@export
`coords<-.CMBWindow` <- function(x,...,value) {
  value <- tolower(value)
  if (rcosmo::coords(x) == value)
  {
    return(x)

  } else {

    return(rcosmo::coords(x, new.coords = value))
  }
}
