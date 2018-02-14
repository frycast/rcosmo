
#' visualise a \code{\link{CMBWindow}}
#'
#'@param win a CMBWindow
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
#'@param eps the geodesic distance between consecutive points to draw
#'on the window boundary
#'@param ... arguments passed to rgl::plot3d
#'
#'@export
plot.CMBWindow <- function(win, add = TRUE, type = "l",
                           col = "red",
                           size = 2, box = FALSE,
                           axes = FALSE, aspect = FALSE,
                           back.col, ...)
{
  if ( coords(win) == "spherical" )
  {
    coords(win) <- "cartesian"
  } else if ( coords(win) == "cartesian" ) {
    # do nothing
  } else {
   stop("'win' must be in either spherical or cartesian coordinates")
  }

  boundary <- polygonBoundary( win, 0.01 )

  if ( !missing(back.col) )
  {
    rgl::open3d()
    rgl::bg3d(back.col)
  }
  rgl::plot3d( boundary, add = add, type = type, col = col, size = size,
               box = box, axes = axes, aspect = aspect, ... )

}


## HELPER FUNCTION FOR plot.CMBWindow
polygonBoundary <- function( vertices.xyz, eps = 0.001 )
{
  rows <- nrow(vertices.xyz)
  if(rows < 3){
    stop("Polygon must have at least 3 vertices.");
  }
  #vertices_xyz_cycled <- rbind(vertices_xyz[2:nrow(vertices_xyz),], vertices_xyz[1,])

  boundary <- data.frame()
  for ( row in 1:rows )
  {
    #Geodesic is the shortest distance between two_points
    V1 <- vertices.xyz[row,]
    V2 <- vertices.xyz[1 + (row %% rows),]

    normal <- vector_cross(V1, V2)
    rotated <- rodrigues(as.matrix(normal)[1,], c(0,0,1), rbind(V1,V2))

    two.longitudes <- atan2(rotated[,2], rotated[,1]) # latitudes are both pi/2

    line.longitudes <- seq(two.longitudes[1], two.longitudes[2], by = eps)
    line <- data.frame(phi = line.longitudes, theta = rep(pi/2,length(line.longitudes)))
    line <- sph2car( line )

    line.rotated <-  as.data.frame(rodrigues(c(0,0,1), as.matrix(normal)[1,], line))
    names(line.rotated) <- c("x","y","z")
    boundary <- rbind(boundary, line.rotated)
  }
  boundary
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
