#' HPDataFrame class
#'
#' HPDataFrames are a type of \code{data.frame} designed to carry
#' data located on the unit sphere. Each row of a \code{HPDataFrame}
#' is associated with a HEALPix pixel index. The \code{HPDataFrame}
#' also holds an attribute called \code{nside} which stores the
#' HEALPix Nside parameter (i.e., the resolution of the HEALPix grid
#' that is being used).
#' Unlike \code{\link{CMBDataFrame}}, HPDataFrames may have
#' repeated pixel indices. They are made this way so that
#' multiple data points falling within a given pixel
#' can be stored in different rows of any given HPDataFrame.
#'
#' @param ... data, can be named vectors or a data.frame
#' @param nside integer, the nside parameter, i.e, resolution
#' @param ordering the HEALPix ordering scheme ("ring" or "nested")
#' @param auto.spix boolean. If TRUE then spix will be found from
#' the coordinates provided in the data. That is, each row of
#' data will be assigned the pixel index of its closest HEALPix
#' pixel center. There must be columns x,y,z for cartesian or
#' theta, phi for spherical colatitude and longitude respectively
#' @param spix a vector of HEALPix pixel indices indicating the
#' pixel locations of the data. Note that \code{spix} is ignored
#' if \code{auto.spix = TRUE}
#'
#' @examples
#' hpdf <- HPDataFrame(I = 1:12, nside = 1)
#' class(hpdf)
#' nside(hpdf)
#' ordering(hpdf)
#'
#' @export
HPDataFrame <- function(..., nside, ordering = "nested",
                        auto.spix = FALSE, spix)
{
  df <- data.frame(...)

  if ( missing(nside) ) {
    stop("nside must be specified")
  }

  if ( !auto.spix )
  {
    if ( missing(spix) )
    {
      pix <- 1:(12*nside^2)
    }
    else
    {
      pix <- spix
    }
  }
  else # auto.spix = TRUE. So, use nestSearch to determine pixel centers
  {

    if (all(c("x","y","z") %in% names(df)))
    {
      pix <- apply(df[,c("x","y","z")], MARGIN = 1, nestSearch,
                   nside = nside, index.only = TRUE)
    }
    else if (all(c("theta","phi") %in% names(df)))
    {
      df.xyz <- coords(df, new.coords = "cartesian")
      pix <- apply(df.xyz[,c("x","y","z")], MARGIN = 1, nestSearch,
                   nside = nside, index.only = TRUE)
    }
    else
    {
      stop(paste0("When auto.spix = TRUE there must be columns ",
                  "x, y, z (cartesian) or theta, phi (spherical)"))
    }

    if ( ordering != "nested" )
    {
      # Development stage: In future we will convert to desired ordering.
      warning(paste0("`auto.spix = TRUE` uses nestSearch, which requires ",
                     "nested ordering. The ordering was changed to nested."))
      ordering <- "nested"
    }
  }

  if (length(pix) != nrow(df))
  {
    stop(paste0("There should be a pixel index assigned to each row. ",
                "Perhaps use spix; perhaps specify nside correctly; ",
                "or perhaps use auto.spix = TRUE."))
  }

  attr(df, "pix") <- pix
  attr(df, "nside") <- nside
  attr(df, "ordering") <- ordering
  class(df) <- unique(c("HPDataFrame", class(df)))
  df
}

#'@export
`pix<-.HPDataFrame` <- function(hpdf,...,value) {
  attr(hpdf, "pix") <- value
  hpdf
}



#' HEALPix pixel indices from \code{\link{HPDataFrame}}
#'
#' If new.pix is unspecified then this function returns the vector of
#' HEALPix pixel indices from a HPDataFrame. If new.pix is specified then
#' this function returns a new HPDataFrame with the same number of rows
#' as \code{hpdf}, but with pix attribute \code{new.pix}. Thus,
#' \code{new.pix} must have length equal to \code{nrow(hpdf)}.
#'
#'
#'@param hpdf a \code{\link{HPDataFrame}}.
#'@param new.pix optional vector of pixel indices with
#'length equal to \code{nrow(hpdf)}
#'
#'@return
#' The vector of HEALPix pixel indices (integers) or,
#' if new.pix is specified,
#' a new HPDataFrame.
#'
#'@examples
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' pix(df)
#'
#'@export
pix.HPDataFrame <- function(hpdf, new.pix)
{
  if ( !missing(new.pix) )
  {
    if (nrow(hpdf) != length(new.pix))
    {
      stop("nrow(hpdf) not equal to length(new.pix)")
    }
    attr(hpdf, "pix") <- new.pix
    return(hpdf)
  }

  return(attr(hpdf, "pix"))
}


#' HEALPix Nside parameter from a \code{\link{HPDataFrame}}
#'
#' This function returns the HEALPix Nside parameter of a \code{\link{HPDataFrame}}
#'
#'@param hpdf a \code{\link{HPDataFrame}}.
#'
#'@return
#' The HEALPix Nside parameter
#'
#'@examples
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' nside(df)
#'
#'@export
nside.HPDataFrame <- function( hpdf )
{
  return( as.integer(attr( hpdf, "nside" )) )
}



#' Plot HPDataFrame
#'
#' This function produces a plot from a \code{\link{HPDataFrame}}.
#' If columns x,y,z (cartesian) or theta,phi (colatitude and longitude
#' respectively) are present in \code{hpdf}, then
#' these will be used as coordinates for plotting. Otherwise, the
#' HEALPix indices as in \code{pix(hpdf)} will be used. If HEALPix
#' indices are used and multiple rows correspond to a single pixel
#' index, then beware that values may be obfuscated in the plot,
#' and all locations are pixel centers.
#'
#'@param hpdf a HPDataFrame.
#'@param add if TRUE then this plot will be added to any existing plot.
#'Note that if \code{back.col} (see below) is specified then a new plot
#'window will be opened and \code{add = TRUE} will have no effect
#'@param sample.size optionally specifies the size of a simple random
#'sample to take before plotting. This can make the plot less
#'computationally intensive
#'@param type a single character indicating the type of item to plot.
#'Supported types are: 'p' for points, 's' for spheres, 'l' for lines,
#''h' for line segments from z = 0, and 'n' for nothing.
#'@param size the size of plotted points
#'@param box whether to draw a box
#'@param axes whether to draw axes
#'@param aspect either a logical indicating whether to adjust the
#'aspect ratio, or a new ratio.
#'@param col specify the colour(s) of the plotted points
#'@param back.col optionally specifies the background colour of
#'the plot. This argument is passed to rgl::bg3d.
#'@param labels optionally specify a vector of labels to plot,
#'such as words or vertex indices. If this is specified then
#'\code{rgl::text3d} is used instead of \code{rgl::plot3d}. Then
#'\code{length(labels)} must equal \code{nrow(hpdf)}
#'@param hp.boundaries integer. If greater than 0 then HEALPix
#'pixel boundaries at \code{nside = hp.boundaries} will be
#'added to the plot
#'@param hpb.col colour for the \code{hp.boundaries}
#'@param ... arguments passed to rgl::plot3d
#'
#'@return
#'A plot of the data locations
#'according to coordinate columns or HEALPix index
#'
#'@examples
#' hpdf <- HPDataFrame(I = rep(0,12), nside = 1)
#' plot(hpdf, size = 5, col = "yellow", back.col = "black",
#'      hp.boundaries = 1)
#'
#'@export
plot.HPDataFrame <- function(hpdf, intensities = "I",
                              add = FALSE, sample.size,
                              type = "p", size = 1, box = FALSE,
                              axes = FALSE, aspect = FALSE,
                              col = "blue", back.col = "black", labels,
                              hp.boundaries = 0, hpb.col = "gray", ...)
{
  pix <- pix(hpdf)

  if ( anyDuplicated(pix) > 0 ) {
    if ( !all(c("x","y","z") %in% names(hpdf))
         && !all(c("theta","phi") %in% names(hpdf)) )
    {
      warning(paste0("Some rows of hpdf share the same pixel index. ",
                   "Quantities may be obfuscated in final plot ",
                   "and sample.size may differ from actual size ",
                   "if pixel indices are relied on for plotting."))
    }
  }

  if ( !missing(sample.size) )
  {
    spix <- sample(pix, sample.size)
    hpdf <- hpdf[pix %in% spix,]
  }

  hpdf.xyz <- coords(hpdf, new.coords = "cartesian")

  ## Do the plotting
  if ( !add )
  {
    rgl::open3d()
    rgl::bg3d(back.col)
  }

  if ( missing(labels) )
  {
    rgl::plot3d(hpdf.xyz$x, hpdf.xyz$y, hpdf.xyz$z,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }
  else
  {
    rgl::text3d(hpdf.xyz$x, hpdf.xyz$y, hpdf.xyz$z, labels,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }

  if ( hp.boundaries > 0 )
  {
    plotHPBoundaries(nside = hp.boundaries, col = hpb.col)
  }
}




#' HEALPix ordering scheme from a HPDataFrame
#'
#' This function returns the HEALPix ordering scheme from a HPDataFrame.
#' The ordering scheme is either "ring" or "nested". If a new ordering
#' is specified, using e.g. \code{new.ordering = "ring"}, the
#' ordering scheme of the HPDataFrame will be converted.
#'
#'@param hpdf a \code{\link{HPDataFrame}}.
#'@param new.ordering specifies the new ordering ("ring" or "nest")
#'if a change of ordering scheme is desired.
#'
#'@return
#' The name of the HEALPix ordering scheme that is used in the
#' HPDataFrame hpdf, or a new hpdf with the desired new.ordering
#'
#'@examples
#' ## Plot using indices
#' df <- HPDataFrame(I = rep(0,12), nside = 1, ordering = "nested")
#' ordering(df)
#' ordering(df, new.ordering = "ring")
#'
#' ## Plot using coordinates
#' hp1 <- HPDataFrame(x = c(1,0,0),
#'                    y = c(0,1,0),
#'                    z = c(0,0,1),
#'                    nside = 1,
#'                    auto.spix = TRUE)
#' plot(hp, size = 5, hp.boundaries = 1)
#'
#'@export
ordering.HPDataFrame <- function( hpdf, new.ordering )
{
  if ( missing(new.ordering) )
  {
    return(attr( hpdf, "ordering" ))

  } else {
    new.ordering <- as.character(tolower(new.ordering))

    if (identical(as.character(attr(hpdf, "ordering")), new.ordering))
    {
      # Nothing to do

    } else if ( identical(new.ordering, "nested") ) {

      message("Converting to nested ordering...\n")
      pix(hpdf) <- rcosmo:::ring2nest(nside = nside(hpdf),
                                      pix = pix(hpdf))
      attr(hpdf, "ordering") <- "nested"

    } else if ( identical(new.ordering, "ring") ) {

      message("Converting to ring ordering...\n")
      pix(hpdf) <- nest2ring(nside = nside(hpdf), pix = pix(hpdf))
      attr(hpdf, "ordering") <- "ring"

    } else {

      stop("new.ordering must be either 'ring' or 'nested'")

    }

    return(hpdf)
  }
}

#' Assign new ordering scheme to HPDataFrame
#' @export
`ordering<-.HPDataFrame` <- function(hpdf,...,value) {
  rcosmo:::ordering(hpdf, new.ordering = value)
  hpdf
}





#' Coordinate system from a \code{\link{HPDataFrame}}
#'
#' Add or change coordinates in a \code{\link{HPDataFrame}}.
#' This does not affect the argument object \code{hpdf}.
#' Instead it returns a new \code{\link{HPDataFrame}}
#' with the desired coordinates. To change \code{hpdf}
#' directly see \code{\link{coords<-.HPDataFrame}}.
#'
#' If columns exist labelled x,y,z (cartesian) or theta, phi
#' (colatitude and longitude respectively), then these will be
#' treated as the coordinates of \code{hpdf} and converted
#' accordingly.
#' If columns x,y,z or theta,phi are not present then the healpix
#' pixel indices as given by \code{pix(hpdf)} are used for
#' assigning coordinates.
#'
#'@param hpdf a HPDataFrame.
#'@param new.coords specifies the new coordinate system
#'("spherical" or "cartesian")
#'@param healpix.only boolean. If TRUE then columns x,y,z
#'or theta, phi will be ignored and removed if present.
#'This forces the coordinates to be found from HEALPix
#'pixel indices only
#'
#'@return
#' A \code{\link{HPDataFrame}} with columns x,y,z (cartesian)
#' or theta, phi (colatitude and longitude respectively)
#'
#'@examples
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' coords(df, new.coords = "cartesian")
#' # Notice that df is unchanged
#' df
#'
#' # Instead, change df directly
#' coords(df) <- "spherical"
#'
#' ## specify cartesian coordinates then convert to spherical
#' hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
#'                    nside = 1, auto.spix = TRUE)
#' hp1 <- coords(hp1, new.coords = "spherical")
#'
#' ## Instead, ignore/drop existing coordinates and use HEALPix only
#' hp2 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
#'                    nside = 1, auto.spix = TRUE)
#' hp2 <- coords(hp1, new.coords = "spherical", healpix.only = TRUE)
#'
#'@export
coords.HPDataFrame <- function( hpdf, new.coords, healpix.only = FALSE )
{

  ns <- rcosmo:::nside(hpdf)
  od <- rcosmo:::ordering(hpdf)
  pix <- rcosmo:::pix(hpdf)
  nc <- (new.coords == "cartesian")

  if ( healpix.only == TRUE )
  {
    crd.cols <- which( names(hpdf) %in%
                  c("x","y","z","theta","phi") )
    if ( length(crd.cols) > 0 )
    {
      hpdf <- hpdf[ , -crd.cols]
    }

  }

  if ( new.coords == "spherical" )
  {
    if ( all(c("theta","phi") %in% names(hpdf)) )
    {
      return(df)
    }

    if ( !all(c("x","y","z") %in% names(hpdf)) )
    {
      # Use pix to convert to spherical
      sph <- rcosmo:::pix2coords_internal(nside = ns,
                nested = (od == "nested"),
                spix = pix,
                cartesian = FALSE)
      sph <- as.data.frame(sph)
      names(sph) <- c("theta","phi")
      hpdf <- cbind(hpdf, as.data.frame(sph))
    }
    else
    {
      x.i <- which(names(hpdf) == "x")
      y.i <- which(names(hpdf) == "y")
      z.i <- which(names(hpdf) == "z")

      crds <- hpdf[,c(x.i, y.i, z.i)]
      crds <- rcosmo:::car2sph(crds)
      other <- hpdf[,-c(x.i, y.i, z.i), drop = FALSE]
      hpdf <- cbind(crds, other)
    }

  } else if ( new.coords == "cartesian" ) {

    if ( all(c("x","y","z") %in% names(hpdf)) )
    {
      return(hpdf)
    }

    if ( !all(c("theta","phi") %in% names(hpdf)) )
    {
      # Use pix to convert to cartesian
      xyz <- rcosmo:::pix2coords_internal(nside = ns,
               nested = (od == "nested"),
               spix = pix,
               cartesian = TRUE)
      xyz <- as.data.frame(xyz)
      names(xyz) <- c("x","y","z")
      hpdf <- cbind(hpdf, as.data.frame(xyz))
    }
    else
    {
      theta.i <- which(names(hpdf) == "theta")
      phi.i <- which(names(hpdf) == "phi")

      crds <- hpdf[,c(theta.i, phi.i)]
      crds <- rcosmo:::sph2car(crds)
      other <- hpdf[,-c(theta.i, phi.i), drop = FALSE]
      hpdf <- cbind(crds, other)
    }
  }

  class(hpdf) <- unique(c("HPDataFrame", class(hpdf)))
  attr(hpdf, "ordering") <- od
  attr(hpdf, "nside") <- ns
  attr(hpdf, "pix") <- pix

  return(hpdf)
}



#' Assign new coordinate system to a \code{\link{HPDataFrame}}
#'
#'@keywords internal
#'
#' @seealso \code{\link{coords.HPDataFrame}}
#'
#' @examples
#'
#' ## Create df with no coords, then create df2 with cartesian coords
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' df
#' df2 <- coords(df, new.coords = "cartesian")
#' df2
#' df
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' df
#'
#' @export
`coords<-.HPDataFrame` <- function(hpdf,...,value) {
  return(coords(hpdf, new.coords = value))
}


#' Check if an object is of class \code{\link{HPDataFrame}}
#'
#' @param hpdf Any R object
#'
#' @return TRUE if \code{hpdf} is a HPDataFrame, otherwise FALSE
#'
#'@examples
#'
#' df <- CMBDataFrame(nside = 16)
#' is.HPDataFrame(df)
#'
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' is.HPDataFrame(df)
#'
#' @export
is.HPDataFrame <- function(hpdf)
{
  identical(as.numeric(sum(class(hpdf) == "HPDataFrame")), 1)
}



#' Print a \code{\link{HPDataFrame}}
#'
#' This function neatly prints the contents of a HPDataFrame.
#'
#'@param hpdf a HPDataFrame.
#'@param ... arguments passed to \code{\link{print.tbl_df}}
#'
#'@return
#'Prints contents of the HPDataFrame to the console.
#'
#'@examples
#' df <- HPDataFrame(I = rep(0,12), nside = 1, ordering = "nested")
#' print(df)
#' df
#'
#'@export
print.HPDataFrame <- function(hpdf,...)
{
  cat("A HPDataFrame\n")
  print(tibble::as.tibble(hpdf), ...)
}




#' Geodesic area covered by a \code{\link{HPDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{hpdf}
#'
#'@param hpdf a HPDataFrame
#'
#'@return the sum of the areas of all pixels (rows) in hpdf
#'
#'@examples
#'
#' ## At low resolution, a few data points can
#' ## occupy a large pixel area, e.g.:
#' hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
#'                    nside = 1, auto.spix = TRUE)
#' pix(hp1)
#' geoArea(hp1) # pi = 1/4*(surface area of unit sphere)
#' plot(hp1, size = 5, hp.boundaries = 1)
#'
#'@export
geoArea.HPDataFrame <- function(hpdf)
{
  nside <- rcosmo:::nside(hpdf)
  if ( !is.numeric(nside) ) stop("problem with hpdf nside attribute")
  return(pi/(3*nside^2)*length(unique(rcosmo:::pix((hpdf)))))
}







#' Get a sub window from a \code{\link{HPDataFrame}}
#'
#' This function returns a
#' HPDataFrame containing the data in \code{hpdf} restricted to the
#' CMBWindow \code{new.window}. If the HPDataFrame has columns x,y,z
#' or theta, phi then these will be used to determine locations
#' with priority over the HEALPix indices in \code{pix(hpdf)}
#' unless \code{healpix.only = TRUE} is given. Note that
#' if \code{healpix.only = TRUE} then columns x,y,z or theta, phi
#' will be discarded and replaced with pixel center locations.
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows.
#'
#'If the argument is a list of CMBWindows, then interiors of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned data.frame will
#'be the points of \code{df} in the the intersection of
#'\eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned data.frame
#'consists of the points of \code{df} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned data.frame
#'will be the points of \code{df} in \eqn{B} (resp. \eqn{A}).
#'
#'@param hpdf A HPDataFrame.
#'@param new.window A single \code{\link{CMBWindow}} object or a list of them.
#'@param intersect A boolean that determines
#'the behaviour when \code{win} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param healpix.only A boolean. If the HPDataFrame has columns x,y,z
#' or theta, phi then these will be used to determine locations
#' with priority over the HEALPix indices in \code{pix(hpdf)}
#' unless \code{healpix.only = TRUE} is given. Note that
#' if \code{healpix.only = TRUE} then columns x,y,z or theta, phi
#' will be discarded and replaced with pixel center locations.
#'
#'@return
#' A HPDataFrame containing the data in \code{hpdf} restricted to the
#' CMBWindow \code{new.window}
#'
#'@examples
#'
#'ns <- 16
#'hpdf <- HPDataFrame(nside = ns, I = 1:(12*ns^2))
#'plot(hpdf)
#'
#'
#'
#'@export
window.HPDataFrame <- function(hpdf, new.window, intersect = TRUE, healpix.only = FALSE)
{
  if ( healpix.only )
  {
    hpdf <- rcosmo:::coords(hpdf, new.coords = "cartesian", healpix.only = TRUE)
  }

  return(subWindow(hpdf, win = new.window, intersect = intersect))
}
