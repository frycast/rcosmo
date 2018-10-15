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
#' @param ... Data. Can be named vectors or a data.frame. May
#' include columns (x,y,z) or (theta, phi) representing
#' Cartesian or spherical coordinates
#' of points on the unit sphere.
#' @param nside Integer number \eqn{2^k}, the nside
#' parameter, i.e, resolution. If \code{nside} is unspecified, then
#' the an attempt is made to use columns x,y and z from the provided
#' data, as Cartesian coordinates, to calculate an nside that is
#' sufficient to ensure all points belong to unique pixels.
#' @param ordering The HEALPix ordering scheme ("ring" or "nested").
#' @param auto.spix Boolean. If TRUE then spix will be found from
#' the coordinates provided in the data. That is, each row of
#' data will be assigned the pixel index of its closest HEALPix
#' pixel center. There must be columns x,y,z for cartesian or
#' theta, phi for spherical colatitude and longitude respectively.
#' If \code{auto.spix = FALSE} then \code{nside} must be specified.
#' @param spix A vector of HEALPix pixel indices indicating the
#' pixel locations of the data. Note that \code{spix} is ignored
#' if \code{auto.spix = TRUE}.
#' @param assumedUniquePix A boolean. Sets the \code{assumedUniquePix}
#' attribute of the HPDataFrame. This attribute indicates whether
#' or not the rows of a HPDataFrame can be assumed to belong to
#' unique pixels.
#' @param delete.duplicates Boolean. If TRUE then rows
#' corresponding to duplicate pixel indices will be dropped
#' from the returned HPDataFrame, and assumedUniquePix will
#' be set to TRUE.
#'
#' @details
#' \code{HPDataFrame} with \code{auto.spix = TRUE} can be used to transform any
#' spherical data (not necessarily CMB) to the Healpix representation, see
#' Example 3 below.
#'
#' @examples
#'
#' ##Example 1.
#'
#' hp1 <- HPDataFrame(I=rnorm(5), nside = 1, spix = c(1,1,2,2,3))
#' pix(hp1)
#' coords(hp1, new.coords = "cartesian")
#' class(hp1)
#' assumedUniquePix(hp1)
#'
#' ##Example 2.
#'
#' # Where nside is not specified
#' sky <- CMBDataFrame(nside = 32, coords = "cartesian", ordering = "nested")
#' sky.s <- CMBDataFrame(sky, sample.size = 100)
#' hpdf <- HPDataFrame(sky.s, auto.spix = TRUE)
#' class(hpdf)
#' assumedUniquePix(hpdf)
#'
#' ##Example 3.
#' #
#' ## With earth data.
#' ## Download World Cities Database from
#' ## https://simplemaps.com/static/data/world-cities/basic/simplemaps_worldcities_basicv1.4.zip
#' ## unpack the file worldcities.csv
#' #
#' # worldcities <- read.csv("worldcities.csv")
#' # uscities <- worldcities[worldcities$country == "United States",]
#' #
#' # Prepare a data frame with cities' coordinates
#' # usdf <- data.frame(phi = pi/180*uscities$lng, theta = pi/2 - pi/180*uscities$lat,
#' #                  I=rep(1,length(uscities$lng)))
#' #
#' # Select k cities with different coordinates
#' # k <- 1000
#' # usdf <- usdf[sample(nrow(usdf), k), ]
#' # plot(usdf$phi, usdf$theta)
#' # usdf[duplicated(usdf), ]
#' # usdf<- usdf[!duplicated(usdf), ]
#' # usdf[duplicated(usdf), ]
#' # usdf <- coords(usdf, new.coords = "cartesian")
#' #
#' # Create and plot the corresponding HPDataFrame with unique pixels
#' # ushp <- HPDataFrame(usdf, auto.spix = TRUE)
#' # plot(ushp, size = 2)
#'
#' @export
HPDataFrame <- function(..., nside, ordering = "nested",
                        auto.spix = FALSE, spix,
                        assumedUniquePix = FALSE,
                        delete.duplicates = FALSE) {

  if ( !auto.spix ) {

    if ( missing(nside) ) {

      stop("If auto.spix = FALSE, then nside must be specified")
    }

    if ( missing(spix) ) {

      pix <- 1:(12*nside^2)
    } else {

      pix <- spix
    }

    args <- list(...)
    nargs <- length(args)
    if ( nargs == 0 ) {

      df <- data.frame(I = rep(NA, length(pix)))
    } else {

      df <- data.frame(args)
    }
  } else { # auto.spix = TRUE. So, use nestSearch to determine pixel centers

    df <- data.frame(...)

    if ( missing(nside) ) {

      nside <- separatingNside(df)

      if (is.infinite(nside)) {

        stop(paste0("The resolution required to separate pixels was infinite. ",
                    "Perhaps you need to remove duplicate rows or some locations ",
                    "are extremely close to eachother."))
      }

      assumedUniquePix <- TRUE

      if ( nside >  2^12) {

        ### FIX ME :
        # This is to prevent numeric overflow in the pix2coords_internal
        ##
        warning(paste0("The resolution required to seperate pixels is too large. ",
                    "Perhaps this is unnecessary and you can remove some locations ",
                    "that are very close to eachother. The resolution was capped ",
                    "at nside = 4096, so some indices may now be duplicates."))

        assumedUniquePix <- FALSE
        nside <- 2^12
      }


    }

    if (all(c("x","y","z") %in% names(df))) {

      cart <- TRUE

    } else if (all(c("theta","phi") %in% names(df))) {

      df <- coords(df, new.coords = "cartesian")
      cart <- FALSE

    } else {

      stop(paste0("When auto.spix = TRUE there must be columns ",
                  "x, y, z (cartesian) or theta, phi (spherical)"))
    }

    pix <- nestSearch(df[,c("x","y","z")], nside = nside,
                      index.only = TRUE)

    if ( !cart ) coords(df) <- "spherical"

    if ( ordering != "nested" ) {

      # Development stage: In future we will convert to desired ordering.
      warning(paste0("`auto.spix = TRUE` uses nestSearch, which requires ",
                     "nested ordering. The ordering was changed to nested."))
      ordering <- "nested"
    }
  }

  if ( length(pix) != nrow(df) ) {

    stop(paste0("There should be a pixel index assigned to each row. ",
                "Perhaps use spix; perhaps specify nside correctly; ",
                "or perhaps use auto.spix = TRUE."))
  }

  if ( delete.duplicates )
  {
    df <- df[!duplicated(pix),]
    assumedUniquePix <- TRUE
  }

  attr(df, "pix") <- pix
  attr(df, "nside") <- nside
  attr(df, "ordering") <- ordering
  attr(df, "healpixCentered") <- FALSE
  attr(df, "assumedUniquePix") <- assumedUniquePix
  class(df) <- c("HPDataFrame", "data.frame")
  df
}

#'@export
`pix<-.HPDataFrame` <- function(x,...,value) {
  attr(x, "pix") <- value
  x
}



#' HEALPix pixel indices from \code{\link{HPDataFrame}}
#'
#' If new.pix is unspecified then this function returns the vector of
#' HEALPix pixel indices from a HPDataFrame. If new.pix is specified then
#' this function returns a new HPDataFrame with the same number of rows
#' as \code{x}, but with pix attribute \code{new.pix}. Thus,
#' \code{new.pix} must have length equal to \code{nrow(x)}.
#'
#'
#'@param x a \code{\link{HPDataFrame}}.
#'@param new.pix optional vector of pixel indices with
#'length equal to \code{nrow(x)}
#'@param ... Unused arguments.
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
pix.HPDataFrame <- function(x, new.pix, ...)
{
  if ( !missing(new.pix) ) {

    if (nrow(x) != length(new.pix)) {

      stop("nrow(hpdf) not equal to length(new.pix)")
    }
    attr(x, "pix") <- new.pix
    return(x)
  }

  return(attr(x, "pix"))
}


#' HEALPix Nside parameter from a \code{\link{HPDataFrame}}
#'
#' This function returns the HEALPix Nside parameter
#' of a \code{\link{HPDataFrame}}
#'
#'@param x A \code{\link{HPDataFrame}}.
#'
#'@return
#' The HEALPix Nside parameter.
#'
#'@examples
#' df <- HPDataFrame(I = rep(0,12), nside = 1)
#' nside(df)
#'
#'@export
nside.HPDataFrame <- function( x ) {

  return( as.integer(attr( x, "nside" )) )
}



#' Plot HPDataFrame
#'
#' This function produces a plot from a \code{\link{HPDataFrame}}.
#' If columns x,y,z (cartesian) or theta,phi (colatitude and longitude
#' respectively) are present in \code{x}, then
#' these will be used as coordinates for plotting. Otherwise, the
#' HEALPix indices as in \code{pix(x)} will be used. If HEALPix
#' indices are used and multiple rows correspond to a single pixel
#' index, then beware that values may be obfuscated in the plot,
#' and all locations are pixel centers.
#'
#'@param x A HPDataFrame.
#'@param intensities The column name for the data in \code{x}
#'that is to be treated
#'as intensities for plotting.
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
#'\code{length(labels)} must equal \code{nrow(x)}
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
plot.HPDataFrame <- function(x, intensities = "I",
                              add = FALSE, sample.size,
                              type = "p", size = 1, box = FALSE,
                              axes = FALSE, aspect = FALSE,
                              col = "blue", back.col = "black", labels,
                              hp.boundaries = 0, hpb.col = "gray", ...) {

  hpdf <- x
  pix <- pix(hpdf)

  if ( anyDuplicated(pix) > 0 ) {

    if ( !all(c("x","y","z") %in% names(hpdf))
         && !all(c("theta","phi") %in% names(hpdf)) ) {

      warning(paste0("Some rows of hpdf share the same pixel index. ",
                   "Quantities may be obfuscated in final plot ",
                   "and sample.size may differ from actual size ",
                   "if pixel indices are relied on for plotting."))
    }
  }

  if ( !missing(sample.size) ) {

    spix <- sample(pix, sample.size)
    hpdf <- hpdf[pix %in% spix,]
  }

  hpdf.xyz <- coords(hpdf, new.coords = "cartesian")

  ## Do the plotting
  if ( !add ) {

    rgl::open3d()
    rgl::bg3d(back.col)
  }

  if ( missing(labels) ) {

    rgl::plot3d(hpdf.xyz$x, hpdf.xyz$y, hpdf.xyz$z,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  } else {

    rgl::text3d(hpdf.xyz$x, hpdf.xyz$y, hpdf.xyz$z, labels,
                col = col, type = type, size = size,
                box = box, axes = axes, add = add, aspect = aspect, ...)
  }

  if ( hp.boundaries > 0 ) {

    rcosmo::displayPixelBoundaries(nside = hp.boundaries, col = hpb.col)
  }
}




#' HEALPix ordering scheme from a HPDataFrame
#'
#' This function returns the HEALPix ordering scheme from a HPDataFrame.
#' The ordering scheme is either "ring" or "nested". If a new ordering
#' is specified, using e.g. \code{new.ordering = "ring"}, the
#' ordering scheme of the HPDataFrame will be converted.
#'
#'@param x a \code{\link{HPDataFrame}}.
#'@param new.ordering Specifies the new ordering ("ring" or "nest")
#'if a change of ordering scheme is desired.
#'@param ... Unused arguments.
#'
#'@return
#' The name of the HEALPix ordering scheme that is used in the
#' HPDataFrame x, or a new HPDataFrame
#' with the desired new.ordering
#'
#'@examples
#'
#' df <- HPDataFrame(I = rep(0,12), nside = 1, ordering = "nested")
#' ordering(df)
#' df1 <- ordering(df, new.ordering = "ring")
#' ordering(df1)
#'
#'@export
ordering.HPDataFrame <- function( x, new.ordering, ... ) {

  hpdf <- x

  if ( missing(new.ordering) ) {

    return(attr( hpdf, "ordering" ))

  } else {
    new.ordering <- as.character(tolower(new.ordering))

    if (identical(as.character(attr(hpdf, "ordering")), new.ordering))
    {
      # Nothing to do

    } else if ( identical(new.ordering, "nested") ) {

      message("Converting to nested ordering...\n")
      pix(hpdf) <- rcosmo::ring2nest(nside = nside(hpdf),
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
#' @keywords internal
#' @export
`ordering<-.HPDataFrame` <- function(x,...,value) {
  rcosmo::ordering(x, new.ordering = value)
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
#'@param x a HPDataFrame, \code{hpdf}.
#'@param new.coords specifies the new coordinate system
#'("spherical" or "cartesian")
#'@param healpixCentered boolean. If TRUE then columns x,y,z
#'or theta, phi will be ignored and removed if present.
#'This forces the coordinates to be found from HEALPix
#'pixel indices only. Then the HEALPixCentered
#'attribute of \code{hpdf} will be set to \code{TRUE}.
#'@param ... Unused arguments.
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
#' hp2 <- coords(hp1, new.coords = "spherical", healpixCentered = TRUE)
#'
#'@export
coords.HPDataFrame <- function( x, new.coords, healpixCentered = FALSE, ... ) {

  hpdf <- x

  ns <- rcosmo::nside(hpdf)
  od <- rcosmo::ordering(hpdf)
  pix <- rcosmo::pix(hpdf)

  if ( healpixCentered == TRUE ) {

    # Delete any coordinates so they will be created again
    # from HEALPix only
    crd.cols <- which( names(hpdf) %in%
                  c("x","y","z","theta","phi") )
    if ( length(crd.cols) > 0 ) {

      hpdf <- hpdf[ , -crd.cols, drop = FALSE]
    }

  }

  if ( new.coords == "spherical" ) {

    if ( all(c("theta","phi") %in% names(hpdf)) ) {

      return(hpdf)
    }

    if ( !all(c("x","y","z") %in% names(hpdf)) ) {

      # Use pix to convert to spherical
      sph <- pix2coords_internal(nside = ns,
                nested = (od == "nested"),
                spix = pix,
                cartesian = FALSE)
      sph <- as.data.frame(sph)
      names(sph) <- c("theta","phi")
      hpdf <- cbind(hpdf, as.data.frame(sph))


    } else {

      x.i <- which(names(hpdf) == "x")
      y.i <- which(names(hpdf) == "y")
      z.i <- which(names(hpdf) == "z")

      crds <- hpdf[,c(x.i, y.i, z.i)]
      crds <- car2sph(crds)
      other <- hpdf[,-c(x.i, y.i, z.i), drop = FALSE]
      hpdf <- cbind(crds, other)

    }

  } else if ( new.coords == "cartesian" ) {

    if ( all(c("x","y","z") %in% names(hpdf)) ) {

      return(hpdf)
    }

    if ( !all(c("theta","phi") %in% names(hpdf)) ) {

      # Use pix to convert to cartesian
      xyz <- pix2coords_internal(nside = ns,
               nested = (od == "nested"),
               spix = pix,
               cartesian = TRUE)
      xyz <- as.data.frame(xyz)
      names(xyz) <- c("x","y","z")
      hpdf <- cbind(hpdf, as.data.frame(xyz))

    } else {

      theta.i <- which(names(hpdf) == "theta")
      phi.i <- which(names(hpdf) == "phi")

      crds <- hpdf[,c(theta.i, phi.i)]
      crds <- sph2car(crds)
      other <- hpdf[,-c(theta.i, phi.i), drop = FALSE]
      hpdf <- cbind(crds, other)
    }
  }

  class(hpdf) <- unique(c("HPDataFrame", class(hpdf)))
  attr(hpdf, "ordering") <- od
  attr(hpdf, "nside") <- ns
  attr(hpdf, "pix") <- pix
  attr(hpdf, "healpixCentered") <- healpixCentered

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
`coords<-.HPDataFrame` <- function(x,...,value) {
  return(coords(x, new.coords = value))
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
is.HPDataFrame <- function(hpdf) {

  identical(as.numeric(sum(class(hpdf) == "HPDataFrame")), 1)
}



#' Print a \code{\link{HPDataFrame}}
#'
#' This function neatly prints the contents of a HPDataFrame.
#'
#'@param x A HPDataFrame.
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
print.HPDataFrame <- function(x,...) {

  cat("A HPDataFrame\n")
  print(tibble::as.tibble(x), ...)
}




#' Geodesic area covered by a \code{\link{HPDataFrame}}
#'
#' Gives the surface on the unit sphere
#' that is encompassed by all pixels in \code{x}.
#'
#'@param x A HPDataFrame.
#'
#'@return The sum of the areas of all pixels (rows) in \code{x}.
#'
#'@examples
#'
#' ## Generate random I for HPDataFrame
#' hp1 <- HPDataFrame(I=rnorm(5), nside = 1, spix = c(1,1,2,2,3))
#' pix(hp1)
#'
#' ## The total number of Healpix points at nside=1 equals 12. As hp1 has five
#' ## I values at 3 Healpix points, then the occupied area is
#' ## pi = 1/4*(surface area of unit sphere)
#'
#' geoArea(hp1)
#' plot(hp1, size = 5, hp.boundaries = 1)
#'
#'@export
geoArea.HPDataFrame <- function(x) {

  nside <- rcosmo::nside(x)
  if ( !is.numeric(nside) ) stop("problem with hpdf nside attribute")
  return(pi/(3*nside^2)*length(unique(rcosmo::pix((x)))))
}




#' Get a sub window from a \code{\link{HPDataFrame}}
#'
#' This function returns a
#' HPDataFrame containing the data in \code{hpdf} restricted to the
#' CMBWindow \code{new.window}. If the HPDataFrame has columns x,y,z
#' or theta, phi then these will be used to determine locations
#' with priority over the HEALPix indices in \code{pix(hpdf)}
#' unless \code{healpixCentered = TRUE} is given. Note that
#' if \code{healpixCentered = TRUE} then columns x,y,z or theta, phi
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
#'@param x A \code{\link{HPDataFrame}}.
#'@param new.window Optional.
#'A single \code{\link{CMBWindow}} object or a list of them.
#'@param intersect A boolean that determines
#'the behaviour when \code{new.window} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param healpixCentered A boolean. If the HPDataFrame has columns x,y,z
#' or theta, phi then these will be used to determine locations
#' with priority over the HEALPix indices in \code{pix(x)}
#' unless \code{healpixCentered = TRUE} is given. Note that
#' if \code{healpixCentered = TRUE} then columns x,y,z or theta, phi
#' will be discarded and replaced with pixel center locations.
#' @param ... Unused arguments.
#'
#'@return
#' A HPDataFrame containing the data in \code{x} restricted to the
#' CMBWindow \code{new.window}. Or, if \code{new.window} is
#' unspecified, then the window attribute of \code{x}
#' is returned instead (and may be NULL).
#'
#'@examples
#' ns <- 16
#' hpdf <- HPDataFrame(nside = ns, I = 1:(12*ns^2))
#' hpdf
#'
#' win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' plot(hpdf); plot(win1)
#'
#' hpdf.win <- window(hpdf, new.window = win1)
#' plot(hpdf.win, col = "yellow", size = 4, add = TRUE)
#' attributes(hpdf.win)
#' window(hpdf.win)
#' hpdf.win
#'
#'
#'@export
window.HPDataFrame <- function(x, new.window, intersect = TRUE,
                               healpixCentered = FALSE, ...) {

  if ( missing(new.window) ) {

    return(attr(x, "window"))
  }

  if ( healpixCentered ) {

    hpdf <- rcosmo::coords(x, new.coords = "cartesian",
                            healpixCentered = TRUE)
  }

  return(subWindow(x, win = new.window, intersect = intersect))
}



#' Summarise a \code{\link{HPDataFrame}}
#'
#' This function produces a summary from a HPDataFrame.
#'
#'@param object A HPDataFrame.
#'@param intensities the name of a column specifying intensities
#'(or potentially another numeric quantity of interest)
#'@param ... Unused arguments.
#'
#'@return
#'A summary includes window's type and area,
#' total area covered by observations,
#' and main statistics for intensity values
#'
#' @examples
#' ns <- 2
#' hpdf <- HPDataFrame(I = rnorm(12*ns^2), nside = 2,
#'                     ordering = "nested")
#' win <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#' hpdf.win <- window(hpdf, new.window = win)
#' summary(hpdf.win)
#'
#'@export
summary.HPDataFrame <- function(object, intensities = "I", ...) {

  hpdf <- object

  ans <- list(intensities = summary(hpdf[,intensities, drop = TRUE]))

  if ( is.null(window(hpdf)) ) {

    ans$window <- "full sky"
  } else {

    ans$window <- window(hpdf)
  }

  ans$names <- names(hpdf)
  ans$ordering <- ordering(hpdf)
  ans$HEALPixCentered <- attr(hpdf, "HEALPixCentered")
  ans$nside <- nside(hpdf)
  ans$pix <- pix(hpdf)
  ans$n <- nrow(hpdf)
  ans$area <- geoArea(hpdf)
  class(ans) <- "summary.HPDataFrame"
  return(ans)
}




#'Print a summary of a HPDataFrame
#'
#'@keywords internal
#'
#'@param x a \code{summary.HPDataFrame} object, i.e.,
#'the output of \code{\link{summary.HPDataFrame}}
#'
#'@export
#'
#'
print.summary.HPDataFrame <- function(x, ...) {
  cat(
    cli::rule(center = " HPDataFrame Object ", line = "bar4"), "\n",
    sep = ""
  )

  # Window loop details here in boxes
  if ( any(x$window != "full sky") ) {

    cat("Number of CMBWindows: ", length(x$window), "\n" )
    if ( length(x$window) <= 5 ) {

      lapply(x$window, function(x) { print(summary(x)); cat("\n\n") } )
    } else {

      cat("Too many windows to print them all here", "\n\n")
    }
  } else {

    cat("Full sky map\n")
  }
  cat("HEALPix Centered: ", x$HEALPixCentered, "\n")

  cat("Total area covered by all pixels: ", x$area, "\n")


  # Summary of intensities
  cat(cli::rule(line = "~"), "\n", sep = "")
  cli::cat_line("Intensity quartiles")
  print(x$intensities)

  # Finishing line
  cat(
    cli::rule(line = "="),
    sep = ""
  )
}



#'separatingNside
#'
#'@param df A data.frame. Must have columns x,y,z
#'for Cartesian coordinates that must represent points
#'on the unit sphere.
#'
#'@return An nside (power of 2) such that all points
#'in \code{df} will belong to unique pixels.
#'
#'@examples
#' sky <- CMBDataFrame(nside = 16, coords = "cartesian", ordering = "nested")
#' separatingNside(sky)
#'
#'@keywords internal
#'
#'@export
separatingNside <- function(df) {

  dist <- minDist(df)
  return(minDist2nside(dist))
}


#'Check if object was assumed to have unique HEALPix indices
#'
#'The function checks object's attribute assumedUniquePix. The attribute is True if
#'the  object was assumed to have rows that correspond to unique
#'HEALPix pixel indices.
#'
#'@param obj Any object
#'
#'@return A boolean. This is TRUE if \code{obj}
#'is a \code{\link{CMBDataFrame}} or a \code{\link{HPDataFrame}} whose
#'rows were assumed to correspond to unique HEALPix pixel indices.
#'
#'@examples
#'
#' hp1 <- HPDataFrame(I=rnorm(5), nside = 1, spix = c(1,1,2,2,3))
#' pix(hp1)
#' coords(hp1, new.coords = "cartesian")
#' assumedUniquePix(hp1)
#'
#' sky <- CMBDataFrame(nside = 32, coords = "cartesian", ordering = "nested")
#' sky.s <- CMBDataFrame(sky, sample.size = 100)
#' hpdf <- HPDataFrame(sky.s, auto.spix = TRUE)
#' assumedUniquePix(hpdf)
#'
#'@export
assumedUniquePix <- function(obj) {

  if (is.CMBDataFrame(obj)) return(TRUE)
  if (is.HPDataFrame(obj)) return(attr(obj, "assumedUniquePix"))
  return(FALSE)
}





