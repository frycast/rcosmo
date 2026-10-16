#' Coordinate conversion generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' coords functions \code{\link{coords.CMBDataFrame}}, \code{\link{coords.CMBWindow}},
#' \code{\link{coords.HPDataFrame}}, \code{\link{coords.data.frame}}
#'
#' @param x An object.
#' @param ... Unused arguments.
#'
#' @seealso
#' \code{\link{coords.CMBDataFrame}}
#' \code{\link{coords.CMBWindow}}
#' \code{\link{coords.HPDataFrame}}
#' \code{\link{coords.data.frame}}
#'
#'
#'@export
#'
#' @return A character string naming the coordinate system when it is queried, or an object of the input class with converted coordinates when a new system is supplied. See the individual methods for supported classes.
coords <- function(x, ...) UseMethod("coords", x)

#' Coordinate conversion generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' coords functions \code{\link{coords.CMBDataFrame}}, \code{\link{coords.CMBWindow}},
#' \code{\link{coords.HPDataFrame}}, \code{\link{coords.data.frame}}
#'
#' @keywords internal
#'
#'@export
#'
#' @return An object of the same class as \code{x}, with coordinates converted to the requested system.
#' @param x The object to update.
#' @param ... Extra arguments passed to the corresponding method.
#' @param value The replacement value; see the corresponding accessor method.
`coords<-` <- function(x, ..., value) UseMethod("coords<-", x)

#'geoArea generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' geoArea functions \code{\link{geoArea.CMBDataFrame}}, \code{\link{geoArea.HPDataFrame}},
#' \code{\link{geoArea.CMBWindow}}
#'
#'@param x An object.
#'
#'@seealso
#'\code{\link{geoArea.CMBDataFrame}}
#'\code{\link{geoArea.HPDataFrame}}
#'\code{\link{geoArea.CMBWindow}}
#'
#'@export
#'
#' @return A numeric scalar giving the area covered by the object on the unit sphere, in steradians.
geoArea <- function(x) UseMethod("geoArea", x)


#'pix generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' pix functions \code{\link{pix.CMBDataFrame}}, \code{\link{pix.HPDataFrame}}
#'
#'@param x An object.
#'@param ... Extra arguments.
#'
#'@seealso
#'\code{\link{pix.CMBDataFrame}}
#'\code{\link{pix.HPDataFrame}}
#'
#'
#'@export
#'
#' @return A numeric vector of HEALPix pixel indices when queried, or an object of the input class with updated pixel indices when new indices are supplied.
pix <- function(x, ...) UseMethod("pix", x)

#'pix generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' pix functions \code{\link{pix.CMBDataFrame}}, \code{\link{pix.HPDataFrame}}
#'
#'
#'@seealso
#'\code{\link{pix.CMBDataFrame}}
#'\code{\link{pix.HPDataFrame}}
#'
#'@keywords internal
#'
#'@export
#'
#' @return An object of the same class as \code{x}, with updated HEALPix pixel indices.
#' @param x The object to update.
#' @param ... Extra arguments passed to the corresponding method.
#' @param value The replacement value; see the corresponding accessor method.
`pix<-` <- function(x, ..., value) UseMethod("pix<-", x)

#'nside generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' nside functions \code{\link{nside.CMBDataFrame}}, \code{\link{nside.HPDataFrame}}
#'
#'@param x An object.
#'
#'@seealso
#'\code{\link{nside.CMBDataFrame}}
#'\code{\link{nside.HPDataFrame}}
#'
#'
#'@export
#'
#' @return A numeric scalar giving the HEALPix nside parameter. The corresponding full-sky grid contains \code{12*nside^2} pixels.
nside <- function(x) UseMethod("nside", x)

#'ordering generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' ordering functions \code{\link{ordering.CMBDataFrame}}, \code{\link{ordering.HPDataFrame}}
#'
#'@param x An object.
#'@param ... Extra arguments.
#'
#'@seealso
#'\code{\link{ordering.CMBDataFrame}}
#'\code{\link{ordering.HPDataFrame}}
#'
#'
#'@export
#'
#' @return A character string, \code{"nested"} or \code{"ring"}, when queried, or an object of the input class reordered to the requested scheme.
ordering <- function(x, ...) UseMethod("ordering", x)

#'ordering generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' ordering functions \code{\link{ordering.CMBDataFrame}}, \code{\link{ordering.HPDataFrame}}
#'
#'
#'@seealso
#'\code{\link{ordering.CMBDataFrame}}
#'\code{\link{ordering.HPDataFrame}}
#'
#'@keywords internal
#'
#'@export
#'
#' @return An object of the same class as \code{x}, with rows reordered to the requested HEALPix ordering.
#' @param x The object to update.
#' @param ... Extra arguments passed to the corresponding method.
#' @param value The replacement value; see the corresponding accessor method.
`ordering<-` <- function(x, ..., value) UseMethod("ordering<-", x)


#'window generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' window functions \code{\link{window.CMBDataFrame}}, \code{\link{window.HPDataFrame}},
#' \code{\link{window.data.frame}}
#'
#'@param x An object.
#'@param ... Extra arguments.
#'
#'@seealso
#'\code{\link{window.CMBDataFrame}}
#'\code{\link{window.HPDataFrame}}
#'\code{\link{window.data.frame}}
#'
#'
#'@export
#'
#' @return The stored window (a \code{CMBWindow}, a list of windows, or \code{NULL} for full sky) when queried, or an object of the input class restricted to the requested sky region. See the individual methods.
window <- function(x, ...) UseMethod("window", x)


#'window generic
#'
#' Detailed descriptions and and examples can be found in documentation for specific
#' window functions \code{\link{window.CMBDataFrame}}, \code{\link{window.HPDataFrame}},
#' \code{\link{window.data.frame}}
#'
#'@seealso
#'\code{\link{window.CMBDataFrame}}
#'\code{\link{window.HPDataFrame}}
#'\code{\link{window.data.frame}}
#'
#'@keywords internal
#'
#'@export
#'
#' @return An object of the same class as \code{x}, with observations restricted to the specified window.
#' @param x The object to update.
#' @param ... Extra arguments passed to the corresponding method.
#' @param value The replacement value; see the corresponding accessor method.
`window<-` <- function(x, ..., value) UseMethod("window<-", x)


