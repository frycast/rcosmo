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
`window<-` <- function(x, ..., value) UseMethod("window<-", x)


