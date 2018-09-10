#' Coordinate conversion generic
#'
#' @param x An object
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
#' @keywords internal
#'
#'@export
`coords<-` <- function(x, ..., value) UseMethod("coords<-", x)

#'geoArea generic
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

#'maxDist generic
#'
#'@param x An object.
#'
#'@seealso
#'\code{\link{maxDist.CMBDataFrame}}
#'\code{\link{maxDist.CMBWindow}}
#'
#'@export
maxDist <- function(x) UseMethod("maxDist", x)

#'pix generic
#'
#'@param x An object.
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
#'@param x An object.
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
#'@param x An object.
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
#'@seealso
#'\code{\link{window.CMBDataFrame}}
#'\code{\link{window.HPDataFrame}}
#'\code{\link{window.data.frame}}
#'
#'@keywords internal
#'
#'@export
`window<-` <- function(x, ..., value) UseMethod("window<-", x)
