#' Handling and Analysing CMB data
#'
#'The package \code{\link{rcosmo}} offers various tools for
#'\itemize{
#' \item{Downloading and transforming Cosmic Microwave Background
#' radiation (CMB) and spherical data}
#' \item{Working with Hierarchical Equal Area
#' isoLatitude Pixelation of a sphere (Healpix)}
#' \item{Spherical geometry}
#' \item{Statistical analysis of CMB and spherical data}
#' \item{Visualisation of Healpix data}
#' }
#'
#' Most of \code{\link{rcosmo}} features were developed for CMB,
#' but it can also be used for other spherical data. It contains
#' tools for transforming spherical data in cartesian and geographic
#' coordinates to the Healpix representation.
#'
#'
#'@section Update: Current updates are available through
#' URL: https://github.com/VidaliLama/rcosmo
#'
#'@section  BugReports: https://github.com/VidaliLama/rcosmo/issues
#'
#'@description Handling and Analysing Spherical, Healpix and Cosmic Microwave Background
#' data on a HEALPix grid.
#'
#'@author Daniel Fryer \email{d.fryer@latrobe.edu.au},
#'Andriy Olenko \email{a.olenko@latrobe.edu.au}, Ming Li \email{Ming.Li@latrobe.edu.au},
#'Yuguang Wang.
#'@docType package
#'@name rcosmo
#'@importFrom Rcpp evalCpp
#'@useDynLib rcosmo, .registration = TRUE
#'@aliases rcosmo-package
NULL

