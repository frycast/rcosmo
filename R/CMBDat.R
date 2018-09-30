#' CMBDat class
#'
#' The function \code{CMBDat} creates objects of class \code{CMBDat}.
#' These are lists containing header information and other metadata as well
#' as an element called data, whose columns may include, for example, the
#' intensity (I), polarisation (Q, U), PMASK and TMASK. It also may contain an
#' \code{\link{mmap}} object that points to the CMB map data table in the FITS
#' file.
#'
#'@param filename The path to the fits file.
#'@param mmap A boolean indicating whether to use memory mapping.
#'@param spix The sample pixels (rows) to read from the FITS file
#'binary data table (optional)
#'@return A list containing header information and other metadata
#'as well as an element called \code{data} where:
#'If \code{mmap = FALSE} then a \code{data.frame} is
#'included, named \code{data}, whose columns may include, for
#'example, the intensity (I), polarisation (Q, U), PMASK and TMASK.
#'If \code{mmap = TRUE} then a \code{\link{mmap}} object is returned
#'that points to the CMB map data table in the FITS file.
#'
#'@examples
#'
#' ## Ensure you have a FITS file with the correct path
#' ## before running the example:
#' ## download a FITS file and use real data
#' # downloadCMBMap()
#' # cmbdat <- CMBDat("CMB_map_smica1024.fits", mmap = TRUE)
#' # class(cmbdat)
#' # str(cmbdat)
#'
#'## View metadata
#'# cmbdat$header1
#'# cmbdat$header2
#'# cmbdat$resoln
#'# cmbdat$method
#'# cmbdat$coordsys
#'# cmbdat$nside
#'# cmbdat$hdr
#'
#'@export
CMBDat <- function(filename, mmap = FALSE, spix) {

  # FITS standard: 2880byte blocks, 80char keyword strings
  chars <- 80L
  bytes <- 2880L

  zz <- file(filename, "rb")
  header1 <- FITSio::readFITSheader(zz)
  header2 <- FITSio::readFITSheader(zz)
  hdr <- FITSio::parseHdr(header2)

  blocks <- ceiling(length(header1) * chars / bytes) +
            ceiling(length(header2) * chars / bytes)

  ### Info needed for reading in binary data
  # Number of columns (e.g. I,Q,U,PMASK,TMASK)
  tfields <- as.numeric(hdr[which(hdr == "TFIELDS") + 1])

  ## Full map metadata
  # Resolution (arcmin)
  resoln <- as.numeric(hdr[which(hdr == "RESOLN") + 1])
  # The method e.g. SMICA
  method <- hdr[which(hdr == "METHOD") + 1]
  # Coordinate system e.g. Galactic
  coordsys <- tolower(hdr[which(hdr == "COORDSYS") + 1])
  # The Nside parameter
  nside <- as.integer(hdr[which(hdr == "NSIDE") + 1])
  # Value representing bad data
  baddata <- hdr[which(hdr == "BAD_DATA") + 1]
  # HEALPix ordering scheme
  ordering <- tolower(hdr[which(hdr == "ORDERING") + 1])

  # Column map metadata
  TFORMn <- character(tfields)
  TTYPEn <- character(tfields)
  TUNITn <- character(tfields)
  for (i in 1:tfields) {

    # FITS Format Code, e.g. E, B, ...
    tmp <- gsub(" ", "", hdr[which(hdr == paste("TFORM",
                                                i, sep = "")) + 1])
    TFORMn[i] <- tmp

    # column names
    tmp <- gsub(" ", "", hdr[which(hdr == paste("TTYPE",
                                                i, sep = "")) + 1])
    TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)

    # Column units e.g. K_CMB
    tmp <- gsub(" ", "", hdr[which(hdr == paste("TUNIT",
                                                i, sep = "")) + 1])
    TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
  }


  # Prepare FITS Format code details
  btype <- rep(NA, tfields)
  for (i in 1:tfields) {

    switch(tolower(TFORMn[i]),
           b = {
             btype[i] <- 1
             },
           e = {
             btype[i] <- 2
             },
           stop("Unknown TFORMn, contact rcosmo package developers"))
  }

  mystruct <- do.call(mmap::struct, CTypeExpression(TTYPEn, btype))
  map <- mmap::mmap(file = filename,
                    mode = mystruct,
                    off = blocks * bytes,
                    endian = "big")
  mmap::extractFUN(map) <- function(X) do.call(data.frame, X)



  if ( mmap == FALSE ) {

    if ( missing(spix) ) {

      spix <- 1:(12 * nside ^ 2)
    }

    col <- map[spix]
    mmap::munmap(map)

  } else {

    col <- map
  }

  close(zz)

  cmbdat <- list(data = col, colnames = TTYPEn, hdr = hdr,
                 resoln = resoln, method = method,
                 coordsys = coordsys, ordering = ordering,
                 nside = nside, baddata = baddata,
                 header1 = header1, header2 = header2,
                 mmap = mmap)
  class(cmbdat) <- c("CMBDat", "list")
  return(cmbdat)
}


# # # # # # # # # # # # # # # # # # # # # # # # # #
#       Helper functions for mmap C_types         #
CTypeSwitch <- function(x) {
  sapply(x, switch, quote(mmap::int8()), quote(mmap::real32()))
}

simplifyNames <- function(x) {

  sapply(x, function(x) {

      if (toupper(x) == "I_STOKES") return("I")
      if (toupper(x) == "Q_STOKES") return("Q")
      if (toupper(x) == "U_STOKES") return("U")
      return(x)
  })
}

CTypeExpression <- function(TTYPEn, btype) {

  names <- simplifyNames(TTYPEn)
  CTypes <- CTypeSwitch(btype)
  names(CTypes) <- names
  return(CTypes)
}
# # # # # # # # # # # # # # # # # # # # # # # # # #




#' Get a sub window from a \code{\link{CMBDat}} object
#'
#' This function returns a
#' data.frame containing the data in \code{x} restricted to the
#' CMBWindow \code{new.window}
#'
#'Windows that are tagged with \code{set.minus} (see \code{\link{CMBWindow}})
#'are treated differently from other windows.
#'
#'If the argument is a list of CMBWindows, then interious of all windows whose
#'winType does not include "minus" are united (let \eqn{A} be their union) and
#'exteriors of all windows whose winType does include "minus" are intersected,
#'(let \eqn{B} be their intersection). Then, provided that
#'\code{intersect = TRUE} (the default), the returned data.frame will
#'be the points of \code{cmbdat$data} in the the intersection of
#'\eqn{A} and \eqn{B}.
#'Otherwise, if \code{intersect = FALSE}, the returned data.frame
#'consists of the points of \code{x$data} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned data.frame
#'will be the points of \code{x} in \eqn{B} (resp. \eqn{A}).
#'
#'@param x a \code{\link{CMBDat}} object.
#'@param new.window A single \code{\link{CMBWindow}} object or a list of them.
#'@param intersect A boolean that determines
#'the behaviour when \code{new.window} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'@param ... Unused arguments.
#'
#'@return
#' A CMBDataFrame containing the data in \code{x} restricted to the
#' CMBWindow \code{new.window}
#'
#'@examples
#'
#'win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#'
#'## Ensure you have a FITS file with correct path
#'## before uncommenting and running the rest of the example:
#'# cmbdat <- CMBDat("CMB_map_smica1024.fits", mmap = TRUE)
#'# class(cmbdat)
#'# cmbdat.win <- window(cmbdat, new.window = win1)
#'# class(cmbdat.win)
#'
#'@export
window.CMBDat <- function(x, new.window, intersect = TRUE, ...) {

  return(subWindow(x, win = new.window, intersect = intersect))
}




#' Check if an object is of class CMBDat
#'
#' @param cmbdf Any R object
#'
#' @return TRUE if \code{cmbdf} is a CMBDat object, otherwise FALSE
#'
#'@examples
#' ## First download the map
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' # cmbdat <- CMBDat("CMB_map_smica1024.fits", mmap = TRUE)
#' # class(cmbdat)
#' # is.CMBDat(cmbdat)
#'
#' @export
is.CMBDat <- function(cmbdf)
{
  identical(as.numeric(sum(class(cmbdf) == "CMBDat")), 1)
}
