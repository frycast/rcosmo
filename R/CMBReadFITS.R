






############################################################################
## WARNING: THIS FUNCTION IS INCOMPLETE: THE mmap PARAMETER REQUIRES WORK ##
############################################################################
#' Read CMB data from a FITS file.
#'
#' \code{CMBReadFITS} is adapted from the \code{\link{readFITS}}
#' function in package
#'   \href{https://cran.r-project.org/web/packages/FITSio/index.html}{FITSio}.
#'   \code{CMBReadFITS} is in development stage and will only work with
#'   'CMB_map_smica1024.fits'.
#'   When it works, \code{CMBReadFITS} is much faster than
#'   \code{\link{readFITS}}.
#'   However, \code{\link{readFITS}} is more general and so is more
#'   likely to work.
#'
#'
#' The function \code{CMBReadFITS} creates objects of class \code{CMBDat}.
#' These are lists containing header information and other metadata as well
#' as an element called data, whose columns may include, for example, the
#' intensity (I), polarisation (Q, U), PMASK and TMASK. It also may contain an
#' \code{\link{mmap}} object that points to the CMB map data table in the FITS
#' file.
#'
#'@aliases CMBDat
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
#' cmbdat <- CMBReadFITS("CMB_map_smica1024.fits", mmap = TRUE)
#' class(cmbdat)
#' str(cmbdat)
#'
#'# View metadata
#'dat$header1
#'dat$header2
#'dat$resoln
#'dat$method
#'dat$coordsys
#'dat$nside
#'dat$hdr
#'
#'@export
CMBReadFITS <- function(filename, mmap = FALSE, spix) {

  # FITS standard: 2880byte blocks, 80char keyword strings
  chars <- 80L
  bytes <- 2880L

  zz <- file(filename, "rb")
  header1 <- FITSio::readFITSheader(zz)
  header2 <- FITSio::readFITSheader(zz)
  hdr <- FITSio::parseHdr(header2)

  blocks <- ceiling(length(header1)*chars/bytes) +
            ceiling(length(header2)*chars/bytes)

  # Info needed for reading in binary data
  naxis1 <- as.numeric(hdr[which(hdr == "NAXIS1") + 1]) # Number of bytes per row
  naxis2 <- as.numeric(hdr[which(hdr == "NAXIS2") + 1]) # Number of rows
  tfields <- as.numeric(hdr[which(hdr == "TFIELDS") + 1]) # Number of columns (e.g. I,Q,U,PMASK,TMASK)

  # Full map metadata
  resoln <- as.numeric(hdr[which(hdr == "RESOLN") + 1]) # Resolution (arcmin)
  method <- hdr[which(hdr == "METHOD") + 1] # The method e.g. SMICA
  coordsys <- tolower(hdr[which(hdr == "COORDSYS") + 1]) # Coordinate system e.g. Galactic
  nside <- as.integer(hdr[which(hdr == "NSIDE") + 1]) # The Nside parameter
  baddata <- hdr[which(hdr == "BAD_DATA") + 1] # Value representing bad data
  ordering <- tolower(hdr[which(hdr == "ORDERING") + 1]) # HEALPix ordering scheme

  # Column map metadata
  TFORMn <- character(tfields)
  TTYPEn <- character(tfields)
  TUNITn <- character(tfields)
  for (i in 1:tfields) {

    tmp <- gsub(" ", "", hdr[which(hdr == paste("TFORM", # FITS Format Code, e.g. E, B, ...
                                                i, sep = "")) + 1])
    TFORMn[i] <- tmp

    tmp <- gsub(" ", "", hdr[which(hdr == paste("TTYPE", # column names
                                                i, sep = "")) + 1])
    TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)

    tmp <- gsub(" ", "", hdr[which(hdr == paste("TUNIT", # Column units e.g. K_CMB
                                                i, sep = "")) + 1])
    TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
  }


  # Prepare FITS Format code details
  bsize <- rep(NA, tfields)
  btype <- rep(NA, tfields)
  bsign <- rep(NA, tfields)
  for (i in 1:tfields)
  {
    switch(tolower(TFORMn[i]),
           b = {
             btype[i] <- 1},
           e = {
             btype[i] <- 2},
           stop("Unknown TFORMn, contact rcosmo package developers"))
  }

  swap <- "big" != .Platform$endian

  mystruct <- do.call(mmap::struct, CTypeExpression(TTYPEn, btype))
  map <- mmap::mmap(file = filename,
                    mode = mystruct,
                    off = blocks*bytes,
                    endian = "big")
  mmap::extractFUN(map) <- function(X) do.call(data.frame, X)

  if ( mmap == FALSE )
  {
    if ( missing(spix) )
    {
      spix <- 1:(12*nside^2)
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
      return(x)})
}

CTypeExpression <- function(TTYPEn, btype)
{
  names <- simplifyNames(TTYPEn)
  CTypes <- CTypeSwitch(btype)
  names(CTypes) <- names
  return(CTypes)
}
# # # # # # # # # # # # # # # # # # # # # # # # # #








#' Get a sub window from a \code{\link{CMBDat}} object
#'
#' This function returns a
#' data.frame containing the data in \code{cmbdat} restricted to the
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
#'consists of the points of \code{cmbdat$data} in the union of
#'\eqn{A} and \eqn{B}.
#'
#'Note that if \eqn{A} (resp. \eqn{B}) is empty then the returned data.frame
#'will be the points of \code{cmbdat} in \eqn{B} (resp. \eqn{A}).
#'
#'@param cmbdat a \code{\link{CMBDat}} object.
#'@param new.window A single \code{\link{CMBWindow}} object or a list of them.
#'@param intersect A boolean that determines
#'the behaviour when \code{win} is a list containing BOTH
#'regular type and "minus" type windows together (see details).
#'
#'@return
#' A CMBDataFrame containing the data in \code{cmbdat} restricted to the
#' CMBWindow \code{new.window}
#'
#'@examples
#'
#'win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
#'cmbdat <- CMBReadFITS("CMB_map_smica1024.fits", mmap = TRUE)
#'class(cmbdat)
#'cmbdat.win <- window(cmbdat, new.window = win1)
#'class(cmbdat.win)
#'
#'@export
window.CMBDat <- function(cmbdat, new.window, intersect = TRUE)
{
  return(subWindow(cmbdat, win = new.window, intersect = intersect))
}




