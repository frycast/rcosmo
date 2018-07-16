






############################################################################
## WARNING: THIS FUNCTION IS INCOMPLETE: THE mmap PARAMETER REQUIRES WORK ##
############################################################################
#' Read CMB data from a FITS file.
#'
#' \code{CMBReadFITS} is adapted from the \code{\link{readFITS} }function in package
#'   \href{https://cran.r-project.org/web/packages/FITSio/index.html}{FITSio}.
#'   \code{CMBReadFITS} is in development stage and will only work with 'CMB_map_smica1024.fits'.
#'   When it works, \code{CMBReadFITS} is much faster than \code{\link{readFITS}}.
#'   However, \code{\link{readFITS}} is more general and so is more likely to work.
#'
#'@param filename The path to the fits file.
#'@param mmap A boolean indicating whether to use memory mapping.
#'@return A list containing header information and other metadata
#'as well as an element called \code{data} where:
#'If \code{mmap = FALSE} then a \code{data.frame} is
#'included, named \code{data}, whose columns may include, for
#'example, the intensity (I), polarisation (Q, U), PMASK and TMASK.
#'If \code{mmap = TRUE} then a \code{\link{mmap}} object is returned
#'that points to the CMB map data table in the FITS file.
#'@examples
#'dat <- CMBReadFITS("CMB_map_smica1024.fits")
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
CMBReadFITS <- function(filename, mmap = FALSE) {

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



  # if ( mmap == FALSE )
  # {
  #   # col <- as.numeric(rep(btype, naxis2))
  #   # len <- naxis2*tfields
  #   # for (i in 1:len){
  #   #     col[i] <- switch(col[i],
  #   #                      .Internal(readBin(zz, "integer", 1L,
  #   #                                        1L, FALSE, swap)), # PMASK,TMASK
  #   #                      .Internal(readBin(zz, "double", 1L,
  #   #                                        4L, TRUE, swap)) ) # I,Q,U
  #   # }
  #   # col <- matrix(col, nrow = naxis2, byrow = TRUE)
  #   # col <- as.data.frame(col)
  # }
  # else
  # {}

  mystruct <- do.call(mmap::struct, CTypeExpression(TTYPEn, btype))
  map <- mmap::mmap(file = filename,
                    mode = mystruct,
                    off = blocks*bytes,
                    endian = "big")
  mmap::extractFUN(map) <- function(X) do.call(data.frame, X)

  if ( mmap == FALSE )
  {
    col <- map[1:(12*nside^2)]
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
