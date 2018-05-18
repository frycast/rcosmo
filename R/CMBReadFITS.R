#' Read CMB data from a FITS file.
#'
#' \code{CMBReadFITS} is adapted from the \link{readFITS} function in package
#'   \href{https://cran.r-project.org/web/packages/FITSio/index.html}{FITSio}.
#'   \code{CMBReadFITS} is in development stage and will only work with 'CMB_map_smica1024.fits'.
#'   When it works, \code{CMBReadFITS} is much faster than \code{readFITS}.
#'   However, \code{readFITS} is more general and so is more likely to work.
#'
#'@param filename The path to the fits file.
#'@return A list containing the intensity (I), polarisation (Q, U), PMASK, TMASK and metadata.
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
CMBReadFITS <- function(filename = "CMB_map_smica1024.fits", sample.size) {

zz <- file(filename, "rb")
header1 <- FITSio::readFITSheader(zz)
header2 <- FITSio::readFITSheader(zz)
hdr <- FITSio::parseHdr(header2)

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


# Prepare column lists
col <- rep(list(array(NA, dim = c(naxis2, 1))), tfields)


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
for (i in 1:naxis2) {
  for (j in 1:tfields) {
    col[[j]][i, ] <- switch(btype[j],
                            .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)),   # PMASK,TMASK
                            .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) )    # I,Q,U
  }
}

close(zz)


names(col) <- TTYPEn

# Return a list and close the file  -----------------------------------------------------
cmbdat <- list(col = col, hdr = hdr, resoln = resoln, method = method, coordsys = coordsys,
               ordering = ordering, nside = nside, baddata = baddata, header1 = header1,
               header2 = header2)


return(cmbdat)
}











#' Read CMB data from a FITS file.
#'
#' \code{CMBReadFITS} is adapted from the \link{readFITS} function in package
#'   \href{https://cran.r-project.org/web/packages/FITSio/index.html}{FITSio}.
#'   \code{CMBReadFITS} is in development stage and will only work with 'CMB_map_smica1024.fits'.
#'   When it works, \code{CMBReadFITS} is much faster than \code{readFITS}.
#'   However, \code{readFITS} is more general and so is more likely to work.
#'
#'@param filename The path to the fits file.
#'@param sample.size If specified then a simple random sample of the data will loaded
#'@return A list containing the intensity (I), polarisation (Q, U), PMASK, TMASK and metadata.
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
CMBReadFITS2 <- function(filename = "CMB_map_smica1024.fits", sample.size) {

  zz <- file(filename, "rb")
  header1 <- FITSio::readFITSheader(zz)
  header2 <- FITSio::readFITSheader(zz)
  hdr <- FITSio::parseHdr(header2)

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

  if (!missing(sample.size))
  {
    pix <- sample(1:naxis2, size = sample.size)
    col <- rep(list(array(NA, dim = c(sample.size, 1))), tfields)

    for (i in 1:naxis2) {
      for (j in 1:tfields) {
        element <- switch(btype[j],
                          .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)),   # PMASK,TMASK
                          .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) )    # I,Q,U
        if (i %in% pix)
        {
          col[[j]][i, ] <- element
        }
      }
    }
  }
  else
  {
    col <- rep(list(array(NA, dim = c(naxis2, 1))), tfields)

    for (i in 1:naxis2) {
      for (j in 1:tfields) {
        col[[j]][i, ] <- switch(btype[j],
                                .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)),   # PMASK,TMASK
                                .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) )    # I,Q,U
      }
    }
  }


  close(zz)


  names(col) <- TTYPEn

  # Return a list and close the file  -----------------------------------------------------
  if (!missing(sample.size))
  {
    cmbdat <- list(col = col, hdr = hdr, resoln = resoln, method = method, coordsys = coordsys,
                   ordering = ordering, nside = nside, baddata = baddata, header1 = header1,
                   header2 = header2)
  }
  else
  {
        cmbdat <- list(col = col, pix = pix, hdr = hdr, resoln = resoln, method = method,
                       coordsys = coordsys, ordering = ordering, nside = nside, baddata = baddata,
                       header1 = header1, header2 = header2)
  }


  return(cmbdat)
}
