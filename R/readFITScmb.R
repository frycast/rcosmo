#' Read CMB data from a FITS file.
#'
#' \code{readFITScmb} is adapted from the \link{readFITS} function in package
#'   \href{https://cran.r-project.org/web/packages/FITSio/index.html}{FITSio}.
#'   \code{readFITScmb} is in development stage and will only work with 'CMB_map_smica1024.fits'.
#'   When it works, \code{readFITScmb} is much faster than \code{readFITS}.
#'   However, \code{readFITS} is more general and so is more likely to work.
#'
#' @param filename The path to the fits file.
#' @return A list containing the intensity (I), polarisation (Q, U), PMASK, TMASK and metadata.
#' @examples
#'   readFITScmb("CMB_map_smica1024.fits")
#'

readFITScmb <- function(filename = "CMB_map_smica1024.fits") {

zz <- file(filename, "rb")

header <- readFITSheader(zz)
hdr <- parseHdr(header)
header <- readFITSheader(zz)
hdr <- parseHdr(header)
naxis1 <- as.numeric(hdr[which(hdr == "NAXIS1") + 1])
naxis2 <- as.numeric(hdr[which(hdr == "NAXIS2") + 1])
tfields <- as.numeric(hdr[which(hdr == "TFIELDS") + 1])




# Extract metadata ---------------------------------------------------
TFORMn <- character(tfields)
TTYPEn <- character(tfields)
TUNITn <- character(tfields)
TNULLn <- integer(tfields)
TSCALn <- numeric(tfields)
TZEROn <- numeric(tfields)
TDISPn <- character(tfields)
THEAPn <- integer(tfields)
TDIMn <- character(tfields)
for (i in 1:tfields) {
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TFORM",
                                              i, sep = "")) + 1])
  TFORMn[i] <- tmp
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TTYPE",
                                              i, sep = "")) + 1])
  TTYPEn[i] <- ifelse(length(tmp) != 1, "", tmp)
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TUNIT",
                                              i, sep = "")) + 1])
  TUNITn[i] <- ifelse(length(tmp) != 1, "", tmp)
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TNULL",
                                              i, sep = "")) + 1])
  TNULLn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TSCAL",
                                              i, sep = "")) + 1])
  TSCALn[i] <- ifelse(length(tmp) != 1, 1, as.numeric(tmp))
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TZERO",
                                              i, sep = "")) + 1])
  TZEROn[i] <- ifelse(length(tmp) != 1, 0, as.numeric(tmp))
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TDISP",
                                              i, sep = "")) + 1])
  TDISPn[i] <- ifelse(length(tmp) != 1, "", tmp)
  tmp <- gsub(" ", "", hdr[which(hdr == paste("THEAP",
                                              i, sep = "")) + 1])
  THEAPn[i] <- ifelse(length(tmp) != 1, NA, as.numeric(tmp))
  tmp <- gsub(" ", "", hdr[which(hdr == paste("TDIM",
                                              i, sep = "")) + 1])
  TDIMn[i] <- ifelse(length(tmp) != 1, "", tmp)
}




# Prepare column lists--------------------------------------------------------
col <- vector("list", tfields)
btype <- c(4,4,4,3,3)
mult <- c(1,1,1,1,1)
for (i in 1:tfields) {
  col[[i]] <- switch(btype[i],
                     array("", dim = c(naxis2, 1)),
                     array(FALSE, dim = c(naxis2, mult[i])),
                     array(NA, dim = c(naxis2, mult[i])),    # Used for PMASK, TMASK
                     array(NA, dim = c(naxis2, mult[i])),    # Used for I,Q,U
                     array(NA, dim = c(naxis2, mult[i])),
                     array(NA, dim = c(2 * naxis2, mult[i])))
}






# The guts: read the binary data --------------------------------------------------
swap <- "big" != .Platform$endian
for (i in 1:naxis2) {
  col[[1]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # I
  col[[2]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # Q
  col[[3]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # U
  col[[4]][i, ] <- .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)) # TMASK
  col[[5]][i, ] <- .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)) # PMASK
}




# Return a neat list and close the file  -----------------------------------------------------
cmbdat <- list(col = col, hdr = hdr, colNames = TTYPEn, colUnits = TUNITn,
               TNULLn = TNULLn, TSCALn = TSCALn, TZEROn = TZEROn, TDISPn = TDISPn)
close(zz)
cmbdat$header <- header
return(cmbdat)
}
