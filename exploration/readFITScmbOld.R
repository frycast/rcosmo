readFITScmb <- function(filename = "CMB_map_smica1024.fits") {

zz <- file(filename, "rb")
# The repetition is necessary:
header <- readFITSheader(zz)
hdr <- parseHdr(header)
# tmp <- hdr[which(hdr == "NAXIS") + 1] #Check this is always "0"
header <- readFITSheader(zz)
hdr <- parseHdr(header)
# hdr[which(hdr == "XTENSION") + 1] #Check this is always 'BINTABLE'
## The following is executed in "D <- readFITSbintable(zz, hdr)":
# nchar(hdr[1]) check this != 80
naxis1 <- as.numeric(hdr[which(hdr == "NAXIS1") + 1])
naxis2 <- as.numeric(hdr[which(hdr == "NAXIS2") + 1])
tfields <- as.numeric(hdr[which(hdr == "TFIELDS") + 1])
# ifelse(length(hdr[which(hdr == "PCOUNT") + 1]) != 1, 0, as.numeric(tmp)) #Check this == 0
# ifelse(length(hdr[which(hdr == "GCOUNT") + 1]) != 1, 1, as.numeric(tmp)) #Check this != 1






# Extract metadata ---------------------------------------------------
TFORMn <- character(tfields)
TTYPEn <- character(tfields) # This one can let us select I,Q and U.
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




# Set up byte flags------------------------------------------------------
# bsize <- integer(tfields)    #Check: 4 4 4 1 1
# btype <- integer(tfields)    #Check: 4 4 4 3 3
# bsign <- logical(tfields)    #Check: all FALSE
# mult <- integer(tfields)     #Check: 1 1 1 1 1 
# for (i in 1:tfields) {
#   nc <- nchar(TFORMn[i])
#   tmp <- substr(TFORMn[i], 1, nc - 1)
#   mult[i] <- ifelse(tmp == "", 1, as.numeric(tmp))
#   form <- tolower(substr(TFORMn[i], nc, nc))
#   switch(form, l = {
#     bsize[i] <- 1
#     btype[i] <- 2
#     bsign[i] <- FALSE
#   }, b = {              # Check TFORMn always 'B' for PMASK, TMASK
#     bsize[i] <- 1
#     btype[i] <- 3
#     bsign[i] <- FALSE
#   }, i = {
#     bsize[i] <- 2
#     btype[i] <- 3
#     bsign[i] <- TRUE
#   }, j = {
#     bsize[i] <- 4
#     btype[i] <- 3
#     bsign[i] <- TRUE
#   }, k = {
#     bsize[i] <- 8
#     btype[i] <- 6
#     bsign[i] <- TRUE
#   }, a = {
#     bsize[i] <- 1
#     btype[i] <- 1
#   }, e = {             # Check TFORMn is always 'E' for I,Q,U.
#     bsize[i] <- 4
#     btype[i] <- 4
#   }, d = {
#     bsize[i] <- 8
#     btype[i] <- 4
#   }, stop("Unknown TFORMn ", toupper(form), " (X, C, M, P not yet implemented)\n"))
# }





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






# The guts: read the data --------------------------------------------------
# if (btype[j] <= 5)... #Check btype <= 5

# ORIGINAL SLOW VERSION
# for (i in 1:naxis2) {
#   for (j in 1:tfields) {
#     if (btype[j] <= 5) {
#       col[[j]][i, ] <- switch(btype[j], 
#                               readChar(zz, nchars = mult[j]), 
#                               readChar(zz, nchars = mult[j]), 
#                               readBin(zz, what = integer(), n = mult[j], size = bsize[j], signed = bsign[j], endian = "big"), 
#                               readBin(zz, what = numeric(), n = mult[j], size = bsize[j], endian = "big"), # Used for TMASK, PMASK
#                               readBin(zz, what = complex(), n = mult[j], size = bsize[j], endian = "big")) # Used for I,Q,U
#     }
#   }
# }
# FASTER VERSION TAYLORED TO CMB
# seek(zz, where = 8640) reset position (debugging)
swap <- "big" != .Platform$endian
for (i in 1:naxis2) {
  col[[1]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # I
  col[[2]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # Q
  col[[3]][i, ] <- .Internal(readBin(zz, "double", 1, 4, TRUE, swap)) # U
  col[[4]][i, ] <- .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)) # TMASK
  col[[5]][i, ] <- .Internal(readBin(zz, "integer", 1, 1, FALSE, swap)) # PMASK
}





# Scaling etc -------------------------------------------------------------
# nrow(col[[i]]) == 1 || ncol(col[[i]]) == 1 || btype[i] <= 2)   #Check this is FALSE
# btype[i] >= 3 && (TSCALn[i] != 1 || TZEROn[i] != 0)            #Check this is FALSE
# done (no scaling)




# Return a neat list and close the file  -----------------------------------------------------
cmbdat <- list(col = col, hdr = hdr, colNames = TTYPEn, colUnits = TUNITn, 
               TNULLn = TNULLn, TSCALn = TSCALn, TZEROn = TZEROn, TDISPn = TDISPn)
close(zz)
cmbdat$header <- header
return(cmbdat)
}