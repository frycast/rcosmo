# R function for creating a CMBDataFrame

# CMBData can be the output of readFITScmb or it can be a string location of FITS file.

# Coords can be "spherical," "HEALPix" or "cartesian". If healpix then spix must be NULL.

# If the spix argument takes a vector of sample pixel indices then map is assumed not fullSky.
# spix can also take a path to csv file of indices.
# Ideally specifying spix will allow user to subset the data and we could include an extra argument
# that allows the user to subset the data using coordinate subset.


#' CMB Data Frames
#'
#' The function \code{CMBDataFrame} creates CMB Data Frames. These are a special type of
#' \link{data.frame} that carry extra information about the HEALPix ordering scheme,
#' coordinate system, and Nside parameter.
#'
#'@param CMBData can be the output of readFITScmb or it can be a string location of FITS file.
#'@param coords can be "spherical," "HEALPix" or "cartesian". If "healpix" then spix must be NULL.
#'@param includePolar TRUE if polarisation data is required, otherwise FALSE
#'@param includeMasks TRUE if TMASK and PMASK are required, otherwise FALSE
#'@param spix a vector of sample pixel indices
#'@param sampleSize if a positive integer is given a simple random sample of size sampleSize
#'will be taken from CMBData
#'@param Nside optionally specify the Nside parameter manually
#'@param ordering optionally specify the HEALPix ordering scheme (RING or NESTED) manually
#'
#'@return
#'A data frame whose columns contain the pixel center coordinates (longitude, lattitude)
#'or (x,y,z), CMB intensities (I), and optionally polarisation (Q,U) and masks (TMASK, PMASK).
#'The row.names attribute of the resulting CMB Data Frame contains HEALPix indices.
#'
#'@examples
#' ## Method 1: Read the data while constructing the CMBDataFrame
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#'
#' ## Method 2: Read the data first separately using readFITScmb
#' cmbdat <- readFITScmb("../../CMB_map_smica1024.fits")
#'
#' # Specify a sample size for a random sample
#' df2 <- CMBDataFrame(CMBData = cmbdat, sampleSize = 800000)
#' plotCMB(df2)
#'
#' # Specify a vector of pixel indices to keep, using spix
#' df3 <- CMBDataFrame(CMBData = cmbdat, spix = c(2,4,6))
#'
#' # Take a look at the attributes
#' attributes(df2)$Nside
#' attributes(df2)$coords
#' attributes(df2)$ordering
#' # The row.names attribute contains HEALPix pixel indices
#' attributes(df3)$row.names
#' head(attributes(df2)$row.names, n = 10)
#' head(attributes(df)$row.names, n = 10)
#'
#'@export
CMBDataFrame <- function(CMBData = "../../CMB_map_smica1024.fits",
                         coords = "spherical",
                         includePolar = FALSE,
                         includeMasks = FALSE,
                         spix = NULL,
                         sampleSize = NULL,
                         Nside = NULL,
                         ordering = NULL) {

  # If CMBData is a string then assume it is a path to file:
  if (is.character(CMBData)) {
    CMBData <- readFITScmb(CMBData)
  }

  # If spix is a string then assume it is a path to file:
  if (is.character(spix)) {
    cat("Reading sample pixel indices from file...")
    cat("\n")
    spix <- read.table(spix, sep = ",")[,1]
  }

  # Check that Nside has not been given explicitly, if not get it from header:
  Nside <- ifelse(is.null(Nside), as.numeric(CMBData$hdr[which(CMBData$hdr == "NSIDE")+1]), Nside)
  # Check that Nside is an integer and greater than 0:
  try(if(Nside %% 1 != 0 || (Nside <= 0)) stop("Failed to obtain valid Nside from FITS header"))

  # The number of pixels is defined to be 12*Nside^2, unless a sample is specified:
  Npix <- ifelse(!is.null(spix), length(spix), 12*Nside^2)

  if (!is.null(sampleSize)) {
    try(if(!is.null(spix)) stop("At least one of spix or sampleSize should be NULL"))
    spix <- sample(seq(1,12*Nside^2),sampleSize)
    Npix <- sampleSize
  }

  if (!is.null(spix)) {
    try(if(any(spix %% 1 != 0)) stop("Sample pixel indices must be integers"))
    spix <- as.integer(spix)
  }

  ordering <- ifelse(is.null(ordering), tolower(CMBData$hdr[which(CMBData$hdr == "ORDERING")+1]), ordering)
  try(if(ordering != "ring" && ordering != "nested") stop("Failed to obtain valid ordering scheme from FITS header"))

  coords <- tolower(coords)
  try(if(coords != "spherical"
         && coords != "cartesian"
         && coords != "healpix") stop("Invalid argument coords"))

  try(if(coords == "healpix" && !is.null(spix)) stop("When spix is not null coords must be either cartesian or spherical"))

  if (coords == "spherical" || coords == "cartesian"){

    if (ordering == "nested"){
      if(is.null(spix)) {
        cat("Generating coordinates from HEALPix ordering...")
        cat("\n")
        sph <- pix2angC(Nside, Nest = TRUE)
      } else {
        cat("Generating coordinates from HEALPix ordering...")
        cat("\n")
        sph <- pix2angC(Nside, Nest = TRUE, spix = spix)
      }
    } else {
      # Else ordering is RING:
      if(is.null(spix)) {
        cat("Generating coordinates from HEALPix ordering...")
        cat("\n")
        sph <- pix2angC(Nside, Nest = FALSE)
      } else {
        cat("Generating coordinates from HEALPix ordering...")
        cat("\n")
        sph <- pix2angC(Nside, Nest = FALSE, spix = spix)
      }
    }

    # Reorder the coordinates to match spix if spix was not null
    if (!is.null(spix)) {
      sph <- sph[match(spix, sph[,5]),]
    }

    if (coords == "spherical"){
      if (!is.null(spix)){
        cmbdf <- data.frame(long = sph[,2], lat = pi/2 - sph[,1], I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(long = sph[,2], lat = pi/2 - sph[,1], I = CMBData$col[[1]])
      }
    } else {
      # Else coords are cartesian:
      xyz <- matrix(c(sph[,2], pi/2 - sph[,1], rep(1,Npix)), nrow = Npix)
      xyz <- sphereplot::sph2car(xyz, deg = FALSE)
      if (!is.null(spix)){
        cmbdf <- data.frame(x = xyz[,1], y = xyz[,2], z = xyz[,3], I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(x = xyz[,1], y = xyz[,2], z = xyz[,3], I = CMBData$col[[1]])
      }
    }

  } else {
    # Else coords are HEALPix
    if (!is.null(spix)){
      cmbdf <- data.frame(I = CMBData$col[[1]][spix])
    } else {
      cmbdf <- data.frame(I = CMBData$col[[1]])
    }
  }

  if (includePolar == TRUE) {
    cmbdf$Q = CMBData$col[[2]]
    cmbdf$U = CMBData$col[[3]]
  }

  if (includeMasks == TRUE) {
    cmbdf$TMASK = CMBData$col[[4]]
    cmbdf$PMASK = CMBData$col[[5]]
  }

  cat("Adding CMB Data Frame attributes...")
  cat("\n")

  if (!is.null(spix)){
    attr(cmbdf, "row.names") <- spix
  }
  attr(cmbdf, "Nside") <- Nside
  attr(cmbdf, "ordering") <- ordering
  attr(cmbdf, "coords") <- coords
  class(cmbdf) <- c("CMBDataFrame","data.frame")
  return(cmbdf)
}
