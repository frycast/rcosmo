# R function for creating a CMBDataFrame

# CMBData can be the output of CMBReadFITS or it can be a string location of FITS file.

# Coords can be "spherical," "healpix only" or "cartesian".

### ---- NOT SURE WHY THE FOLLOWING USED TO BE HERE -------- ###
### ---- BUT PRETTY SURE THAT IT IS UNNECESSARY ------------ ###
### ---- BECAUSE pix2coords is NOT USED WHEN coords = unspecified --- ###
### ---- (pix2coords WAS REORDERING COORDINATES) ------------- ###
# If healpix then spix must be unspecified.
# try(if(coords == "healpix" && !missing(spix))
#   stop("When spix is not missing coords must be either cartesian or spherical"))
### --------------------------------------------------- ###

############### TO DO: #################################
# 1. If ordering is given and CMBData is path to FITS then
#     we should convert the odering to the desired ordering after
#     importing in whatever order the FITS specifies
# 2. Remove NULL arguments and use missing() instead

# NOTE: THE VALUE IN NOT FORCING spix TO BE ORDERED IS THAT
# SPIX CAN BE USED TO CONVERT THE ORDERING SCHEME OF THE WHOLE
# CMBDataFrame?? THE CODE IS STILL THERE COMMENTED OUT BUT
# I DONT THINK WE NEED spix TO BE UNORDERED

# If the spix argument takes a vector of sample pixel indices then map is assumed not fullSky.
# spix can also take a path to csv file of indices.
# Ideally specifying spix will allow user to subset the data and we could include an extra argument
# that allows the user to subset the data using coordinate subset.


#' CMB Data Frames
#'
#' The function \code{CMBDataFrame} creates CMB Data Frames. These are
#' a special type of \link{data.frame} that carry extra information
#' about the HEALPix ordering scheme, coordinate system, and nside parameter.
#'
#'@param CMBData can be a string location of FITS file,
#'another CMBDataFrame, or nothing.
#'@param coords can be "spherical," "cartesian", or unspecified (HEALPix only).
#'@param win optional \code{\link{CMBWindow}} object that specifies a
#'spherical polygon within which to subset the full sky CMB data
#'@param include.polar TRUE if polarisation data is required, otherwise FALSE
#'@param include.masks TRUE if TMASK and PMASK are required, otherwise FALSE
#'@param spix optional vector of sample pixel indices or a path to a file
#'containing comma delimited sample pixel indices. The ordering scheme is
#'given by \code{ordering}. If \code{ordering} is unspecified then
#'CMBData must be either a CMBDataFrame or a FITS file and the ordering
#'scheme is then assumed to match that of CMBData.
#'@param sample.size if a positive integer is given, a simple random
#'sample of size equal to sample.size will be taken from CMBData. If
#'spix is specified then \code{sample.size} must be unspecified.
#'@param nside optionally specify the nside parameter manually
#'(usually 1024 or 2048)
#'@param ordering specifies the desired HEALPix ordering scheme
#'("ring" or "nested") for the output CMBDataFrame.
#'If \code{ordering} is unspecified then the ordering
#'scheme will be taken from the CMBData object, which must then be
#'either a CMBDataFrame or a path to a FITS file. This parameter also specifies
#'the ordering scheme of \code{spix}.
#'@param intensities a vector of intensities to be included
#'if \code{CMBData} is unspecified. Note that \code{length(intensities)}
#'must equal \eqn{12*nside^2} if spix and
#'sample.size are unspecified, otherwise \code{length(intensities)}
#'must equal \code{length(spix)} or \code{length(sample.size)}
#'
#'@return
#'A data frame whose columns contain the pixel center coordinates
#'theta, phi (meaning colatitude in range \eqn{[0,pi]} and longitude
#'in range \eqn{[0,2pi)} respectively)
#'or (x,y,z), CMB intensities (I), and
#'optionally polarisation (Q,U) and masks (TMASK, PMASK).
#'The row.names attribute of the resulting CMB Data Frame contains
#'HEALPix indices.
#'
#'@examples
#' ## Method 1: Read the data while constructing the CMBDataFrame
#' df <- CMBDataFrame("CMB_map_smica1024.fits")
#'
#' # Specify a sample size for a random sample
#' df.sample <- CMBDataFrame(df, sample.size = 800000)
#' plot(df.sample)
#'
#' # Specify a vector of pixel indices to keep, using spix
#' df.subset <- CMBDataFrame(df, spix = c(2,4,6))
#'
#' # Take a look at the summary
#' summary(df)
#'
#' # Access HEALPix pixel indices using pix function
#' # (these are stored in the row.names attribute)
#' pix(df)
#'
#'@export
CMBDataFrame <- function(CMBData,
                         coords,
                         win,
                         include.polar = FALSE,
                         include.masks = FALSE,
                         spix,
                         sample.size,
                         nside,
                         ordering,
                         intensities) {

  ### --- PREPARATION AND CHECKING ARGUMENTS ARE VALID --- ###

  if ( !missing(win) )
  {
    if ( !rcosmo::is.CMBWindow(win) ) {

      if ( !is.list(win) )
      {
        stop("'win' must be a CMBWindow or list of CMBWindows")
      }

      if (!all(sapply(win, rcosmo::is.CMBWindow)))
      {
        stop("'win' must be a CMBWindow or list of CMBWindows")
      }
    }
  }

  # If spix is a string then assume it is a path to file:
  if (!missing(spix) && is.character(spix)) {
    message("Reading sample pixel indices from file...\n")
    spix <- read.table(spix, sep = ",")[,1]
  }

  if (!missing(spix)) {
    try(if(any(spix %% 1 != 0))
      stop("Sample pixel indices must be integers"))
    try(if(!missing(intensities) && length(intensities) != length(spix))
      stop("Intensities parameter must have same length as spix parameter"))
    try(if(!missing(sample.size))
      stop("At least one of spix or sample.size should be unspecified"))

    spix <- sort(as.integer(spix))
  }

  if (!missing(sample.size)) {
    try(if(!missing(spix))
      stop("At least one of spix or sample.size should be unspecified"))
    try(if(!missing(intensities)
           && length(intensities) != length(sample.size))
      stop(paste("Intensities parameter must have",
           "same length as sample.size parameter")))
  }

  if ( !missing(coords) ) {
    coords <- tolower(coords)
  }

  try(if( !missing(coords)
         && coords != "spherical"
         && coords != "cartesian")
        stop(paste("Invalid argument coords must be unspecified,",
             "spherical or cartesian")))

  ### ---------------------------------------------------  ###


  ########################################################################
  ##### If CMBData is a string then assume it is a path to file: #########
  ########################################################################
  if (!missing(CMBData) && is.character(CMBData)) {

    try(if (!missing(nside))
      stop("nside must be unspecified if CMBData is a path to a FITS file"))

    try(if (!missing(intensities))
      stop(paste("intensities parameter must be unspecified when",
           "CMBData is a path to FITS file")))

    CMBData <- CMBReadFITS(CMBData)

    # Get Nside from FITS header:
    nside <- as.numeric(CMBData$hdr[which(CMBData$hdr == "NSIDE")+1])
    # Check that nside is an integer and greater than 0:
    try(if(nside %% 1 != 0 || (nside <= 0))
      stop("Failed to obtain valid nside from FITS header"))


    # Get ordering from FITS header:
    orderFITS <- tolower(CMBData$hdr[which(CMBData$hdr == "ORDERING")+1])
    try(if(orderFITS != "ring" && orderFITS != "nested")
      stop(paste("Failed to obtain valid ordering scheme from FITS header,",
           "instead obtained: ", orderFITS)))


    if ( !missing(ordering) ) {

      ordering <- tolower(ordering)

      if ( ordering != orderFITS && !missing(spix) )
      {
        # The ordering scheme of spix is assumed to be specified
        # by the ordering parameter, and the output CMBData frame
        # will also have ordering scheme given by the ordering parameter.
        # So, Here is where we must
        # do a conversion of spix into the orderFITS ordering scheme.
        # Further down, when it is safe to do so,
        # the entire output CMBDataFrame can be converted to
        # have ordering scheme specified by the ordering parameter.
        stop(paste("ordering (",ordering,") does not match the ordering scheme",
             "(",orderFITS,") given in the FITS file.\n",
             "This error is temporary (rcosmo development stage) and is",
             "only here because the necessary conversions have not been",
             "implemented"))
      }

    }


    # The total number of pixels is equal to 12*nside^2,
    # unless a sample is specified:
    Npix <- ifelse(!missing(spix), length(spix), 12*nside^2)

    if (!missing(sample.size)) {
      spix <- sample(seq(1,12*nside^2), sample.size)
      Npix <- sample.size
    }

    if (!missing(coords)) {

      message("Generating coordinates from HEALPix ordering...\n")

      nest <- ifelse(orderFITS == "nested", TRUE, FALSE)
      cartesian <- ifelse(coords == "cartesian", TRUE, FALSE)

      # generate the coordinates from HEALPix indices
      if (missing(spix))
      {
        coordinates <- rcosmo::pix2coords_internal(nside = nside, nested = nest,
                                          spix = NULL, cartesian = cartesian)
      } else {

        coordinates <- rcosmo::pix2coords_internal(nside = nside, nested = nest,
                                          spix = spix, cartesian = cartesian)
      }

      # Put the coordinates in a data.frame
      if (coords == "spherical"){
        cmbdf <- data.frame(theta = coordinates[,1], phi = coordinates[,2])
      } else {
        cmbdf <- data.frame(x = coordinates[,1], y = coordinates[,2], z = coordinates[,3])
      }

      # Add the corresponding intensities from CMBData into the data.frame
      if (!missing(spix)){
        cmbdf <- data.frame(cmbdf, I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(cmbdf, I = CMBData$col[[1]])
      }

    # Else coords are unspecified (HEALPix)
    } else {

      if (!missing(spix)){
        cmbdf <- data.frame(I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(I = CMBData$col[[1]])
      }

    }

    if (include.polar == TRUE) {
      if (!missing(spix)) stop(paste("(development stage) include.polar must",
                                    "be FALSE if spix is specified"))
      cmbdf$Q <- CMBData$col[[2]]
      cmbdf$U <- CMBData$col[[3]]
    }

    if (include.masks == TRUE) {
      if (!missing(spix)) stop(paste("(development stage) include.masks must",
                                     "be FALSE if spix is specified"))
      cmbdf$TMASK <- CMBData$col[[4]]
      cmbdf$PMASK <- CMBData$col[[5]]
    }

    message("Adding CMB Data Frame attributes...\n")

    if (!missing(spix)){
      attr(cmbdf, "row.names") <- as.integer(spix)
    }
    attr(cmbdf, "nside") <- nside
    class(cmbdf) <- c("CMBDataFrame","data.frame")
    attr(cmbdf, "ordering") <- orderFITS
    if (missing(coords)) {

      attr(cmbdf, "coords") <- NULL

    } else {

      attr(cmbdf, "coords") <- coords

    }


    # Here is where the cmbdf ordering will be converted if needed
    if ( !missing(ordering) ) { ordering(cmbdf) <- ordering }

  ##############################################################
  ###### Otherwise, is CMBData is a CMBDataFrame? ##############
  ##############################################################
  } else if ( !missing(CMBData) && is.CMBDataFrame(CMBData) ) {

    try(if ( include.polar != FALSE || include.masks != FALSE )
      stop(paste("include.polar and include.masks must be",
           "FALSE unless a FITS file is given")))

    try(if ( !missing(intensities) )
      stop(paste("when CMBData is specified the intensities",
           "parameter must be unspecified")))

    try(if (!missing(nside))
      stop(paste("the nside parameter must be unspecified",
           "if CMBData is a CMBDataFrame")))

    if ( !missing(sample.size) )
    {
      message("Taking a simple random sample...\n")
      cmbdf <- sampleCMB(CMBData, sample.size = sample.size)
    } else
    {
      cmbdf <- CMBData
    }

    if (!missing(ordering))
    {
      ordering(cmbdf) <- ordering
    }

    if (!missing(spix))
    {
      cmbdf <- cmbdf[spix,]
    }

    if (!missing(coords))
    {
      coords(cmbdf) <- coords
    }

  ################################################################
  ###### Otherwise, check if CMBData is unspecified ##############
  ################################################################
  } else if ( missing(CMBData) ) {

    if ( missing(nside) || missing(ordering)  ) {
      stop(paste("If CMBData is unspecified then nside",
           "and ordering must both be specified"))
    }

    if ( !missing(spix) )
    {
      len <- length(spix)
      try( if ( !missing(intensities) && length(intensities) != len )
        stop(paste("intensities parameter and spix",
             "parameter must have same length")))

    } else if ( !missing(sample.size) ) {

      spix <- NULL
      len <- length(sample.size)
      try( if ( !missing(intensities) &&  length(intensities) != len )
        stop(paste("intensities parameter and sample.size",
             "parameter must have same length")))

      stop(paste("(development stage) if CMBData is unspecified",
           "then sample.size must be unspecified"))


    } else {

      spix <- NULL
      len <- 12*nside^2
      try( if ( !missing(intensities) && length(intensities) != len )
        stop("intensities parameter must be of length 12*nside^2"))

    }

    if ( missing(intensities) )
    {
      intensities <- rep(NA, len)
    }


    if ( !missing(coords) )
    {
      message("Generating coordinates from HEALPix ordering...\n")

      cmbdf <- rcosmo::pix2coords(nside = nside, ordering = ordering,
                                  coords = coords, spix = spix)

      cmbdf <- cbind(cmbdf, I = intensities)

    } else {

      cmbdf <- data.frame(I = intensities)

    }


    message("Adding CMB Data Frame attributes...\n")

    if ( !missing(spix) && !is.null(spix) ){
      attr(cmbdf, "row.names") <- as.integer(spix)
    }
    attr(cmbdf, "nside") <- nside
    class(cmbdf) <- c("CMBDataFrame","data.frame")
    attr(cmbdf, "ordering") <- ordering
    if ( missing(coords) ) {
      attr(cmbdf, "coords") <- NULL
    } else {
      attr(cmbdf, "coords") <- coords
    }



  } else {

    stop("CMBData must be a CMBDataFrame, a path to a FITS file, or unspecified")

  }

  if ( !missing(win) )
  {
    window(cmbdf) <- win
  }

  return(cmbdf)
}

