# R function for creating a CMBDataFrame

############### TO DO: #################################
# 1. If ordering is given and CMBData is path to FITS then
#     we should convert the odering to the desired ordering after
#     importing in whatever order the FITS specifies
# 2. Remove NULL arguments and use missing() instead



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
#'must equal \eqn{12*nside^2} if either spix or
#'sample.size are unspecified.
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
#' # Specify a vector of pixel indices using spix
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

  ### ----- PREPARATION AND CHECKING ARGUMENTS ARE VALID ----- ###

  if ( !missing(win) )
  {
    if ( !rcosmo::is.CMBWindow(win) ) {

      if ( !is.list(win) )
      {
        stop("'win' must be a CMBWindow or list of CMBWindows")
      }

      if ( !all(sapply(win, rcosmo::is.CMBWindow)) )
      {
        stop("'win' must be a CMBWindow or list of CMBWindows")
      }
    }
  }

  if ( !missing(ordering) )
  {
    ordering <- tolower(ordering)
  }

  if ( !missing(spix) && !missing(sample.size) )
  {
    stop("At least one of spix or sample.size should be unspecified")
  }

  if ( !missing(spix) )
  {
    if ( is.character(spix) )
    {
      # If spix is a string then assume it is a path to file:
      message("Reading sample pixel indices from file...\n")
      spix <- read.table(spix, sep = ",")[,1]
    }

    if( any(spix %% 1 != 0) )
    {
      stop("Sample pixel indices must be integers")
    }

    spix <- sort(as.integer(spix))
  }


  if ( !missing(coords) ) {
    coords <- tolower(coords)

    if ( coords != "spherical" && coords != "cartesian" )
    {
      stop(paste("Invalid argument coords must be unspecified,",
                 "spherical or cartesian"))
    }
  }

  if ( missing(spix) )
  {
    spix <- NULL
  }
  else
  {
    len <- length(spix)
    spix <- as.integer(spix)
  }
  if ( !missing(sample.size) )
  {
    len <- sample.size
  }

  ## CASE 1: CMBData is a path
  CMBData.is.path <- FALSE
  if (!missing(CMBData) && is.character(CMBData))
  {
    CMBData.is.path <- TRUE

    if ( !missing(nside) )
    {
      stop("nside must be unspecified when 'CMBData' is specified")
    }

    if ( !missing(intensities) )
    {
      stop("intensities must be unspecified when 'CMBData' is specified")
    }
  }

  ## CASE 2: CMBData is a CMBDataFrame
  CMBData.is.cmbdf <- FALSE
  if ( !missing(CMBData) && is.CMBDataFrame(CMBData) )
  {
    CMBData.is.cmbdf <- TRUE

    if ( !missing(nside) )
    {
      stop("nside must be unspecified when 'CMBData' is specified")
    }

    if ( !missing(intensities) )
    {
      stop("intensities must be unspecified when 'CMBData' is specified")
    }

    if ( include.polar != FALSE || include.masks != FALSE )
    {
      stop(paste("include.polar and include.masks must be",
                 "FALSE unless a FITS file is given"))
    }

  }

  ## CASE 3: CMBData is unspecified
  if ( missing(CMBData) )
  {
    if ( missing(nside) || missing(ordering) )
    {
      stop(paste("If 'CMBData' is unspecified then 'nside'",
                 "and 'ordering' must both be specified"))
    }

    if ( !missing(intensities) )
    {
      if ( is.null(spix) && missing(sample.size) )
      {
        if ( length(intensities) != 12*nside^2 )
        {
          stop(paste("The intensities parameter must have length 12*nside^2",
                     "unless spix or sample.size is specified"))
        }
      }
    }
  }

  ### ---------------------------------------------------  ###



  ########################################################################
  ##### CASE 1: CMBData is a path to a FITS file                 #########
  ########################################################################
  if (CMBData.is.path) {

    CMBData <- CMBReadFITS(CMBData)

    # Get Nside from FITS header:
    nside <- as.numeric(CMBData$hdr[which(CMBData$hdr == "NSIDE")+1])
    # Check that nside is an integer and greater than 0:
    if(nside %% 1 != 0 || (nside <= 0))
    {
      stop("Failed to obtain valid nside from FITS header")
    }

    if ( !missing(sample.size) )
    {
      spix <- sort(sample(seq(1,12*nside^2), sample.size))
    }

    # Get ordering from FITS header:
    orderFITS <- tolower(CMBData$hdr[which(CMBData$hdr == "ORDERING")+1])
    if(orderFITS != "ring" && orderFITS != "nested")
    {
      stop(paste("Failed to obtain valid ordering scheme from FITS header,",
           "instead obtained: ", orderFITS))
    }

    if (!missing(coords)) {

      message("Generating coordinates from HEALPix ordering...\n")

      nest <- ifelse(orderFITS == "nested", TRUE, FALSE)
      cartesian <- ifelse(coords == "cartesian", TRUE, FALSE)

      # generate the coordinates from HEALPix indices
      coordinates <- rcosmo::pix2coords_internal(nside = nside, nested = nest,
                                          spix = spix, cartesian = cartesian)

      # Put the coordinates in a data.frame
      if (coords == "spherical"){
        cmbdf <- data.frame(theta = coordinates[,1], phi = coordinates[,2])
      } else {
        cmbdf <- data.frame(x = coordinates[,1], y = coordinates[,2], z = coordinates[,3])
      }

      # Add the corresponding intensities from CMBData into the data.frame
      if ( !is.null(spix) ){
        cmbdf <- data.frame(cmbdf, I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(cmbdf, I = CMBData$col[[1]])
      }

    # Else coords are unspecified (HEALPix)
    } else {

      if ( !is.null(spix) ){
        cmbdf <- data.frame(I = CMBData$col[[1]][spix])
      } else {
        cmbdf <- data.frame(I = CMBData$col[[1]])
      }

    }

    if (include.polar == TRUE) {
      if (!is.null(spix)) stop(paste("(development stage) include.polar must",
                                    "be FALSE if spix is specified"))
      cmbdf$Q <- CMBData$col[[2]]
      cmbdf$U <- CMBData$col[[3]]
    }

    if (include.masks == TRUE) {
      if (!is.null(spix)) stop(paste("(development stage) include.masks must",
                                     "be FALSE if spix is specified"))
      cmbdf$TMASK <- CMBData$col[[4]]
      cmbdf$PMASK <- CMBData$col[[5]]
    }

    message("Adding CMB Data Frame attributes...\n")

    if ( is.null(spix) ) spix <- 1:(12*nside^2)
    attr(cmbdf, "row.names") <- spix
    attr(cmbdf, "nside") <- nside
    class(cmbdf) <- c("CMBDataFrame","data.frame")
    attr(cmbdf, "ordering") <- orderFITS
    if (missing(coords)) coords <- NULL
    attr(cmbdf, "coords") <- coords

    if (!missing(ordering))
    {
      ordering(cmbdf) <- ordering
    }

  ##############################################################
  ###### CASE 2: CMBData is a CMBDataFrame             #########
  ##############################################################
  } else if ( CMBData.is.cmbdf ) {

    nside <- rcosmo::nside(CMBData)
    n <- nrow(CMBData)

    if (( !missing(sample.size) || !is.null(spix) ) && len > n)
    {
      stop("sample.size or length(spix) exceeds number of rows of 'CMBData'")
    }

    if ( !missing(sample.size) )
    {
      spix <- sort(sample(seq(1,n), sample.size))
    }

    if (!is.null(spix))
    {
      cmbdf <- CMBData[spix,]
    } else {
      cmbdf <- CMBData
    }

    if (!missing(ordering))
    {
      ordering(cmbdf) <- ordering
    }

    if (!missing(coords))
    {
      coords(cmbdf) <- coords
    }


  ################################################################
  ###### CASE 3: CMBData is unspecified             ##############
  ################################################################
  } else if ( missing(CMBData) ) {

    if ( !missing(sample.size) )
    {
      spix <- sort(sample(seq(1,12*nside^2), sample.size))
    }
    else if ( is.null(spix) )
    {
      len <- 12*nside^2
    }

    if ( !is.null(spix) && !missing(intensities)
         && length(intensities) != length(spix) )
    {
      intensities <- intensities[spix]
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

    if ( is.null(spix) ) spix <- as.integer(1:(12*nside^2))
    attr(cmbdf, "row.names") <- spix
    attr(cmbdf, "nside") <- nside
    class(cmbdf) <- c("CMBDataFrame","data.frame")
    attr(cmbdf, "ordering") <- ordering
    if ( missing(coords) ) coords <- NULL
    attr(cmbdf, "coords") <- coords


  } else {

    stop("CMBData must be a CMBDataFrame, a path to a FITS file, or unspecified")

  }

  if ( !missing(win) )
  {
    window(cmbdf) <- win
  }

  return(cmbdf)
}

