#' Download CMB Maps from Planck Public Data Release.
#'
#' The function \code{downloadCMBMap} downloads CMB maps from
#' \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}.
#'
#' CMB maps have been produced by the COMMANDER, NILC, SEVEM, and SMICA
#' pipelines, respectively.
#'
#' For each pipeline, the intensity maps are provided at Nside = 2048, at
#' 5 arcmin resolution, and the polarization maps are provided at Nside = 1024
#' at 10 arcmin resolution.
#'
#' @param foreground A string naming the foreground separation method pipeline.
#' Please choose one of "COMMANDER", "NILC", "SEVEM" or "SMICA"
#' (not case sensitive).
#' @param nside An integer. The nside parameter (resolution) required.
#' The available options are \code{1024} or \code{2048}.
#' @param destfile  An optional character string with the path and file
#' name for the downloaded file to be saved. Defaults to the working
#' directory. Tilde-expansion is performed.
#'
#' @return CMB Map FITS File (Flexible Image Transport System). The
#' FITS file can be loaded into a \code{\link{CMBDataFrame}} using
#' the \code{\link{CMBDataFrame}} function (see examples).
#'
#' @examples
#' ## Download SMICA with \code{nside = 1024}
#' ## and save in working directory
#' ## as "CMB_map_smica1024.fits"
#' # downloadCMBMap(foreground = "smica", nside = 1024)
#' ## Load the downloaded map into a CMBDataFrame
#' # sky <- CMBDataFrame("CMB_map_smica1024.fits")
#'
#' ## Download SMICA with Nside=2048 and save in the working directory
#' ## as "CMB_map_smica2048.fits"
#' # downloadCMBMap(foreground = "smica", nside = 2048)
#'
#' ## Download COMMANDER with Nside=1024 and save in a specified folder,
#' ## for example,
#' # dest <- "CMB_map_commander1024.fits"
#' # downloadCMBMap(foreground = "commander", nside = 1024, destfile = dest)
#
#' @keywords CMB Maps
#' @references Planck Public Data Release 2 Maps
#' \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}
#' @references Other fits maps can also be downloaded
#' using the general command \code{\link{download.file}}.
#'
#' @export
downloadCMBMap <- function(foreground = "smica", nside = 1024, destfile){

  webpath <- paste0("http://irsa.ipac.caltech.edu/",
                    "data/Planck/release_2/all-sky-maps/",
                    "maps/component-maps/cmb/")
  prefix <- "COM_CMB_IQU-"
  suffix <- "full.fits"
  ns1024 <- "_1024_R2.02_"
  ns2048 <- "_2048_R2.01_"
  foregrounds <- c("commander","nilc","sevem","smica")

  if ( nside == 1024 )
  {
    ns <- ns1024
  }
  else if ( nside == 2048 )
  {
    ns <- ns2048
  }
  else
  {
    stop(paste0("Please specify nside parameter as ",
                "either nside = 1024 or nside = 2048"))
  }

  foreground <- tolower(foreground)
  if ( !(foreground %in% foregrounds) )
  {
    stop(paste0("foreground parameter was not one of: ",
                paste0(foregrounds, collapse = ", ")))
  }

  url <- paste0(webpath, prefix, foreground, ns, suffix)

  if ( missing(destfile) ) {
    destfile <- paste0("CMB_map_", foreground, nside, ".fits")
  }

  utils::download.file(url, destfile, mode = "wb")
}














#' Download CMB Power Spectra from Planck Legacy Archive.
#'
#' The function \code{downloadCMBPS} downloads
#' CMB power spectra components from
#'  \url{http://pla.esac.esa.int/pla/#cosmology}.
#'
#'
#' \code{link = 1}: Best-fit LCDM CMB power spectra
#' from the baseline Planck
#'  TT, TE, EE+lowE+lensing (2 <= ell <= 2508).
#'
#' \code{link = 2}: Baseline high-ell
#' Planck TT power spectra (2 <= ell <= 2508).
#'
#' \code{link = 3}: Baseline high-ell
#' Planck EE power spectra (2 <= ell <= 1996).
#'
#' \code{link = 4}: Baseline high-ell
#' Planck TE power spectra (2 <= ell <= 1996).
#'
#' \code{link = 5}: Low-ell
#' Planck EB power spectra (2 <= ell <= 29).
#'
#' \code{link = 6}: Low-ell
#' Planck BB power spectra (2 <= ell <= 29).
#'
#' @param link  The link code (an integer from 1 to 6) for the URL to
#' download the file. See code details in this help file.
#' @param save A boolean indicating whether to save or not
#' (since the downloaded
#' data is returned anyway).
#' @param destfile  A character string with the file name for the downloaded
#' file  to be saved. Tilde-expansion is performed.
#'
#' @return The Data Frame with CMB Power Spectra and,
#' if \code{save = TRUE} a txt file is saved in \code{destfile}
#'
#' @examples
#' ## Download the Low-ell Planck BB power spectra (2 <= ell <= 29) and
#' ## save it to C:/PW.txt
#' # downloadCMBPS(link=6, destfile="C:/PW.txt")
#'
#' ## Download the Best-fit LCDM CMB power spectra
#' ## and plot it
#' # CMBPS <- downloadCMBPS(link=1, save = FALSE)
#' # plot(CMBPS$L,CMBPS$TT, type="o",col="red",cex=0.3,
#' #      main="CMB Angular Power Spectra",xlab=expression(l),
#' #      ylab=expression(paste(D[l],"(",mu,K^2,")")))
#'
#' @keywords CMB Power Spectra
#'
#' @references Planck Legacy Archive
#' \url{http://pla.esac.esa.int/pla/#cosmology}
#'
#' @export
downloadCMBPS <- function(link = 1, destfile, save = TRUE) {

  webpath <- "http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID="
  prefix <- "COM_PowerSpect_CMB"
  types <- c("-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory",
            "-TT-full", "-EE-full", "-TE-full",
            "-low-ell-EB-full", "-low-ell-BB-full")
  suffix <- "_R3.01.txt"
  items <- list(c("L","TT","TE","EE","BB","PP"),
                c("L","DL","Minus_dDL","Plus_dDL"),
                c("L","DL","dDL"))

  type <- types[link]

  target <- paste0(webpath, prefix, type, suffix)

  if (missing(destfile)) {
    destfile <- paste0(prefix, type, suffix)
  }

  i <- switch(link, 1, 2, 2, 2, 3, 3)
  CMB_PowerSpectra <- utils::read.table(target, quote = "\"",
                                 col.names = items[[i]])

  if ( save ) utils::write.table(CMB_PowerSpectra, file = destfile)

  invisible(CMB_PowerSpectra)
}
