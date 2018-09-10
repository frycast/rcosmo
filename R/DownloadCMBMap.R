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
#' @param pipeline A string naming the foreground separation method pipeline.
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
#' the \code{\link{CMBDat}} function.
#'
#' @examples
#' ## Download SMICA with \code{nside = 1024}
#' ## and save in working directory
#' ## as "CMB_map_smica1024.fits"
#' downloadCMBMap(foreground = "smica", nside = 1024)
#' ## Load the downloaded map into a CMBDataFrame
#' sky <- CMBDat("CMB_map_smica1024.fits")
#'
#' ## Download SMICA with Nside=2048 and save in the working directory
#' ## as "CMB_map_smica2048.fits"
#' downloadCMBMap(foreground = "smica", nside = 2048)
#'
#' ## Download COMMANDER with Nside=1024 and save in a specified folder,
#' ## for example,
#' dest <- "C:/CMB_map_commander1024.fits"
#' downloadCMBMap(foreground = "commander", nside = 1024, destfile = dest)
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

  download.file(url, destfile, mode = "wb")
}

