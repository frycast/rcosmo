#' Download CMB Maps from Planck Public Data Release.
#'
#' The function \code{downloadCMBMap} downloads CMB maps from
#' \url{https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}.
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
#' @param destfile A required character string with an explicit path and file
#' name for the downloaded file. No default destination is used.
#' Tilde-expansion is performed.
#' @param release An integer. Indicates the Planck map release
#' to download. Currently only 2 and 3 are supported.
#'
#' @return The integer status code from \code{utils::download.file},
#' returned invisibly (0 indicates success). The FITS map is saved to
#' \code{destfile} and can be read with \code{CMBDataFrame}.
#'
#' @examples
#' ## Requires internet access and downloads a large FITS file.
#' if (interactive()) {
#'   dest <- tempfile(fileext = ".fits")
#'   downloadCMBMap(foreground = "smica", nside = 1024, destfile = dest)
#'   sky <- CMBDataFrame(dest, spix = 1:10)
#'   unlink(dest)
#' }
#'
#' @references Planck Public Data Release 2 Maps
#' \url{https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}
#' @references Other fits maps can also be downloaded
#' using the general command \code{\link{download.file}}.
#'
#' @export
downloadCMBMap <- function(
  foreground = "smica", nside = 1024, destfile, release = 2){

  if (missing(destfile) || !is.character(destfile) ||
      length(destfile) != 1L || is.na(destfile) || !nzchar(destfile)) {
    stop("destfile must be an explicit, non-empty file path")
  }

  webpath <- paste0("http://irsa.ipac.caltech.edu/",
                    "data/Planck/release_",release,
                    "/all-sky-maps/",
                    "maps/component-maps/cmb/")
  prefix <- "COM_CMB_IQU-"
  suffix <- "full.fits"
  ns1024 <- "_1024_R2.02_"
  ns2048r2 <- "-field-Int_2048_R2.01_"
  ns2048r3 <- "_2048_R3.00_"
  foregrounds <- c("commander","nilc","sevem","smica")

  if (! release %in% c(2,3) ) {
    stop("downloadCMBMap only supports release=2 or release=3")
  }
  if (release == 3 & nside != 2048) {
    stop("when release=3, downloadCMBMap only supports nside=2048")
  }

  if ( nside == 1024 )
  {
    ns <- ns1024
  }
  else if ( nside == 2048 )
  {
    if (release == 2) {
      ns <- ns2048r2
    }
    else if (release == 3) {
      ns <- ns2048r3
    }
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


  utils::download.file(url, destfile, mode = "wb")
}














#' Download CMB Power Spectra from Planck Legacy Archive.
#'
#' The function \code{downloadCMBPS} downloads
#' CMB power spectra components from
#'  \url{https://www.cosmos.esa.int/web/planck/pla}.
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
#' @param destfile A character string with an explicit path and file name.
#' Required when \code{save = TRUE}; unused when \code{save = FALSE}.
#' No default destination is used. Tilde-expansion is performed.
#'
#' @return A data frame returned invisibly, containing multipoles in column
#' \code{L} and the selected power spectrum and uncertainty columns described
#' above. When \code{save = TRUE}, it is also written to \code{destfile}.
#'
#' @examples
#' ## Requires internet access.
#' if (interactive()) {
#'   dest <- tempfile(fileext = ".txt")
#'   spectrum <- downloadCMBPS(link = 6, destfile = dest)
#'   unlink(dest)
#'   spectrum <- downloadCMBPS(link = 1, save = FALSE)
#'   plot(spectrum$L, spectrum$TT, type = "l")
#' }
#'
#' @references Planck Legacy Archive
#' \url{https://www.cosmos.esa.int/web/planck/pla}
#'
#' @export
downloadCMBPS <- function(link = 1, destfile, save = TRUE) {

  if (save) {
    if (missing(destfile) || !is.character(destfile) ||
        length(destfile) != 1L || is.na(destfile) || !nzchar(destfile)) {
      stop("destfile must be an explicit, non-empty file path")
    }
  }

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


  i <- switch(link, 1, 2, 2, 2, 3, 3)
  CMB_PowerSpectra <- utils::read.table(target, quote = "\"",
                                 col.names = items[[i]])

  if ( save ) utils::write.table(CMB_PowerSpectra, file = destfile)

  invisible(CMB_PowerSpectra)
}
