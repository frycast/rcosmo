#' Download CMB Power Spectra from Planck Legacy Archive.
#'
#' The function \code{downloadCMBPS} downloads CMB power spectra components from
#'  \url{http://pla.esac.esa.int/pla/#cosmology}.
#'
#'
#' \code{link = 1}: Best-fit LCDM CMB power spectra from the baseline Planck
#'  TT, TE, EE+lowE+lensing (2 <= ell <= 2508).
#'
#' \code{link = 2}: Baseline high-ell Planck TT power spectra (2 <= ell <= 2508).
#'
#' \code{link = 3}: Baseline high-ell Planck EE power spectra (2 <= ell <= 1996).
#'
#' \code{link = 4}: Baseline high-ell Planck TE power spectra (2 <= ell <= 1996).
#'
#' \code{link = 5}: Low-ell Planck EB power spectra (2 <= ell <= 29).
#'
#' \code{link = 6}: Low-ell Planck BB power spectra (2 <= ell <= 29).
#'
#' @param link  The URL to download the file
#' @param destfile  A character string with the file name for the downloaded
#' file  to be saved. Tilde-expansion is performed.
#'
#' @return The Data Frame with CMB Power Spectra and a txt file in destfile
#'
#' @examples
#' ## Download the Low-ell Planck BB power spectra (2 <= ell <= 29) and
#' ## save it to C:/PW.txt
#' downloadCMBPS(link=6, destfile="C:/PW.txt")
#'
#' ## Download the Best-fit LCDM CMB power spectra to the working directory
#' ## and plot it
#' CMBPS<- downloadCMBPS(link=1)
#' plot(CMBPS$L,CMBPS$TT, type="o",col="red",cex=0.3,
#'      main="CMB Angular Power Spectra",xlab=expression(l),
#'      ylab=expression(paste(D[l],"(",mu,K^2,")")))
#'
#' @keywords CMB Power Spectra
#'
#' @references Planck Legacy Archive \url{http://pla.esac.esa.int/pla/#cosmology}
#'
#' @export
downloadCMBPS <- function(link=1, destfile){

  items <- switch(link,
                  items1,
                  items2,
                  items3,
                  items4,
                  items5,
                  items6)

  link.str <- switch(link,
                     l1,
                     l2,
                     l3,
                     l4,
                     l5,
                     l6)
  if ( missing(destfile) ) {
    destfile<- switch(link,
                      destfile11,
                      destfile21,
                      destfile31,
                      destfile41,
                      destfile51,
                      destfile61)
  }

  CMB_PowerSpectra <- read.table(link.str, quote="\"", col.names = items)
  write.table(CMB_PowerSpectra, file=destfile)
  invisible(CMB_PowerSpectra)
}

l1 <- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
            "COM_PowerSpect_CMB-base-plikHM-",
            "TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",sep = "")
items1 <- c("L","TT","TE","EE","BB","PP")
destfile11 <- paste(getwd(),"/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",sep = "")

l2<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
           "COM_PowerSpect_CMB-TT-full_R3.01.txt",sep = "")
items2<- c("L","DL","Minus_dDL","Plus_dDL")
destfile21 <- paste(getwd(),"/COM_PowerSpect_CMB-TT-full_R3.01.txt",sep = "")


l3<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
           "COM_PowerSpect_CMB-EE-full_R3.01.txt",sep = "")
items3<- c("L","DL","Minus_dDL","Plus_dDL")
destfile31 <- paste(getwd(),"/COM_PowerSpect_CMB-EE-full_R3.01.txt",sep = "")


l4<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
           "COM_PowerSpect_CMB-TE-full_R3.01.txt",sep = "")
items4<- c("L","DL","Minus_dDL","Plus_dDL")
destfile41 <- paste(getwd(),"/COM_PowerSpect_CMB-TE-full_R3.01.txt",sep = "")

l5<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
           "COM_PowerSpect_CMB-low-ell-EB-full_R3.01.txt",sep = "")
items5<- c("L","DL","dDL")
destfile51 <- paste(getwd(),"/COM_PowerSpect_CMB-low-ell-EB-full_R3.01.txt",sep = "")

l6<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
           "COM_PowerSpect_CMB-low-ell-BB-full_R3.01.txt",sep = "")
items6<- c("L","DL","dDL")
destfile61 <- paste(getwd(),"/COM_PowerSpect_CMB-low-ell-BB-full_R3.01.txt",sep = "")
