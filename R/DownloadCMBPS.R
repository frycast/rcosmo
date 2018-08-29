#' Download CMB Power Spectra from Planck Legacy Archive.
#'
#' The function \code{downloadCMBPS} download CMB power spectra components from \url{http://pla.esac.esa.int/pla/#cosmology}.
#'
#' \code{link = 1}: Best-fit LCDM CMB power spectra from the baseline Planck TT,TE,EE+lowE+lensing (2 <= ell <= 2508).
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
#' @return The CMB Power Spectra Data Frame
#' @examples
#' ## Download the Best-fit LCDM CMB power spectra and plot
#' CMBPS<- downloadCMBPS(link=1)
#' plot(CMBPS$L,CMBPS$TT, type="o",col="red",cex=0.3,
#'      main="CMB Angular Power Spectra",xlab=expression(l),ylab=expression(paste(D[l],"(",mu,K^2,")")))
#'
#' @keywords CMB Power Spectra
#' @references Planck Legacy Archive \url{http://pla.esac.esa.int/pla/#cosmology}
#' @export
#'
downloadCMBPS <- function(link=1){

  items <- switch(link,
                  items1,
                  items2,
                  items3,
                  items4,
                  items5,
                  items6)

  link.str <- switch(link,
                     link1,
                     link2,
                     link3,
                     link4,
                     link5,
                     link6)

  CMB_PowerSpectra <- read.table(link.str, quote="\"", col.names = items)
  return(CMB_PowerSpectra)
}

link1 <- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
                "COM_PowerSpect_CMB-base-plikHM-","TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",sep = "")
items1 <- c("L","TT","TE","EE","BB","PP")

link2<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
              "COM_PowerSpect_CMB-TT-full_R3.01.txt",sep = "")
items2<- c("L","DL","Minus_dDL","Plus_dDL")

link3<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
              "COM_PowerSpect_CMB-EE-full_R3.01.txt",sep = "")
items3<- c("L","DL","Minus_dDL","Plus_dDL")

link4<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
               "COM_PowerSpect_CMB-TE-full_R3.01.txt",sep = "")
items4<- c("L","DL","Minus_dDL","Plus_dDL")

link5<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
              "COM_PowerSpect_CMB-low-ell-EB-full_R3.01.txt",sep = "")
items5<- c("L","DL","dDL")

link6<- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
              "COM_PowerSpect_CMB-low-ell-BB-full_R3.01.txt",sep = "")
items6<- c("L","DL","dDL")

