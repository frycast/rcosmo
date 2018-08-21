#' Download CMB Power Spectra from Planck Legacy Archive.
#'
#' The function \code{DownloadCMBPS} download CMB power spectra components from \url{http://pla.esac.esa.int/pla/#cosmology}.
#'
#' \code{link = 1:} describe
#'
#' \code{link = 2:} describe
#'
#'
#' @param link  The URL to download the file
#' @param items  Column names to specify the obtained data frame
#' @return The CMB power spectra data frame
#' @examples
#' ## Case 1: Download the product paste ('COM_PowerSpect_CMB-base-plikHM-"
#'              ,"TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",sep = "") with 6 columns ("L","TT","TE","EE","BB","PP")
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                     "COM_PowerSpect_CMB-base-plikHM-","TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt",sep = ""),
#'                     items=c("L","TT","TE","EE","BB","PP"))
#' plot(CMBPS$L,CMBPS$TT, type="o",col="red",cex=0.3,
#'      main="CMB Angular Power Spectra",xlab=expression(l),ylab=expression(paste(D[l],"(",mu,K^2,")")))
#'
#' ## Case 2: Download the product 'COM_PowerSpect_CMB-TT-full_R3.01.txt' with 4 columns ("L","DL","-dDL","+dDL");
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                     "COM_PowerSpect_CMB-TT-full_R3.01.txt",sep = ""),items=c("L","DL","Minus_dDL","Plus_dDL"))
#'
#' ## Case 3: Download the product 'COM_PowerSpect_CMB-EE-full_R3.01.txt' with 4 columns ("L","DL","-dDL","+dDL");
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                     "COM_PowerSpect_CMB-EE-full_R3.01.txt",sep = ""),items=c("L","DL","Minus_dDL","Plus_dDL"))
#'
#' ## Case 4: Download the product 'COM_PowerSpect_CMB-TE-full_R3.01.txt' with 4 columns ("L","DL","-dDL","+dDL");
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                     "COM_PowerSpect_CMB-TE-full_R3.01.txt",sep = ""),items=c("L","DL","Minus_dDL","Plus_dDL"))
#'
#' ## Case 5: Download the product 'COM_PowerSpect_CMB-low-ell-EB-full_R3.01.txt' with 3 columns ("L","DL","dDL");
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                     "COM_PowerSpect_CMB-low-ell-EB-full_R3.01.txt",sep = ""),items=c("L","DL","dDL"))
#'
#' ## Case 6: Download the product 'COM_PowerSpect_CMB-low-ell-BB-full_R3.01.txt' with 3 columns ("L","DL","dDL");
#' CMBPS=DownloadCMBPS(link=paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
#'                    "COM_PowerSpect_CMB-low-ell-BB-full_R3.01.txt",sep = ""),items=c("L","DL","dDL"))
#' @keywords CMB Power Spectra
#' @references Planck Legacy Archive \url{http://pla.esac.esa.int/pla/#cosmology}
#' @export
#'
DownloadCMBPS <- function(link = 1){


  items <- switch(link,
                  1 = items1,
                  2 = items2)

  link.str <- switch(link,
                     1 = link1,
                     2 = link2)

  CMB_PowerSpectra <- read.table(link.str, quote="\"", col.names = items)
  return(CMB_PowerSpectra)
}


link1 <- paste("http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=",
               "COM_PowerSpect_CMB-TT-full_R3.01.txt",sep = "")
items1 <- c("L","DL","Minus_dDL","Plus_dDL")


## Choose default items
# DownloadCMBPS(link = 1)
## usage...
