#' Download CMB Maps from Planck Public Data Release.
#'
#' The function \code{DownloadCMBmap} download CMB maps from \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}.
#'
#' CMB maps have been produced by the COMMANDER, NILC, SEVEM, and SMICA pipelines, respectively.
#'
#' For each pipeline, the intensity maps are provided at Nside = 2048, at 5 arcmin resolution, and the polarization maps are provided at Nside = 1024 at 10 arcmin resolution.
#'
#' \code{link = 1:} CMB Maps produced by Commander with Nside=1024;
#'
#' \code{link = 2:} CMB Maps produced by NILC with Nside=1024;
#'
#' \code{link = 3:} CMB Maps produced by SEVEM with Nside=1024;
#'
#' \code{link = 4:} CMB Maps produced by SMICA with Nside=1024;
#'
#' \code{link = 5:} CMB Maps produced by Commander with Nside=2048;
#'
#' \code{link = 6:} CMB Maps produced by NILC with Nside=2048;
#'
#' \code{link = 7:} CMB Maps produced by SEVEM with Nside=2048;
#'
#' \code{link = 8:} CMB Maps produced by SMICA with Nside=2048;
#'
#' @param link  A character string naming the URL of a resource to be downloaded.
#' @param destfile  A character string with the name where the downloaded file is saved. Tilde-expansion is performed.
#' @return The CMB Map Fits File
#' @examples
#' ## Download Commander with Nside=1024
#' DownloadCMBmap(link=1)
#'
#' @keywords CMB Maps
#' @references Planck Public Data Release 2 Maps \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}
#' @references \code{\link{download.file}}.
#' @export
#'
DownloadCMBmap <- function(link = 1,destfile=getwd()){

  url <- switch(link,
                link1,
                link2,
                link3,
                link4,
                link5,
                link6,
                link7,
                like8)
  download.file(url, destfile)

}

link1<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-commander_1024_R2.02_full/index.html"
link2<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-nilc_1024_R2.02_full/index.html"
link3<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-sevem_1024_R2.02_full/index.html"
link4<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-smica_1024_R2.02_full/index.html"
link5<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-commander-field-Int_2048_R2.01_full/index.html"
link6<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-nilc-field-Int_2048_R2.01_full/index.html"
link7<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-sevem-field-Int_2048_R2.01_full/index.html"
link8<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-smica-field-Int_2048_R2.01_full/index.html"

