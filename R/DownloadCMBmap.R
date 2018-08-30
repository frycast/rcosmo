#' Download CMB Maps from Planck Public Data Release.
#'
#' The function \code{downloadCMBmap} download CMB maps from
#' \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}.
#'
#' CMB maps have been produced by the COMMANDER, NILC, SEVEM, and SMICA
#' pipelines, respectively.
#'
#' For each pipeline, the intensity maps are provided at Nside = 2048, at
#' 5 arcmin resolution, and the polarization maps are provided at Nside = 1024
#' at 10 arcmin resolution.
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
#' @param destfile  A character string with the name where the downloaded file
#' is saved. Tilde-expansion is performed.
#' @return The CMB Map Fits File
#' @examples
#' ## Download Commander with Nside=1024
#' downloadCMBmap(link=1, destfile=1) # obtain 'CMB_map_commander1024.fits' and save in the folder
#' ## Download SMICA with Nside=2048
#' downloadCMBmap(link=8, destfile=8) # obtain 'CMB_map_smica2048.fits' and save in the folder
#'
#' @keywords CMB Maps
#' @references Planck Public Data Release 2 Maps
#' \url{http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/matrix_cmb.html}
#' @references \code{\link{download.file}}.
#' @export
#'
downloadCMBmap <- function(link=1,destfile=1){

  url <- switch(link,
                link1,
                link2,
                link3,
                link4,
                link5,
                link6,
                link7,
                link8)
  destfile_location<- switch(destfile,
          destfile1,
          destfile2,
          destfile3,
          destfile4,
          destfile5,
          destfile6,
          destfile7,
          destfile8)
  download.file(url, destfile_location)

}

link1<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-commander_1024_R2.02_full/index.html"
destfile1=paste(getwd(),"/CMB_map_commander1024.fits",sep = "")
link2<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-nilc_1024_R2.02_full/index.html"
destfile2=paste(getwd(),"/CMB_map_nilc1024.fits",sep = "")
link3<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-sevem_1024_R2.02_full/index.html"
destfile3=paste(getwd(),"/CMB_map_sevem1024.fits",sep = "")
link4<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-smica_1024_R2.02_full/index.html"
destfile4=paste(getwd(),"/CMB_map_smica1024.fits",sep = "")
link5<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-commander-field-Int_2048_R2.01_full/index.html"
destfile5=paste(getwd(),"/CMB_map_commander2048.fits",sep = "")
link6<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-nilc-field-Int_2048_R2.01_full/index.html"
destfile6=paste(getwd(),"/CMB_map_nilc2048.fits",sep = "")
link7<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-sevem-field-Int_2048_R2.01_full/index.html"
destfile7=paste(getwd(),"/CMB_map_sevem2048.fits",sep = "")
link8<- "http://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/cmbpreviews/COM_CMB_IQU-smica-field-Int_2048_R2.01_full/index.html"
destfile8=paste(getwd(),"/CMB_map_smica2048.fits",sep = "")

