#' nest_search
#' 
#' \code{nest_search(tp,j_2,j_1,ind_1)} searches for the HEALPix point closest to tp in level j_2 in the subregion
#'     found at level j_1. It emploits and requests the nested structure of HEALPix.
#'     
#' \code{nest_search(tp,j_2)} searches for the HEALPix point closest to tp in level
#'     j_2 among all HEALPix points in level j_2.
#'
#' @param tp is the target point on S^2 in spherical coordinates.
#'  
#' @param j_2 is the high level.
#'
#' @param j_1 is the low level.
#' 
#' @param ind_1 is the initial pix index at level j_1 (in which box to search).
#' 
#' @return the output is the list of the Cartesian coordinates and the index
#' of HEALPix point closest to the target point \code{tp}. That is,
#'    \code{h$p} is the HEALPix point closest to tp in level j_2 in the subregion
#'           ind_1 found at level j_1. 
#'    \code{h$ind} is the pix index of poitn hp at level j_2.
#'    
#' @example 
#' # search for the HEALPix point (subregion) closest to North pole (0,0,1) at level 3
#' tp <- c(0,0,1)
#' j_2 <- 3
#' nest_search(tp,j_2,-1,plot_points=TRUE)
#' 
#' @export
nest_search <- function(tp = c(0,0,1), j_2 = 7, j_1 = 3, pix_1 = 1) {
  
  
  # Nside at level j_2
  nside_2 <- 2^j_2
  
  if (j_1<=0) {
    N_2 <- 12*nside_2^2
    Pix_2 <- 1:N_2
  }
  else {
    # number of points at level j_2 in each box at level j_1
    lev_diff <- 2^((j_2-j_1)*2)
    # pix index at level j_2
    Pix_2 <- (lev_diff*(pix_1-1)+1):(lev_diff*pix_1)
  }
  
  # points in Cartesian coordinates at level j_2 index by Pix_2
  Pix_2[] <- sapply(Pix_2,as.numeric)
  hp_2 <- pix2vec(Nside = nside_2, order_ring = FALSE, Pix = Pix_2)
  # HEALPix point and its relative index in hp_2 at level j_2 closest to
  # point tp in geodesic
  i_0 <- 1
  pnt0 <- c(hp_2$x[i_0],hp_2$y[i_0],hp_2$z[i_0])
  for (i in 1:dim(hp_2)[1]) {
    pnt2 <- c(hp_2$x[i],hp_2$y[i],hp_2$z[i])
    if (l2norm(pnt2-tp)<l2norm(pnt0-tp)) {
      i_0 <- i
      pnt0 <- c(hp_2$x[i_0],hp_2$y[i_0],hp_2$z[i_0])
    }
  }
  i_0 <- sapply(i_0,as.numeric)
  ind_2 <- i_0
  hp <- c(hp_2$x[ind_2],hp_2$y[ind_2],hp_2$z[ind_2])
  # HEALPix index at level j_2 of point hp
  pix <- Pix_2[ind_2]
  h <- list(p=hp,ind=pix)
  
  
  # plot the closest HEALPix point to tp

  d_h <- l2norm(hp-tp)
  
  # # HEALPix points at level j_2
  # scatter3D(hp_2$x,hp_2$y,hp_2$z,pch=".",col="blue",bty = NULL,
  #          cex = 3, colkey = FALSE, axes=FALSE,box=FALSE,
  #          main=paste("Nest search for level", j_2, ", dist(tp,hp) = ", sprintf("%.2e",d_h)))
  # 
  # # add the target point tp
  # scatter3D(tp[1],tp[2],tp[3],pch="o",col="red",bty = NULL,
  #          cex = 1, colkey = FALSE, axes=FALSE,box=FALSE,add=TRUE)
  # 
  # # add the closest HEALPix point to tp
  # scatter3D(hp[1],hp[2],hp[3],pch="*",col="red",bty = NULL,
  #          cex = 2, colkey = FALSE, axes=FALSE,box=FALSE,add=TRUE)
  
  return(h)
}




# HELPER FUNCTION ---------------------------------------------------------

# Euclidean norm (l2-norm)
l2norm <- function(x) {
  return(sqrt(dot(x,x)))
}
