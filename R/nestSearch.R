#' Nest Search
#' 
#' The function \code{nestSearch} computes the nearest HEALPix point at given level 
#' to a point on S^2.
#' 
#' @param tp is a row vector of target point in Cartesian coordinates on S^2
#' @param Nside is the Nside for which the HEALPix points are searched for
#' 
#' @return the output is the HEALPix index only if \code{index_only} TRUE,
#' the output is the list of the HEALPix index and Cartesian coordinates of
#' the HEALPix point closest to \code{tp} if \code{index_only} FALSE.
#' 
#' @examples
#' # Find the pix index and Cartesian coordinates of the HEALPix point 
#' # at Nside closest to the target point tp
#' tp <- c(0,0,1)
#' h <- nestSearch(tp,Nside=2048,index_only=FALSE,plot_points=TRUE)
#' cat("Closest HEALPix point to (",tp,") at Nside = 2048 is (",h$p,")")
#' 
#' @export
nestSearch <- function(tp = c(0,0,1), Nside = 1024, index_only = TRUE){
  
  # # load functions
  # source("nest_search.R")

## multifractal measure for CMB data at level i
# search for the HEALPix point (subregion) closest to tp at level j_1
j_1 <- 3;
h_1 <- nest_search(tp=tp,j_2=j_1,j_1=-1,pix_1=0)
hp_1 <- h_1$p
pix_1 <- h_1$ind

# search for the HEALPix point closest to tp in level j_2 in the subregion
# found at level j_1
j_2 <- 7
h_2 <- nest_search(tp=tp,j_2=j_2,j_1=j_1,pix_1=pix_1)
hp_2 <- h_2$p
pix_2 <- h_2$ind

# search for the HEALPix point closest to tp in level j_3 in the subregion
# found at level j_2
j_3 <- log2(Nside)
h_3 <- nest_search(tp=tp,j_2=j_3,j_1=j_2,pix_1=pix_2)
hp_3 <- h_3$p
pix_3 <- h_3$ind

if (index_only==TRUE) {
  return(pix_3)
} 
else {
  h <- list(p=hp_3,ind=pix_3)
}

# # print results
# tp
# 
# hp_1
# pix_1
# 
# hp_2
# pix_2
# 
# hp_3
# pix_3
}