#' Convert from RING to NEST ordering scheme
#'
#' This is a place holder.
#'
#' @param CMBDataFrame A CMB data frame containing a full sky CMB map
#' in RING HEALPix order.
#' @return This is a place holder
ring2nest <- function(CMBdataFrame = cmbdf){
  
  # Complete this using nestIndex.cpp
  # The first step is to create the CMBDataFrame class which has
  # an attribute that can be set to RING or NEST and another
  # attribute that can be set to Nside - we need to get this Nside.
}



#' Convert from NEST to RING ordering scheme
#'
#' This is a place holder.
#'
#' @param CMBDataFrame A CMB data frame containing a full sky CMB map
#' in NESTED HEALPix order.
#' @return This is a place holder
nest2ring <- function(CMBdataFrame = cmbdf){
  
  # Complete this using nestIndex.cpp as above
  
  # Sketch of how it works
  Nside <- 2
  sph <- pix2angC(Nside, Nest = TRUE)
  cmbdf <- data.frame(long = sph[,2], 
                      lat = pi/2 - sph[,1], 
                      I = rnorm(12*Nside^2))
  sourceCpp("nestIndex.cpp")
  index <- nestIndexC(Nside)
  newIndex <- index[order(index[,1],index[,2]),]
  new_cmbdf <- cmbdf[newIndex[,3],]
  
  
  sph2 <- pix2angC(2, Nest = FALSE)
  compare_cmbdf <- data.frame(long = sph2[,2], 
                              lat = pi/2 - sph2[,1], 
                              I = rnorm(12*Nside^2))
  # PROBLEM: THIS exposed a flaw with our pix2angC.
  # The NESTED AND RING pixels are not in the same
  # positions.
}