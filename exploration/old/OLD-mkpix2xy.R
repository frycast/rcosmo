mkpix2xy <- function() {
# mkpix2xy calculates the vector of x and y in the face from pixel number for in the
# nested order
#
# OUTPUTS:
#  pix2$x  - pix index for x
#
#  pix2$y  - pix index for y
  
nSide <- 1024
pix2x <- matrix(rep(nSide,0),ncol=nSide)
pix2y <- matrix(rep(nSide,0),ncol=nSide)
for (kpix in 0:(nSide-1)) {
  jpix <- kpix
  ix <- 0
  iy <- 0
  # bit position in x and y
  ip <- 1
  while (!(jpix==0)) {
  # bit value in kpix, for ix
  id <- jpix %% 2
  jpix <- trunc(jpix/2)
  ix <- id*ip+ix
  # bit value in kpix, for iy
  id <- jpix %% 2
  jpix <- trunc(jpix/2)
  iy <- id*ip+iy
  # next bit in x and y
  ip <- 2*ip
  }
  # kpix in 0:31
  pix2x[kpix+1] <- ix
  # kpix in 0:31
  print(iy)
  pix2y[kpix+1] <- iy
}

pix2 <- list(x=pix2x,y=pix2y)

return(pix2)
}