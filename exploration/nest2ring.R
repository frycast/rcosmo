nest2ring <- function(nSide,Pix) {
# nest2ring computes the Pix index in the ring order from the Pix index ipnest
# in the nest order at nSide.
# 
# INPUTS:
#  nSide  - Nside for HEALPix
#
#  ipnest - set or subset of Pix index at nSide
#
# OUTPUTS:
#  ipring - corresponding set of ipnest in the ring order

source("mkpix2xy.R")
  
## initialisations
ipnest <- Pix - 1

# number of pix at nSide
nPix <- 12*nSide^2
if (!(all(ipnest>=0) && all(ipnest<nPix))) {
  print("Error in input Pix index!")
}

jrll <- matrix(c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), ncol=1)
jpll <- matrix(c(1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7), ncol=1)
pix2 <- mkpix2xy()
pix2x <- pix2$x
pix2y <- pix2$y

# number of pixels in a face
npface <- nSide^2
nl4 <- 4*nSide 

## find the face number
# face number in 0:11
face_num <- trunc(ipnest/npface)
ipf <- ipnest %% npface

# finds the x,y on the face (starting from the lowest corner)
# from the pixel number
ix <- 0
iy <- 0
scalemlv <- 1
ismax <- 4

n1 <- 1024 
for (i in 0:ismax) {
  ip_low <- ipf %% n1
  ix <- trunc(ix + scalemlv*pix2x[ip_low+1])
  iy <- trunc(iy + scalemlv*pix2y[ip_low+1])
  scalemlv <- scalemlv*32
  ipf   <- trunc(ipf/n1)
} 
ix <- trunc(ix + scalemlv*pix2x[ipf+1])
iy <- trunc(iy + scalemlv*pix2y[ipf+1])

## transform to (horizontal, vertical) coordinates
# 'vertical' in 0:2*(nSide-1)
jrt <- trunc(ix + iy)
# 'horizontal' in -nSide+1:nSide-1
jpt <- trunc(ix - iy)

## Find z coordinate on S^2
# ring number in 1:4*nSide-1
jr <-  trunc(jrll[face_num+1]*nSide - jrt - 1)

# initialisation
nr <- matrix(rep(0,length(ipnest)),ncol=1)
kshift <- matrix(rep(0,length(ipnest)),ncol=1)
n_before <- matrix(rep(0,length(ipnest)),ncol=1)

# north pole area
mask <- (jr < nSide)
jrN <- jr[mask]
nrN <- trunc(jrN)
n_before[mask] <- trunc(2*nrN*(nrN - 1))
kshift[mask] <- 0
nr[mask] <- nrN

# equatorial area
mask <- ((jr >= nSide) & (jr <= 3*nSide))
jrE <- jr[mask]
nrE <- trunc(nSide)
n_before[mask] <- trunc(2*nrE*(2*jrE - nrE - 1))
kshift[mask] <- (jrE - nSide) %% 2
nr[mask] <- nrE

# south pole area
mask <- (jr > 3*nSide)
jrS <- jr[mask]
nrS <- trunc(nl4 - jrS)
n_before[mask] <- trunc(nPix - 2*nrS*(nrS + 1))
kshift[mask] <- 0
nr[mask] <- nrS

# computes the phi coordinate on S^2, in [0,2*pi)
# 'phi' number in the ring in 1:4*nr
jp <- trunc((jpll[face_num+1]*nr + jpt + 1 + kshift)/2)
maskH <- (jp > nl4)
maskL <- (jp < 1)
jp[maskH] <- trunc(jp[maskH] - nl4)
jp[maskL] <- trunc(jp[maskL] + nl4)

# index in 0:nPix-1
ipring <- trunc(n_before + jp - 1)

rPix <- ipring + 1

return(rPix)
}