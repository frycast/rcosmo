#' Nest to Ring.
#'
#' \code{nest2ringR} converts HEALPix pixel indices in the 'ring' ordering scheme to
#' HEALPix pixel indices in the 'nest' ordering scheme.
#'
#' @param nSide is the HEALPix Nside parameter.
#' @param Pix is a vector or matrix of HEALPix pixel indices, in the 'nest' ordering scheme.
#'
#' @return the output is a vector or matrix of HEALPix pixel indices in the 'ring' ordering scheme.
#'
#' @examples
#' # compute HEALPix indices in the ring order of the set Pix given in the nest order at Nside
#' Nside <- 8
#' Pix <-c(1,2,23)
#' nest2ring(Nside,Pix)
#'
#' @export
nest2ringR <- function(nSide, Pix) {

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
mask <- jr < nSide
jrN <- jr[mask]
nrN <- trunc(jrN)
n_before[mask] <- trunc(2*nrN*(nrN - 1))
kshift[mask] <- 0
nr[mask] <- nrN

# equatorial area
mask <- (jr >= nSide) & (jr <= 3*nSide)
jrE <- jr[mask]
nrE <- trunc(nSide)
n_before[mask] <- trunc(2*nrE*(2*jrE - nrE - 1))
kshift[mask] <- (jrE - nSide) %% 2
nr[mask] <- nrE

# south pole area
mask <- jr > 3*nSide
jrS <- jr[mask]
nrS <- trunc(nl4 - jrS)
n_before[mask] <- trunc(nPix - 2*nrS*(nrS + 1))
kshift[mask] <- 0
nr[mask] <- nrS

# computes the phi coordinate on S^2, in [0,2*pi)
# 'phi' number in the ring in 1:4*nr
jp <- trunc((jpll[face_num+1]*nr + jpt + 1 + kshift)/2)
maskH <- jp > nl4
maskL <- jp < 1
jp[maskH] <- trunc(jp[maskH] - nl4)
jp[maskL] <- trunc(jp[maskL] + nl4)

# index in 0:nPix-1
ipring <- trunc(n_before + jp - 1)

rPix <- ipring + 1

return(rPix)
}






# HELPER FUNCTION ---------------------------------------------------------

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
    while ( !(jpix==0) ) {
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
    pix2y[kpix+1] <- iy
  }

  pix2 <- list(x=pix2x,y=pix2y)

  return(pix2)
}
