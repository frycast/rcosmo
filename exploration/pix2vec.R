pix2vec <- function(nSide,ord,Pix) {
# pix2vec computes the Cartesian cooridinates of the HEALPix points 
#    at Pix index at nSide in order ord.
# if nargs() < 3, Pix = 1:nPix, nPix = 12*nSide^2. That is, taking all
#    all HEALPix in the ring order.
# if nargs() < 2, ord = "ring"
# 
# INPUTS:
#   nSide  - Nside of the HEALPix points
#
#   ord    - HEALPix order; in the ring order if ord = "ring", in the nested 
#            order if ord = "nest"
#
#   Pix    - set of Pix index
#
# OUTPUTS:
#   hp     - 3 by length(Pix) matrix of Cartesian coordinates (x,y,z)' of HEALPix points
#            at nSide at Pix index in the order ord
  
source("nest2ring.R")

nPix <- 12*nSide^2
if (nargs() < 3) {
  Pix <- 1:nPix
}

if (nargs() < 2) {
  ord <- "ring"
}

if (ord == "nest") {
  pix <- nest2ring(nSide,Pix)
} else {
  pix <- Pix
}
  
hPix <- pix - 1

# index in 1:n
hPix1 <- trunc(hPix + 1);
nl2 <- trunc(2*nSide);
nl4 <- trunc(4*nSide);

# points in each polar cap, equals 0 for nSide = 1
nCap <- trunc(2*nSide*(nSide-1));
fact1 <- 1.5*nSide;
fact2 <- 3.0*nSide^2;

# associate pixels with areas for north and south caps and equatorial
ncMask <- (hPix1 <= nCap)
eqMask <- ((hPix1 <= nl2*(5*nSide+1)) & (!ncMask));
scMask <- (!(eqMask | ncMask));

# initialisation
nh = length(hPix);
z <- matrix(0,1,nh)
phi <- matrix(0,1,nh)
hDelta_Phi <- matrix(0,1,nh)
z_nv <- matrix(0,1,nh)
z_sv <- matrix(0,1,nh)
phi_nv <- matrix(0,1,nh)
phi_sv <- matrix(0,1,nh)

# North polar cap
if (any(ncMask)) {
  hip <- hPix1[ncMask]/2
  fihip <- trunc(hip)
  iRing <- trunc(sqrt(hip-sqrt(fihip))) + 1
  iPhi <- hPix1[ncMask] - 2*iRing*(iRing-1)

  z[ncMask] <- 1-iRing^2/fact2
phi[ncMask] <- (iPhi-0.5)*pi/(2*iRing)

###
}

# Equatorial region
if (any(eqMask)) {
  ip <- (hPix1[eqMask] - nCap -1)
  iRing <- (trunc(ip/nl4) + nSide)
  iPhi <- (ip %% nl4 + 1)

  fOdd <- 0.5*(1 + ((iRing + nSide) %% 2)) # fOdd = 1 if iRing+nSide is odd; fOdd = 1/2 if iRing+nSide is even
  z[eqMask] <- (nl2-iRing)/fact1
  phi[eqMask] <- (iPhi-fOdd)*pi/(2*nSide)
}

# South polar cap
if (any(scMask)) {
  ip <- (nPix - hPix1[scMask] + 1)
  hip <- ip/2
  fihip <- trunc(hip)
  iRing <- (trunc(sqrt(hip-sqrt(fihip))) + 1)
  iPhi <- (4*iRing + 1 - (ip - 2*iRing*(iRing-1)))

z[scMask] <- (-1+iRing^2/fact2)
phi[scMask] <- ((iPhi - 0.5)*pi/(2*iRing))
}

sth = sqrt(1-z)*sqrt(1+z)

## Pixel centers
# transpose to a column vector
sth <- t(sth)
phi <- t(phi) 
hp <- matrix(0,3,nh)
hp[1,] <- sth*cos(phi)
hp[2,] <- sth*sin(phi)
hp[3,] <- z

return(hp)
}