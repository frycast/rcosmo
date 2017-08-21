#' Pix to Angle
#' 
#'  \code{pix2ang} computes the Cartesian cooridinates of the HEALPix points 
#'    at Pix index at nSide in order ord.
#'    
#' @param  nSide is the Nside of the HEALPix points.
#'
#' @param order_ring TRUE if the HEALPix order is in the ring order; 
#'        FALSE if in the nested order.
#'
#' @param Pix is the set of Pix index.
#'
#' @return the output is a 2 by length(Pix) matrix of spherical coordinates (theta,phi)' of HEALPix points
#'            at nSide at Pix index in the order ord.
#'            
#' @examples
#' # if nargs() < 3, Pix = 1:nPix, nPix = 12*nSide^2. That is, taking all
#'    all HEALPix in the ring order.
#' Nside <- 8
#' pix2ang(Nside,TRUE)
#' 
#' # comput the HEALPix points in spherical coordinates at Nside in the nest order at pix indices from Pix
#' Nside <- 8
#' Pix <- c(1,2,23)
#' pix2ang(Nside,FALSE,Pix)
#' 
#' @export
pix2ang <- function(nSide,order_ring,Pix) {
  
source("nest2ring.R")

nPix <- 12*nSide^2
if (nargs() < 3) {
  Pix <- 1:nPix
}

if (order_ring==TRUE) {
  pix <- Pix
} else {
  pix <- nest2ring(nSide,Pix)
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

## Pixel centers
# transpose to a column vector
z <- t(z)
phi <- t(phi) 
hp <- matrix(0,2,nh)
hp[1,] <- acos(z)
hp[2,] <- phi

return(hp)
}