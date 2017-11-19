#' pix2vec
#' 
#' \code{pix2vec} computes the Cartesian cooridinates of the HEALPix points 
#'    at Pix index at Nside in the specified order.
#'
#' @param Nside is the Nside of the HEALPix points.
#'
#' @param order_ring TRUE if the HEALPix order is in the ring order; 
#'        FALSE if in the nested order.
#'
#' @param Pix is the set of Pix index.
#'
#' @return output is a 3 by length(Pix) matrix of Cartesian coordinates (x,y,z)' of HEALPix points
#'            at Nside at Pix index in the specified order.
#'
#' @example 
#' # if nargs() < 3, Pix = 1:nPix, nPix = 12*Nside^2. That is, taking all
#'    all HEALPix in the ring order.
#' Nside <- 8
#' pix2vec(Nside,TRUE)
#' 
#' # comput the HEALPix points in Cartesian coordinates at Nside in the nest order at pix indices from Pix
#' Nside <- 8
#' Pix <- c(1,2,23)
#' pix2vec(Nside,FALSE,Pix)
#'
#' @export
pix2vec <- function(Nside = 1024, order_ring = TRUE, Pix = 1) {
  
# source("nest2ring.R")

nPix <- 12*Nside^2
if (nargs() < 3) {
  Pix <- 1:nPix
}

if (order_ring==TRUE) {
  pix <- Pix
} else {
  pix <- nest2ring(Nside,Pix)
}
  
hPix <- pix - 1

# index in 1:n
hPix1 <- trunc(hPix + 1);
nl2 <- trunc(2*Nside);
nl4 <- trunc(4*Nside);

# points in each polar cap, equals 0 for Nside = 1
nCap <- trunc(2*Nside*(Nside-1));
fact1 <- 1.5*Nside;
fact2 <- 3.0*Nside^2;

# associate pixels with areas for north and south caps and equatorial
ncMask <- hPix1 <= nCap
eqMask <- (hPix1 <= nl2*(5*Nside+1)) & (!ncMask)
scMask <- !(eqMask | ncMask)

# initialisation
nh = length(hPix);
z <- rep(0,nh)
phi <- rep(0,nh)
hDelta_Phi <- rep(0,nh)
z_nv <- rep(0,nh)
z_sv <- rep(0,nh)
phi_nv <- rep(0,nh)
phi_sv <- rep(0,nh)

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
  ip <- hPix1[eqMask] - nCap -1
  iRing <- trunc(ip/nl4) + Nside
  iPhi <- ip %% nl4 + 1

  fOdd <- 0.5*(1 + ((iRing + Nside) %% 2)) # fOdd = 1 if iRing+Nside is odd; fOdd = 1/2 if iRing+Nside is even
  z[eqMask] <- (nl2-iRing)/fact1
  phi[eqMask] <- (iPhi-fOdd)*pi/(2*Nside)
}

# South polar cap
if (any(scMask)) {
  ip <- nPix - hPix1[scMask] + 1
  hip <- ip/2
  fihip <- trunc(hip)
  iRing <- trunc(sqrt(hip-sqrt(fihip))) + 1
  iPhi <- 4*iRing + 1 - (ip - 2*iRing*(iRing-1))

z[scMask] <- -1+iRing^2/fact2
phi[scMask] <- (iPhi - 0.5)*pi/(2*iRing)
}

sth <- sqrt(1-z^2)

## Pixel centers
hp <- data.frame( x = sth*cos(phi),
                  y = sth*sin(phi),
                  z = z)

return(hp)
}