#' pix2vec
#'
#' Compute the Cartesian cooridinates of the HEALPix points
#' at indices \code{spix} at \code{Nside} in the specified
#' ordering scheme.
#'
#' @param Nside is the Nside of the HEALPix points.
#' @param order specifies the ordering scheme used for the input
#' pixels \code{spix}. If \code{order = "nested"} then the input
#' pixels are converted to ring ordering. If \code{order = "ring"}
#' then no conversion is done.
#' @param spix is a vector of one or more pixel indices whose ordering
#' scheme is given by \code{order} and whose Nside parameter is given
#' by \code{Nside}. If spix = 0 then all pixel indices at
#' \code{Nside} will be used.
#'
#' @return Output is a data.frame with \code{length(spix)} rows, or
#' \eqn{12*Nside^2} rows if \code{spix = 0},
#' and 3 columns of Cartesian coordinates (x,y,z) of HEALPix points
#' at \code{Nside} in ring ordering scheme.
#'
#' @examples
#' # Taking all HEALPix points at Nside = 8, using ring order.
#' pix2vec(Nside = 8, order = "ring")
#'
#' # Compute the HEALPix points in Cartesian coordinates at Nside = 8,
#' # in nested ordering scheme, at the indices from spix
#' Pix <- c(1,2,23)
#' pix2vec(Nside = 8, order = "nested", spix = Pix)
#'
#' @export
pix2vec <- function(Nside = 16, order = "ring", spix = 0) {

  # source("nest2ring.R")
  nPix <- 12*Nside^2
  if ( length(spix) == 1 && spix == 0) {
    spix <- 1:nPix
  }

  if (order == "ring") {
    pix <- spix
  } else if (order == "nested") {
    pix <- nest2ring(Nside,spix)
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
  }

  # Equatorial region
  if (any(eqMask)) {
    ip <- hPix1[eqMask] - nCap -1
    iRing <- trunc(ip/nl4) + Nside
    iPhi <- ip %% nl4 + 1

    # fOdd = 1 if iRing+Nside is odd; fOdd = 1/2 if iRing+Nside is even
    fOdd <- 0.5*(1 + ((iRing + Nside) %% 2))
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

  return(data.frame( x = sth*cos(phi),
                     y = sth*sin(phi),
                     z = z))
}
