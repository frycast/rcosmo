#' pix2coords
#'
#' convert HEALPix pixel indices to cartesian or spherical coordinates
#'
#' @param nside the nside parameter
#' @param coords 'cartesian' or 'spherical' coordinates
#' @param ordering 'ring' or 'nested' ordering
#' @param spix optional integer or vector of sample pixel indices
#'
#' @return a data.frame with columns 'x', 'y', 'z' (cartesian) or
#' 'theta', 'phi' (spherical)
#'
#' @export
pix2coords <- function(nside, coords = "cartesian", ordering = "nested", spix)
{
  if ( coords != "cartesian" && coords != "spherical")
  {
    stop("coords must be 'cartesian' or 'spherical'")
  }

  if ( ordering != "ring" && ordering != "nested")
  {
    stop("ordering must be 'ring' or 'nested'")
  }

  cart <- (coords == "cartesian")
  nest <- (ordering == "nested")

  if (missing(spix))
  {
    spix <- NULL
  }

  p2c <- pix2coords_internal(nside = nside, nested = nest,
                             spix = spix, cartesian = cart)

  if ( cart )
  {
    return(data.frame(x = p2c[,1], y = p2c[,2], z = p2c[,3]))
  }
  else
  {
    return(data.frame(theta = p2c[,1], phi = p2c[,2]))
  }
}



#' Ring to Nest.
#'
#' \code{ring2nest} converts HEALPix pixel indices in the 'ring' ordering scheme to
#' HEALPix pixel indices in the 'nested' ordering scheme.
#'
#' @param nside is the HEALPix nside parameter.
#' @param pix is a vector of HEALPix pixel indices, in the 'ring' ordering scheme.
#'
#' @return the output is a vector of HEALPix pixel indices in the 'nested' ordering scheme.
#'
#' @examples
#' # compute HEALPix indices in the ring order of the set pix given in the nest order at nside
#' nside <- 8
#' pix <-c(1,2,23)
#' ring2nest(nside,pix)
#'
#' @export
ring2nest <- function(nside, pix) {

  ## initialisations
  ipring <- pix - 1

  # number of pix at nside
  PixNum <- 12*nside^2
  if (!(all(ipring>=0) && all(ipring<PixNum))) {
    print("Error in input pix index!")
  }

  jrll <- c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
  jpll <- c(1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7)
  pix2 <- mkxy2pix()
  x2pix <- pix2$x
  y2pix <- pix2$y


  nl2 <- 2*nside
  nl4 <- 4*nside
  #points in each polar cap, equal 0 when nside=1
  ncap <- nl2*(nside-1)

  #preallocation
  irn <- rep(0,length(ipring))
  iphi <- rep(0,length(ipring))
  nr <- rep(0,length(ipring))
  face_num <- rep(0,length(ipring))
  kshift <- rep(0,length(ipring))

  ##north polar cap
  mask <- ipring<ncap
  mRing <-ipring[mask]

  # counted from North pole
  mirn <- round(sqrt((mRing+1)/2))
  miphi <- mRing - 2*mirn*(mirn - 1)
  MringPhil <- correct_ring_phi(1, mirn, miphi)
  mirn <- MringPhil$irn
  miphi <- MringPhil$iphi
  kshift[mask] <- 0
  nr[mask] <- trunc(mirn)
  face_num[mask] <- trunc(miphi/mirn)
  iphi[mask] <- miphi
  irn[mask] <- mirn

  ## equatorial region
  mask <- (ncap <= ipring) & (ipring < PixNum - ncap)
  mRing <- ipring[mask]
  mip <- trunc(mRing - ncap)
  ## counted from North pole
  mirn <- trunc(mip/nl4) + nside;
  miphi <- mip %% nl4
  kshift[mask]  <- (mirn+nside) %% 2
  nr[mask] <- nside
  # in {1, 2*nside +1}
  mire <- mirn - nside + 1
  mirm <-  nl2 + 2 - mire

  # face boundary
  mifm <- trunc((miphi - trunc(mire/2) + nside)/nside)
  mifp <- trunc((miphi - trunc(mirm/2) + nside)/nside)

  mface <- matrix(rep(0,length(mRing)),ncol=1)
  smask <- (mifp == mifm)
  mface[smask] <- mifp[smask] %% 4+4 # faces 4 to 7
  smask <- (mifp < mifm)
  mface[smask] <- trunc(mifp[smask]); # (half-)faces 0 to 3
  smask <- (mifp > mifm)
  mface[smask] <- trunc(mifp[smask]+7); #(half-)faces 8 to 11
  face_num[mask] <- mface
  irn[mask] <- mirn
  iphi[mask] <- miphi

  ## south polar cap
  mask <- ipring >= (PixNum - ncap)
  mRing <- ipring[mask]

  mip <- trunc(PixNum - mRing)
  # counted from South pole
  mirs <- round(sqrt(mip/2))
  miphi <- 2*mirs*(mirs + 1) - mip

  MringPhil <- correct_ring_phi(1, mirs, miphi)
  mirs <- MringPhil$irn
  miphi <- MringPhil$iphi
  kshift[mask] <-0
  nr[mask] <- trunc(mirs)
  irn[mask] <- trunc(nl4 - mirs)
  face_num[mask] <- trunc(miphi/mirs) + 8
  iphi[mask] <- miphi

  ## finds the (x,y) on the face
  irt <- irn  - jrll[face_num+1]*nside + 1          # in {-nside+1,0}
  ipt <- 2*iphi - jpll[face_num+1]*nr - kshift + 1   # in {-nside+1,nside-1}
  mask <- ipt >= nl2 # for the face #4
  ipt[mask] <- ipt[mask] - 8*nside
  ix <-  trunc((ipt - irt )/2)
  iy <- -trunc((ipt + irt )/2)
  scalemlv <- 1
  scale_factor <- 16384
  ipf <- 0
  # for nside in [2^14, 2^20]
  ismax <- 1
  if (nside > 1048576){
    ismax <- 3
  }

  for (i in 0:ismax) {
    # last 7 bits
    ix_low <- ix %% 128
    iy_low <- iy %% 128
    ipf <- trunc(ipf +(x2pix[ix_low+1]+y2pix[iy_low+1])*scalemlv)
    scalemlv <- scalemlv*scale_factor
    # truncate out last 7 bits
    ix  <- trunc(ix/128)
    iy  <- trunc(iy/128)
  }

  ipf <-  trunc(ipf +(x2pix[ix+1]+y2pix[iy+1])*scalemlv)
  # in {0, 12*nside**2 - 1}
  ipring <- trunc(ipf + (face_num*(trunc(PixNum/12))))
  nPix <- ipring + 1

  return(as.integer(nPix))
}


## HELPER FUNCTION 1 FOR ring2nest
correct_ring_phi <- function(location,iring,iphi){
  delta <- rep(0,length(iphi))
  delta[iphi<0] <- 1
  delta[iphi>4*iring] <- -1
  nzMask <- delta!=0
  if (any(nzMask)){
    iring[nzMask] <- iring[nzMask] - location[nzMask]*delta[nzMask]
    iphi[nzMask]  <- iphi[nzMask] + delta[nzMask]*(4*iring[nzMask])
  }
  MringPhil <- list(irn=iring,iphi=iphi)
  return(MringPhil)
}

## HELPER FUNCTION 2 FOR ring2nest
mkxy2pix <- function() {
  # sets the array giving the number of the pixel lying in (x,y)
  # x and y are in {1,128}
  # the pixel number is in {0,128**2-1}
  # if  i-1 = sum_p=0  b_p * 2^p
  # then ix = sum_p=0  b_p * 4^p
  # iy = 2*ix
  # ix + iy in {0, 128**2 -1}

  x2pix <- rep(128,0)
  y2pix <- rep(128,0)
  for (i in 1:128) {
    j<- i-1
    k<- 0
    ip<- 1
    ip <- 1
    while ( !(j==0) ) {
      id <- j %% 2
      j <- trunc(j/2)
      k<- id*ip+k
      ip<- 4*ip
    }
    x2pix[i]<- k
    y2pix[i]<- 2*k
  }

  pix2 <- list(x=as.integer(x2pix),y=as.integer(y2pix))
  return(pix2)
}


#' Create a new data.frame with a given coordinate system
#'
#' This does not affect the original object \code{df}.
#' If you wish to affect \code{df} then use
#' \code{\link{coords<-.data.frame}}
#'
#'
#'@param df a data.frame with columns labelled x, y, z (for cartesian)
#' or theta, phi (for spherical colatitude and longitude respectively)
#'@param new.coords specifies the new coordinate system
#'("spherical" or "cartesian").
#'
#'@return
#' A new data.frame whose coordinates are as specified by
#' \code{new.coords}
#'
#' @examples
#'
#' ## Create df with no coords, then create df2 with cartesian coords
#' df <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
#' df2 <- coords(df, new.coords = "spherical")
#' df2
#' df
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' df
#'
#'@export
coords.data.frame <- function(df, new.coords)
{
  if ( new.coords == "spherical" )
  {
    if ( all(c("theta","phi") %in% names(df)) )
    {
      return(df)
    }

    if ( !all(c("x","y","z") %in% names(df)) )
    {
      stop(paste0("df must have columns labelled x, y, z (for cartesian) ",
                  "or theta, phi (for spherical colatitude and ",
                  "longitude respectively)"))
    }

    x.i <- which(names(df) == "x")
    y.i <- which(names(df) == "y")
    z.i <- which(names(df) == "z")

    crds <- df[,c(x.i, y.i, z.i)]
    crds <- rcosmo:::car2sph(crds)
    other <- df[,-c(x.i, y.i, z.i), drop = FALSE]
    df <- cbind(crds, other)

  } else if ( new.coords == "cartesian" ) {

    if ( all(c("x","y","z") %in% names(df)) )
    {
      return(df)
    }

    if ( !all(c("theta","phi") %in% names(df)) )
    {
      stop(paste0("df must have columns labelled x, y, z (for cartesian) ",
                  "or theta, phi (for spherical colatitude and ",
                  "longitude respectively)"))
    }

    theta.i <- which(names(df) == "theta")
    phi.i <- which(names(df) == "phi")

    crds <- df[,c(theta.i, phi.i)]
    crds <- rcosmo:::sph2car(crds)
    other <- df[,-c(theta.i, phi.i), drop = FALSE]
    df <- cbind(crds, other)
  }

  return(df)
}



#' Assign new coordinate system to a \code{\link{data.frame}}
#'
#' @seealso \code{\link{coords.data.frame}}
#'
#' @examples
#'
#' ## Create df with no coords, then create df2 with cartesian coords
#' df <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
#' df2 <- coords(df, new.coords = "cartesian")
#' df2
#' df
#'
#' ## Change the coords of df directly (to spherical)
#' coords(df) <- "spherical"
#' df
#'
#' @export
`coords<-.data.frame` <- function(df,...,value) {
  return(coords(df, new.coords = value))
}

