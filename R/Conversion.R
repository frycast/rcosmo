#' Ring to Nest.
#'
#' \code{ring2nest} converts HEALPix pixel indices in the 'ring' ordering scheme to
#' HEALPix pixel indices in the 'nest' ordering scheme.
#'
#' @param nside is the HEALPix nside parameter.
#' @param pix is a vector or matrix of HEALPix pixel indices, in the 'ring' ordering scheme.
#'
#' @return the output is a vector or matrix of HEALPix pixel indices in the 'nest' ordering scheme.
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

  jrll <- matrix(c(2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4), ncol=1)
  jpll <- matrix(c(1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7), ncol=1)
  pix2 <- mkxy2pix()
  x2pix <- pix2$x
  y2pix <- pix2$y


  nl2 <- 2*nside
  nl4 <- 4*nside
  #points in each polar cap, equal 0 when nside=1
  ncap <- nl2*(nside-1)

  #preallocation
  irn <- matrix(rep(0,length(ipring)),ncol=1)
  iphi <- matrix(rep(0,length(ipring)),ncol=1)
  nr <- matrix(rep(0,length(ipring)),ncol=1)
  face_num <- matrix(rep(0,length(ipring)),ncol=1)
  kshift <- matrix(rep(0,length(ipring)),ncol=1)

  ##north polar cap
  mask <- ipring<ncap
  mRing <-ipring[mask]

  # counted from North pole
  mirn <- round(sqrt((mRing+1)/2))
  miphi <- mRing - 2*mirn*(mirn - 1)
  MringPhil <- correct_ring_phi(1, mirn, miphi)
  mirn <- MringPhil$irn
  miphi <- $iphi
  kshift[mask] <- 0
  nr[mask] <- trunc(mirn)
  face_num[mask] <- trunc(miphi/mirn)
  iphi[mask] <- miphi
  irn[mask] <- mirn

  ## equatorial region
  mask <- (ncap <= ipring) & (ipring < PixNum - ncap)
  mRing <- ipring[mask]
  mip    <- trunc(mRing - ncap)
  ## counted from North pole
  mirn   <- trunc(mip/nl4) + nside;
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
  mface[smask] = mifp[smask] %% 4+4 # faces 4 to 7
  smask <- (mifp < mifm)
  mface[smask] = trunc(mifp[smask]); # (half-)faces 0 to 3
  smask <- (mifp > mifm)
  mface[smask] = trunc(mifp[smask]+7); #(half-)faces 8 to 11
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

  return(nPix)
}






#
# ######## THIS ONLY WORKS FOR FULL LENGTH PIX
# ##' Convert from "ring" to "nested" ordering scheme
# ##'
# ##' @param nside is the HEALPix nside parameter
# ##' @param pix is the set or subset of pixel indices at nside
# ##'
# ##' @return vector of pixel indices in "nested" ordering scheme
# ring2nest <- function(nside, pix){
#   if (pix != 1:(12*nside^2) )
#   {
#     stop("(development stage) ring2nest only works with pix = 1:(12*nside^2)")
#   }
#   pix[nest2ring(2, pix)] <- sort(pix)
#   return(pix)
# }


#### MEALPix CODE
#
# npix = nSide2nPix(nside);
# assert(all(ipring>=0) && all(ipring<npix));
#
# % Initializations
# [x2pix, y2pix] = mk_xy2pix;
# jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]'; % in unit of nside
# jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]'; % in unit of nside/2
#
# nl2 = 2*nside;
# nl4 = 4*nside;
#
# ncap = nl2*(nside-1); % points in each polar cap, =0 for nside =1
#
# %% find ring number, ring position & face number
#
# % Preallocation
# irn = zeros(size(ipring));
# iphi = zeros(size(ipring));
# nr = zeros(size(ipring));
# face_num = zeros(size(ipring));
# kshift = zeros(size(ipring));
#
# %% north polar cap
# mask = ipring < ncap;
# mRing = ipring(mask);
#
# % counted from North pole
# mirn   = round(sqrt((mRing+1)/2));
# miphi  = mRing - 2*mirn.*(mirn - 1);
# [mirn, miphi] = correct_ring_phi(1, mirn, miphi);
#
# kshift(mask) = 0;
# nr(mask) = fix(mirn);
# face_num(mask) = fix(miphi./mirn);
# iphi(mask) = miphi;
# irn(mask) = mirn;
#
# %% equatorial region
# mask = (ncap <= ipring) & (ipring < npix - ncap);
# mRing = ipring(mask);
#
# mip    = fix(mRing - ncap);
# % counted from North pole
# mirn   = fix(mip/nl4) + nside;
# miphi  = mod(mip,nl4);
# kshift(mask)  = mod(mirn+nside,2);
# nr(mask) = nside;
# % in {1, 2*nside +1}
# mire =  mirn - nside + 1;
# mirm =  nl2 + 2 - mire;
#
# % face boundary
# mifm = fix((miphi - fix(mire/2) + nside)/nside);
# mifp = fix((miphi - fix(mirm/2) + nside)/nside);
#
#
# mface = zeros(size(mRing));
# smask = (mifp == mifm);
# mface(smask) = mod(mifp(smask),4)+4; % faces 4 to 7
#
# smask = (mifp < mifm);
# mface(smask) = fix(mifp(smask)); % (half-)faces 0 to 3
#
# smask = (mifp > mifm);
# mface(smask) = fix(mifp(smask)+7); % (half-)faces 8 to 11
#
# face_num(mask) = mface;
# irn(mask) = mirn;
# iphi(mask) = miphi;
#
# %% south polar cap
# mask = ipring >= (npix - ncap);
# mRing = ipring(mask);
#
# mip    = fix(npix - mRing);
# % counted from South pole
# mirs   = round(sqrt(mip/2));
# miphi  = 2*mirs.*(mirs + 1) - mip;
# [mirs, miphi] = correct_ring_phi(1, mirs, miphi);
# kshift(mask) = 0;
# nr(mask) = fix(mirs);
# irn(mask)   = fix(nl4 - mirs);
# face_num(mask) = fix(miphi./mirs) + 8;
# iphi(mask) = miphi;
#
# %% finds the (x,y) on the face
# irt =   irn  - jrll(face_num+1)*nside + 1;          % in {-nside+1,0}
# ipt = 2*iphi - jpll(face_num+1).*nr - kshift + 1;   % in {-nside+1,nside-1}
#
# mask = ipt >= nl2;    % for the face #4
# ipt(mask) = ipt(mask) - 8*nside;
#
# ix =  fix((ipt - irt )/2);
# iy = -fix((ipt + irt )/2);
#
# scalemlv = 1;
# scale_factor = 16384;
# ipf = 0;
# % for nside in [2^14, 2^20]
# ismax = 1;
# if (nside >  1048576)
#   ismax = 3;
# end;
# for i = 0:ismax
# % last 7 bits
# ix_low = mod(ix, 128);
# iy_low = mod(iy, 128);
#
# ipf = fix(ipf +(x2pix(ix_low+1)+y2pix(iy_low+1))*scalemlv);
# scalemlv = scalemlv*scale_factor;
# % truncate out last 7 bits
# ix  = fix(ix/128);
# iy  = fix(iy/128);
# end
# ipf =  fix(ipf +(x2pix(ix+1)+y2pix(iy+1))*scalemlv);
#
# % in {0, 12*nside**2 - 1}
# ipnest = fix(ipf + (face_num.*(fix(npix/12))));
# return


