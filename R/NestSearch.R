#' Nested Search
#'
#' Finds the closest HEALPix pixel center to a given \code{target} point,
#' specified in Cartesian coordinates, using an efficient nested search
#' algorithm. HEALPix indices are all assumed to be in the "nested"
#' ordering scheme.
#'
#' @param target A vector of Cartesian coordinates (x,y,z) for the target
#' point.
#' @param nside An integer, the target resolution at which the
#' resulting pixels are returned.
#' @param j A vector of the resolutions to use in each
#' \code{\link{nestSearch_step}}. This is usually best left as the default.
#' Changing this parameter is for experienced users only.
#'
#' @return if \code{index.only = TRUE} then the output will be a HEALPix index.
#' If \code{index.only} FALSE then the output is the list containing the HEALPix index
#' and Cartesian coordinate vector of the HEALPix point closest to \code{target}
#' at resolution \code{nside}.
#'
#' @examples
#'
#' @export
nestSearch <- function(target, nside,
                       index.only = FALSE,
                       j = 0:log2(nside))
{
  j <- c(j[1], j, j[length(j)]+1)
  h <- 0
  for ( i in 2:length(j) )
  {
    h <- rcosmo::nestSearch_step( target, j2 = j[i],
                                   j1 = j[i-1], pix.j1 = h)
  }

  # Note h is now one level deeper than the target resolution.
  # Convert h to Cartesian coordinates.
  h.xyz <- pix2coords_internal(nside = 2^(j[length(j)]),
                                        nested = TRUE,
                                        cartesian = TRUE,
                                        spix = h)[,1:3]

  dists <-  apply(h.xyz, MARGIN = 1,
                  function(point) {
                    sum((point - target)^2)
                  } )
  index.min <- which.min(dists)

  # Find the parent of h, which is at the target resolution
  h <- parent(h[index.min])
  h.xyz <- pix2coords_internal(nside = nside,
                                        nested = TRUE,
                                        cartesian = TRUE,
                                        spix = h)[,1:3]

  return(list(xyz = as.numeric( h.xyz ),
              pix = h))
}




#' nestSearch_step
#'
#' Search for the closest HEALPix pixel to a \code{target} point,
#' where the search is restricted to within HEALPix pixel,
#' \code{pix.j1}, at resolution j1.
#'
#' j1 and j2 are HEALPix resolution parameters, i.e., \eqn{nside = 2^j}.
#'
#' \code{nestSearch_step(target, j2, j1, pix.j1)} searches within the subregion
#' pix.j1, where pix.j1 is a HEALPix pixel index at resolution j1.
#' The return value is the HEALPix point closest to \code{target}, at
#' resolution j2.
#'
#' Setting \code{pix.j1 = 0} (the default) searches for the HEALPix point closest
#' to \code{target} at resolution j2, among all HEALPix points at resolution j1.
#'
#' @param target is the target point on S^2 in spherical coordinates.
#' @param j1 is the lower resolution, with j1 <= j2.
#' @param j2 is the upper resolution.
#' @param pix.j1 is the initial pix index at resolution j1,
#' i.e., the j1-level pixel to search in. If \code{pix.j1 = 0} then
#' all pixels will be searched (slow).
#'
#' @return A vector containing the HEALPix pixel index, \code{pix},
#' of the closest HEALPix pixel center to the target point,
#' \code{target}, at resolution j2, and its neighbours.
#'
#' @examples
#'
#' @export
nestSearch_step <- function(target, j1 = j2, j2, pix.j1 = 0) {

  if ( j1 > j2 ) stop("j1 must be less than or equal to j2")

  # nside at level j2
  nside.j2 <- 2^j2

  spix.j2 <- rcosmo::pixelWindow(j1, j2, pix.j1)

  # Convert spix.j2 to Cartesian coordinates
  xyz.j2 <- pix2coords_internal(nside = nside.j2,
                                         nested = TRUE,
                                         cartesian = TRUE,
                                         spix = spix.j2)[,1:3]

  dists <-  apply(xyz.j2, MARGIN = 1,
                  function(point) {
                    sum((point - target)^2)
                  } )
  index.min <- which.min(dists)

  return( neighbours(spix.j2[index.min], j2) )
}


#' Find high resolution pixels falling in a lower resolution window
#'
#' Find all pixels in a higher resolution that fall within the specified pixel
#' area at a lower resolution. All pixels are assumed to be in nested ordering.
#'
#'@param j1 An integer. The lower resolution, with j1 =< j2.
#'@param j2 An integer. The upper resolution.
#'@param pix.j1 An integer. The pixel index at resolution j1 within which
#'all pixels from resolution j2 will be returned. \code{pix.j1} can
#'also be a vector of non-zero pixel indices.
#'
#'@return All pixels in resolution j2 that fall within the pixel
#'pix.j1 specified at resolution j1
#'
#'@examples
#'
#' pixelWindow(3, 3, 2)
#' pixelWindow(3, 4, 2)
#' pixelWindow(3, 5, 2)
#'
#'@export
pixelWindow <- function(j1, j2, pix.j1)
{
  if ( j2 < 0 || j1 < 0 || pix.j1 < 0 )
  {
    stop("j1, j2, and pix.j1 must all be non-negative")
  }

  if ( any(j2 < j1) ) stop("j2 cannot be less than j1")

  if ( any(pix.j1 > 12*4^(j1)) ) stop("pix.j1 index out of bounds")

  if ( length(pix.j1) == 1 && pix.j1 == 0 ) {
    # pix indices at level j2
    spix.j2 <- 1:(12*2^(2*j2)) #= 12*nside.j2^2

  } else {
  # Number of pixels at level j2 in each pixel from level j1
    lev.diff <- 4^(j2-j1)
    # pix indices at level j2
    spix.j2 <- unlist(mapply(seq, from = (lev.diff*(pix.j1-1)+1),
                             to = (lev.diff*pix.j1),
                             SIMPLIFY = FALSE))
  }
  return(spix.j2)
}



#######################################################################
###############  HELPER FUNCTIONS                     #################
###############                                       #################
#######################################################################


#' parent
#'
#' Gives the pixel at resolution j - 1 that contains p,
#' where p is specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
parent <- function(p)
{
  (p - p %% 4 + (p %% 4 != 0)*4)/4
}

#' Return children of a pixel
#'
#' Gives four pixels at resolution j + 1 that are contained in p,
#' where p is a pixel specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
#'@examples
#'children(11)
#'
#'@export
children <- function(p)
{
  1:4 + (p-1)*4
}

#' siblings
#'
#' The siblings of pixel p are defined as the
#' children of the parent of p. Note this is resolution independent.
#'
#' @param p Pixel index in nested order.
#'
#'@export
siblings <- function(p) {
  h <- (p - p %% 4 + (p %% 4 != 0)*4)/4
  1:4 + (h-1)*4
}



#' Return neighbours of base pixels
#'
#' The base-resolution comprises twelve pixels. \code{baseNeighbours} returns
#' a map from the base pixel index bp to the vector of base pixels
#' that are neighbours of bp, in counterclockwise order of
#' direction: S,SE,E,NE,N,NW,W,SW. The presence of -1 indicates
#' that the corresponding direction is empty.
#'
#' @param bp The base pixel index
#'
#'@examples
#'## Return neighbours of base pixel 1
#'baseNeighbours(1)
#'
#'## There is no base pixel 14, so baseNeighbours returns NULL
#'baseNeighbours(14)
#'
#'@export
baseNeighbours <- function(bp)
{
  # order: S,SE,E,NE,N,NW,W,SW
  # corners: S,E,N,W
  switch(bp,
         # north pole
         c(9 ,6,-1,2,3,4,-1,5),
         c(10,7,-1,3,4,1,-1,6),
         c(11,8,-1,4,1,2,-1,7),
         c(12,5,-1,1,2,3,-1,8),
         # equatorial
         c(-1,9 ,6,1,-1,4,8,12),
         c(-1,10,7,2,-1,1,5,9),
         c(-1,11,8,3,-1,2,6,10),
         c(-1,12,5,4,-1,3,7,11),
         # south pole
         c(11,10,-1,6,1,5,-1,12),
         c(12,11,-1,7,2,6,-1,9),
         c(9 ,12,-1,8,3,7,-1,10),
         c(10,9 ,-1,5,4,8,-1,11))
}


#' a version of onBPBoundary to use with neighbours
#'
#' @param se The sum of even bits, e.g. sum( f$even )
#' @param so The sum of odd bits, e.g. sum( f$odd )
#' @param j The resolution parameter nside = 2^j
#'
#'@keywords internal
#'
onBPBoundary <- function(se, so, j)
{
  b <- 0L
  # south corner
  if ( se == 0 && so == 0 )
  {
    b <- 1L
  }
  # east corner
  else if ( se == 0 && so == j )
  {
    b <- 3L
  }
  # north corner
  else if ( se == j && so == j )
  {
    b <- 5L
  }
  # west corner
  else if ( se == j && so == 0 )
  {
    b <- 7L
  }
  # south east wall
  else if ( se == 0 )
  {
    b <- 2L
  }
  # north east wall
  else if ( so == j )
  {
    b <- 4L
  }
  # north west wall
  else if ( se == j )
  {
    b <- 6L
  }
  # south west wall
  else if ( so == 0 )
  {
    b <- 8L
  }
  return(b)
}



#' Border pattern is resolution independent
#'
#' @param pype is the output of onBPBoundary
#'
#' @return the output is useful as an index
#' to the output of baseNeighbours. It will
#' return the correct neighbouring base
#' pixel in each direction. The
#' 0 indicates to stay in current BP.
#'
#'@keywords internal
#'
borderPattern <- function(ptype)
{
  switch(ptype + 1,
         # S SE  E NE  N NW  W SW
         c(0, 0, 0, 0, 0, 0, 0, 0),
         c(1, 2, 2, 0, 0, 0, 8, 8),
         c(2, 2, 2, 0, 0, 0, 0, 0),
         c(2, 2, 3, 4, 4, 0, 0, 0),
         c(0, 0, 4, 4, 4, 0, 0, 0),
         c(0, 0, 4, 4, 5, 6, 6, 0),
         c(0, 0, 0, 0, 6, 6, 6, 0),
         c(8, 0, 0, 0, 6, 6, 7, 8),
         c(8, 0, 0, 0, 0, 0, 8, 8))
}




#'neighbours
#'
#'Return the neighbouring pixels to a given pixel p
#'that is specified at resolution j, in the nested order.
#'
#'@param p Pixel index p at resolution j.
#'@param j The resolution parameter with nside = 2^j.
#'
#'@examples
#' demoNeighbours <- function(p,j) {
#'   neighbours(p, j)
#'   displayPixels(boundary.j = j, j = j, plot.j = 5,
#'                 spix = neighbours(p, j),
#'                 boundary.col = "gray",
#'                 boundary.lwd = 1,
#'                 incl.labels = neighbours(p, j),
#'                 col = "blue",
#'                 size = 3)
#'   rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
#' }
#'
#' @export
neighbours <- function(p, j)
{
  if ( j == 0 )
  {
    bs <- baseNeighbours(p)
    return(bs[bs > 0])
  }

  # Get the index in BP and the BP
  ibp <- p2ibp(p, j)
  bp <- p2bp(p, j)

  # Get even and odd binary digit information
  f <- bin2f(dec2bin(ibp-1, digits = 2*j), j = j)
  se <- sum( f$even )
  so <- sum( f$odd )

  # Get the BP border crossing info: target border pixels
  # bdr is the direction of crossing:
  # 0,1,2,3,4,5,6,7,8 = none, S, SE, E, NE, N, NW, W, SW
  bdr <- borderPattern(onBPBoundary(se, so, j))
  target.bp <- rep(bp, 9)
  target.bp[which(bdr != 0)] <- baseNeighbours(bp)[bdr]

  # Separately increment/decrement the odd/even binary reps
  even.dec <- bin2dec(f$even, digits = j)
  odd.dec <- bin2dec(f$odd, digits = j)
  # Order: S,SE, E,NE, N,NW, W,SW,self
  ei <- c(-1,-1,-1, 0, 1, 1, 1, 0,   0)
  oi <- c(-1, 0, 1, 1, 1, 0,-1,-1,   0)
  nbrs <- data.frame(even = I(lapply(rep(even.dec,9) + ei,
                                     dec2bin, digits = j)),
                     odd  = I(lapply(rep(odd.dec, 9) + oi,
                                     dec2bin, digits = j)))

  # Swap some nbrs N<->S depending on region of bp
  if ( bp %in% c(1,2,3,4) ) {
    # North pole, swap S->N in directions 4, 5 and 6 (NE, N, NW)
    i4 <- which(bdr == 4)
    i5 <- which(bdr == 5)
    i6 <- which(bdr == 6)

    #SW->NW
    if (length(i4)!=0) {
      nbrs[i4,]$odd  <- nbrs[i4,]$even
      nbrs[i4,]$even <- list(rep(1,j))
    }
    #S->N
    if (length(i5) != 0) {
      nbrs[i5,]$even <- list(rep(1,j))
      nbrs[i5,]$odd  <- list(rep(1,j))
    }
    #SE->NE
    if (length(i6)!=0) {
      nbrs[i6,]$even <- nbrs[i6,]$odd
      nbrs[i6,]$odd  <- list(rep(1,j))
    }

  } else if ( bp %in% c(9,10,11,12) ) {
    # South pole, swap N->S in directions 1, 2 and 8 (S, SE, SW)
    i1 <- which(bdr == 1)
    i2 <- which(bdr == 2)
    i8 <- which(bdr == 8)

    #N->S
    if (length(i1) != 0) {
      nbrs[i1,]$even <- list(rep(0,j))
      nbrs[i1,]$odd <- list(rep(0,j))
    }
    #NW->SW
    if (length(i2) != 0) {
      nbrs[i2,]$even <- nbrs[i2,]$odd
      nbrs[i2,]$odd <- list(rep(0,j))
    }
    #NE->SE
    if (length(i8) != 0) {
      nbrs[i8,]$odd <- nbrs[i8,]$even
      nbrs[i8,]$even <- list(rep(0,j))
    }
  }

  # Recombine even and odd
  nbrs <- recombineEvenOdd(nbrs, j)

  # Convert indices in BP to actual indices p
  np <- ibp2p(nbrs, target.bp, j)
  return(np[np > 0])
}

# Takes a data.frame with column for even binary and column for odd.
# Returns the recombined digits in decimal as vector. This is
# a helper function for neighbours.
recombineEvenOdd <- function(nbrs, j)
{
  unlist(
    mapply(FUN = function(even,odd) {
      bin <- f2bin(list(even = even, odd = odd), j = j)
      p <- bin2dec(bin, digits = 2*j) + 1
      return(p)
    },
    nbrs$even, nbrs$odd, SIMPLIFY = FALSE, USE.NAMES = FALSE))
}
