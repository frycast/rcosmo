#' Finds the closest pixel center to a point
#'
#' Finds the closest HEALPix pixel center to a given \code{target} point,
#' specified in Cartesian coordinates, using an efficient nested search
#' algorithm. HEALPix indices are all assumed to be in the "nested"
#' ordering scheme.
#'
#' @param target A data.frame, matrix or vector of
#' Cartesian (x,y,z) coordinates for the target point. If
#' a data.frame is used then spherical coordinates can
#' be specified with row names theta and phi.
#' @param nside An integer, the target resolution at which the
#' resulting pixels are returned.
#' @param index.only A boolean indicating whether to return only the
#' pixel index (TRUE), or cartesian coordinates as well (FALSE).
#'
#'
#' @return if \code{index.only = TRUE} then the output will be a HEALPix index.
#' If \code{index.only} FALSE then the output is the list containing the HEALPix index
#' and Cartesian coordinate vector of the HEALPix point closest to \code{target}
#' at resolution \code{nside}.
#'
#' @examples
#'
#' ## Find the closest HEALPix pixel center at resolution j=2 for
#' ## the point (0.6,0.8,0)
#'
#' point <- c(0.6,0.8,0)
#' j <- 2
#' cpoint <- nestSearch(point, nside = 2^j)
#'
#' ## Plot the closest pixel center in blue and the point (0.6,0.8,0) in red
#' displayPixels(j, j, plot.j=j, spix=c(cpoint$pix),
#'               size=5, incl.labels =FALSE)
#' rgl::plot3d(point[1], point[2], point[3],
#'             col="red", size = 5, add = TRUE)
#'
#'
#' ## Repeat the above for 4 points in a data.frame
#' points <- data.frame(x = c(1,0,0,0.6),
#'                      y = c(0,1,0,0.8),
#'                      z = c(0,0,1,0))
#' points
#' j <- 2
#' cpoints <- nestSearch(points, nside = 2^j)
#'
#' ## Plot the closest pixel center in blue and the point (0.6,0.8,0) in red
#' displayPixels(j, j, plot.j=j, spix=c(cpoints$pix),
#'               size=5, incl.labels =FALSE)
#' rgl::plot3d(points[,1], points[,2], points[,3],
#'             col="red", size = 5, add = TRUE)
#'
#' @export
nestSearch <- function(target, nside,
                       index.only = FALSE) {

  # # Convert the target to a list where elements are the row vectors
  # if ( is.numeric(target) && !is.matrix(target) ) { target <- list(target)
  # } else if ( is.matrix(target) ) { target <- as.list(as.data.frame(t(target))) }

  if ( is.matrix(target) ) {

    target <- as.list(as.data.frame(t(target)))
  } else if ( is.data.frame(target) ) {

    if ( all(c("theta","phi") %in% names(target)) ) {

      coords(target) <- "cartesian"
      target <- target[,c("x","y","z")]
    }
    target <- as.list(as.data.frame(t(target)))
  } else if ( is.numeric(target) ) {

    target <- list(target)
  } else {

    stop("Target must be data.frame, matrix or numeric vector")
  }

  j = 0:(log2(nside)+1)
  jlen <- length(j)
  tlen <- length(target)
  h <- rep(list(0), tlen)
  for ( i in 1:jlen ) {

    h <- nestSearch_step( target, j2 = j[i], pix.j1 = h)
  }

  # Note h is now one level deeper than the target resolution.
  # Convert h to Cartesian coordinates.
  h.xyz <- lapply(h, function(x) {
           pix2coords_internal(nside = 2^j[jlen],
                               nested = TRUE,
                               cartesian = TRUE,
                               spix = x) })

  # Get the parent of the closest
  result.h <- list()
  for (i in 1:tlen) {
    dots <- h.xyz[[i]] %*% target[[i]]
    min.h <- h[[i]][max.col(t(dots), ties.method = "first")]
    result.h[[i]] <- parent(min.h)
  }

  h <- unlist(result.h)

  if ( !index.only ) {
    xyz <- pix2coords_internal(nside = nside,
                               nested = TRUE,
                               cartesian = TRUE,
                               spix = h)
    return(list(xyz = xyz, pix = h))
  }

  return(h)
}




#' nestSearch_step
#'
#' This function is inteded for use only with nestSearch.
#'
#' @param target is the target point on S^2 in spherical coordinates.
#' @param j2 is the target resolution.
#' @param pix.j1 is the initial pix index,
#' i.e., the j1-level pixel to search in.
#'
#' @return A vector containing the HEALPix pixel index, \code{pix},
#' of the closest HEALPix pixel center to the target point,
#' \code{target}, at resolution j2, and its neighbours.
#'
#'@keywords internal
#'
#' @export
nestSearch_step <- function(target, j2, pix.j1) {


  # Get the 4 kids
  if (is.list(pix.j1)) { kids <-  lapply(pix.j1, children)
    } else { kids <-  list(children(pix.j1)) }


  nside.j2 <- 2^j2
  # Coordinates of the 4 kids
  xyz.j2 <- lapply(kids, function(x) {
             pix2coords_internal(nside = nside.j2,
                                 nested = TRUE,
                                 cartesian = TRUE,
                                 spix = x) })


  tlen <- length(target)
  result <- list()
  for (i in 1:tlen) {
    dots <- xyz.j2[[i]] %*% target[[i]]
    minkid <- kids[[i]][max.col(t(dots), ties.method = "first")]
    result[[i]] <- rcosmo::neighbours(minkid, j2)
  }


  return( result )
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


#' Return index of parent pixel
#'
#' Gives the pixel at resolution j - 1 that contains p,
#' where p is specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
#' @examples
#'
#'  parent(4)
#'  parent(5)
#'
#' @export
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
  if ( any(p > 0) ) {
    1:4 + rep((p-1)*4, each = 4)
  } else { 1:12 }
}

#' Return siblings of pixel
#'
#' The siblings of pixel p are defined as the
#' children of the parent of p. Note this is resolution independent.
#'
#' @param p Pixel index in nested order.
#'
#' @examples
#'
#' siblings(11)
#' siblings(12)
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
#'@export
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
#'@export
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




#'Return neighbouring pixels
#'
#'Return the neighbouring pixels to a given pixel p
#'that is specified at resolution j, in the nested order.
#'
#'@param p Pixel index p at resolution j.
#'@param j The resolution parameter with nside = 2^j.
#'
#'@examples
#'## Return the neighbouring pixels for base pixel 1
#'neighbours(1, 0)
#'
#'## Plot the neighbouring pixels for base pixel 1
#' demoNeighbours <- function(p,j) {
#'   neighbours(p, j)
#'   displayPixels(boundary.j = j, j = j, plot.j = j+3,
#'                 spix = neighbours(p, j),
#'                 boundary.col = "gray",
#'                 boundary.lwd = 1,
#'                 incl.labels = neighbours(p, j),
#'                 col = "blue",
#'                 size = 3)
#'   rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
#' }
#'
#' demoNeighbours(1,2)
#' demoNeighbours(1,0)
#'
#' @export
neighbours <- function(p, j) {

  # Handle case that p is a vector recursively
  if ( length(p) > 1 )
  {
    np <- length(p)
    #result <- matrix(0, nrow = np, ncol = 9)
    result <- list()
    for ( i in 1:np )
    {
      result[[i]] <- neighbours(p[i], j)
    }
    return(result)
  }

  if ( j == 0 ) {

    bs <- baseNeighbours(p)
    return(bs[bs > 0])
  }

  # Get the index in BP and the BP
  ibp <- p2ibp(p, j)
  bp <- p2bp(p, j)

  # Effectively dec2bin
  digits <- 1:(2*j)
  bin <- as.numeric(intToBits(ibp-1))[digits]

  # Used in bin2dec
  pow <- 2^(0:31)[digits]

  even <- bin[c(FALSE, TRUE)]
  odd <- bin[c(TRUE, FALSE)]

  # Get even and odd binary digit information
  se <- sum( even )
  so <- sum( odd )

  # Get the BP border crossing info: target border pixels
  # bdr is the direction of crossing:
  # 0,1,2,3,4,5,6,7,8 = none, S, SE, E, NE, N, NW, W, SW
  bdr <- borderPattern(onBPBoundary(se, so, j))
  target.bp <- rep(bp, 9)
  target.bp[which(bdr != 0)] <- baseNeighbours(bp)[bdr]

  ## Separately increment/decrement the odd/even binary reps
  # Effectively bin2dec
  even.dec <- sum(pow[as.logical(even)])
  odd.dec <- sum(pow[as.logical(odd)])
  # Order: S,SE, E,NE, N,NW, W,SW,self
  ei <- c(-1,-1,-1, 0, 1, 1, 1, 0,   0)
  oi <- c(-1, 0, 1, 1, 1, 0,-1,-1,   0)
  ed <- rep(even.dec,9) + ei
  od <- rep(odd.dec, 9) + oi

  # even <- lapply(ed, dec2bin, digits = j)
  # odd  <- lapply(od, dec2bin, digits = j)

  dig <- 0:(j-1)
  pos <- dig + rep(c(1, 33, 65, 97, 129, 161, 193, 225, 257), each = j)

  ebits <- as.numeric(intToBits(ed))
  obits <- as.numeric(intToBits(od))

  evenm <- matrix(ebits[pos], ncol = j, byrow = TRUE)
  oddm <- matrix(obits[pos], ncol = j, byrow = TRUE)

  # Swap some nbrs N<->S depending on region of bp
  if ( bp %in% c(1,2,3,4) ) {
    # North pole, swap S->N in directions 4, 5 and 6 (NE, N, NW)
    i4 <- which(bdr == 4)
    i5 <- which(bdr == 5)
    i6 <- which(bdr == 6)

    #SW->NW
    if (length(i4)!=0) {
      oddm[i4,]  <- evenm[i4,]
      evenm[i4,] <- rep(1,j)
    }
    #S->N
    if (length(i5) != 0) {
      evenm[i5,] <- rep(1,j)
      oddm[i5,]  <- rep(1,j)
    }
    #SE->NE
    if (length(i6)!=0) {
      evenm[i6,] <- oddm[i6,]
      oddm[i6,]  <- rep(1,j)
    }

  } else if ( bp %in% c(9,10,11,12) ) {
    # South pole, swap N->S in directions 1, 2 and 8 (S, SE, SW)
    i1 <- which(bdr == 1)
    i2 <- which(bdr == 2)
    i8 <- which(bdr == 8)

    #N->S
    if (length(i1) != 0) {
      evenm[i1,] <- rep(0,j)
      oddm[i1,] <- rep(0,j)
    }
    #NW->SW
    if (length(i2) != 0) {
      evenm[i2,] <- oddm[i2,]
      oddm[i2,] <- rep(0,j)
    }
    #NE->SE
    if (length(i8) != 0) {
      oddm[i8,] <- evenm[i8,]
      evenm[i8,] <- rep(0,j)
    }
  }

  # Recombine even and odd
  bin <- matrix(0, ncol = 2*j, nrow = 9)
  # Effectively f2bin
  bin[,c(FALSE, TRUE)] <- evenm
  bin[,c(TRUE, FALSE)] <- oddm
  # Effectively bin2dec
  nbrs <- (bin %*% pow) + 1

  # Convert indices in BP to actual indices p
  np <- ibp2p(nbrs, target.bp, j)
  return(np[np > 0])
}

# Takes a data.frame with column for even binary and column for odd.
# Returns the recombined digits in decimal as vector. This is
# a helper function for neighbours.
# recombineEvenOdd <- function(nbrs, j)
# {
#   unlist(
#     mapply(FUN = function(even,odd) {
#       bin <- f2bin(list(even = even, odd = odd), j = j)
#       p <- bin2dec(bin, digits = 2*j) + 1
#       return(p)
#     },
#     nbrs$even, nbrs$odd, SIMPLIFY = FALSE, USE.NAMES = FALSE))
# }
