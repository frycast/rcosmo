

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

#' children
#'
#' Gives the pixels at resolution j + 1 that are contained in p,
#' where p is specified at resoution j (notice it does not depend on j).
#'
#' @param p A pixel index specified in nested order.
#'
children <- function(p)
{
  1:4 + (p-1)*4
}

#' siblings
#'
#' The siblings of pixel p are defined as the
#' children of the parent of p. Note this is resolution independent.
#'
#' @p Pixel index in nested order.
#'
siblings <- function(p) {
  h <- (p - p %% 4 + (p %% 4 != 0)*4)/4
  1:4 + (h-1)*4
}

#' displayPixels
#'
#' Display the pixels spix at resolution j by colouring
#' in the grandchildren of spix at resolution plot.j
#'
#' @param j The resolution that spix are specified at.
#' @param boundary.j The resolution to display boundaries at. If
#' this is missing then boundaries will not be plotted.
#' @param plot.j The resolution to plot grandchildren at
#' @param spix Integer vector. The pixel indices to display.
#' These must be in nested order.
#' @param incl.labels Integer vector of pixel indices to label at
#' resolution j.
#' @param boundary.col The boundary colour.
#' @param boundary.lwd The boundary line width.
#' @param col The colour to make the grandchildren.
#' @param size The size to make the grandchildren.
#'
#'
#'@examples
#'
#' demoNeighbours <- function(p,j) {
#'   neighbours(p, j)
#'   displayPixels(boundary.j = j, j = j, plot.j = 5,
#'                 spix = neighbours(p, j),
#'                 boundary.col = "gray",
#'                 boundary.lwd = 1,
#'                 incl.labels = neighbours(p, j),
#'                 col = "blue",
#'                 size = 3)
#'   rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)
#' }
#'
#'
#'
#'
displayPixels <- function(boundary.j, j, plot.j = 5, spix,
                          boundary.col = "gray",
                          boundary.lwd = 1,
                          incl.labels = 1:(12*4^boundary.j),
                          col = "blue",
                          size = 3)
{
  if ( !missing(boundary.j) ) {
    rcosmo::plotHPBoundaries(nside = 2^boundary.j,
                             ordering = "nested",
                             nums.col = "red",
                             col = boundary.col,
                             lwd = boundary.lwd,
                             incl.labels = incl.labels)
  }

  # We do this by plotting grandchildren of the siblings
  gchild <- rcosmo::pixelWindow(j1 = j,
                                j2 = plot.j,
                                pix.j1 = spix)

  hp <- rcosmo::HPDataFrame(nside = 2^plot.j,
                            spix = gchild)
  plot(hp, add = TRUE, col = col, size = size)
}




# We test siblings by colouring the siblings of some pixels at level j
displaySiblings <- function(p, j, boundary.j = j,
                            plot.j = 5, col = "blue",
                            size = 3,
                            label.siblings = TRUE)
{
  spix <- siblings(p)
  if (label.siblings)
  {
    labels <- siblings(p)
  }
  else
  {
    labels <- 0
  }
  displayPixels(boundary.j = boundary.j, j = j,
                plot.j = plot.j, spix = spix, col = col, size = size,
                incl.labels = labels)
}
# Looks good:
# displaySiblings(4, 1)
# displaySiblings(4, 2)
# displaySiblings(27, 1)
# displaySiblings(28, 1)
# displaySiblings(28, 2)
# displaySiblings(127, 3)
# displaySiblings(4, 1)


#' baseSiblings
#'
#' A map from the base pixel index bp to the vector of base pixels
#' that are neighbours of bp, in counterclockwise order of
#' direction: S,SE,E,NE,N,NW,W,SW
#'
#' @param bp The base pixel index
#'
baseSiblings <- function(bp)
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


#' p2ibp
#'
#' Convert a pixel index p to its index within
#' the base pixel to which p belongs
#'
#' @param p The pixel index at resolution j, in nested order.
#' @param j The resolution parameter nside = 2^j
#'
p2ibp <- function(p, j) #indexInBP
{
  (p-1) %% 4^j + 1
}


#' p2bp
#'
#' The base pixel to which pixel p belongs at resolution j
#'
#' @param p The pixel index at resolution j, in nested order.
#' @param j The resolution parameter nside = 2^j
#'
p2bp <- function(p, j)
{
  floor((p-1) / (4^j)) + 1
}

# Find the pixel index p of a given pixel at ibp in bp
ibp2p <- function(ibp, bp, j)
{
  (bp - 1)*4^j + ibp
}

# # Note that for all p we have
# alwaysTrue <- function(p, j) {
#   all((p2bp(p, j) - 1)*4^j + p2ibp(p, j) == p)
# }
# # Check:
# # j <- 3
# # all(alwaysTrue(1:(12*2^j^2), j))

# Convert binary to decimal
bin2dec <- function(x, digits)
{
  pow <- 2^(0:31)[1:digits]
  sum(pow[as.logical(x)])
}

# Convert decimal to binary
dec2bin = function(number, digits) {
  as.numeric(intToBits(number))[1:digits]
}

# Separate a binary number (e.g., output of dec2bin)
# into its even and odd digits
bin2f <- function(bin, j)
{
  even.bits <- bin[seq(2,2*j  , by = 2)]
  odd.bits  <- bin[seq(1,2*j-1, by = 2)]

  return(list(even = even.bits, odd = odd.bits))
}

# Recombine even and odd digits into a binary number
f2bin <- function(f, j)
{
  bin <- vector(mode = "integer", length = 2*j)
  bin[seq(2,2*j  , by = 2)] <- f$even
  bin[seq(1,2*j-1, by = 2)] <- f$odd
  return(bin)
}


#' onBPBoundary
#'
#' Check if the pixel p at resolution j is on the boundary of a base pixel,
#' note \code{nside = 2^j}. If \code{p} is on the boundary then a number
#' is returned corresponding to the boundary position,
#' in order from 1 to 8 this is S, SE, E, NE, N, NW, W, SW,
#' where S,E,N,W are corners and SE, NE, NW, SW are edges.
#' If \code{p} is not on the boundary then 0 is returned.
#'
#' Pixel p is on the boundary of a base pixel if and only if
#' the binary representation of \eqn{p - 1} has either
#' its even bits or its odd bits (or both) all zeros or all ones.
#'
#' @param p Integer vector. The pixel index (or indices) at resolution j,
#' in nested order.
#' @param j The resolution within which to specify p.
#'
#' @return
#' If \code{p} is on the boundary then a number
#' is returned corresponding to the boundary position,
#' in order from 1 to 8 this is S, SE, E, NE, N, NW, W, SW,
#' where S,E,N,W are corners and SE, NE, NW, SW are edges.
#' If \code{p} is not on the boundary then 0 is returned.
#'
#'@examples
#'
#' # Visualise to see which are on boundary
#' j <- 2
#' spix <- 1:16
#' rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)
#' displayPixels(boundary.j = j, j = j, spix = spix,
#'               boundary.col = "gray", incl.labels = spix)
#'
#' # Corner pixels counter-clockwise S, E, N, W
#' onBPBoundary(c(1,6,16,11),j)
#'
#' # Edge pixels counter-clockwise SE, NE, NW, SW
#' onBPBoundary(c(2,8,16,9),j)
#'
#' # Non-boundary pixels
#' onBPBoundary(c(4,7,10,13),j)
#'
onBPBoundary_old <- function(p,j)
{
  result <- integer(1)
  # The number of bits in the binary representation is 2*j (giving 4^j pixels).
  vapply(p, function(x) {
    bin <- dec2bin(p2ibp(x, j) - 1, digits = 2*j)
    f <- bin2f(bin, j = j)
    se <- sum( f$even )
    so <- sum( f$odd )

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
  }, FUN.VALUE = result)
}



#' a version of onBPBoundary to use with neighbours
#'
#' @param se The sum of even bits, e.g. sum( f$even )
#' @param so The sum of odd bits, e.g. sum( f$odd )
#' @param j The resolution parameter nside = 2^j
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
#' to the output of baseSiblings. It will
#' return the correct neighbouring base
#' pixel in each direction. The
#' 0 indicates to stay in current BP.
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


neighbours <- function(p, j)
{
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
  target.bp[which(bdr != 0)] <- baseSiblings(bp)[bdr]

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


demoNeighbours <- function(p,j) {
  neighbours(p, j)
  displayPixels(boundary.j = j, j = j, plot.j = 5,
                spix = neighbours(p, j),
                boundary.col = "gray",
                boundary.lwd = 1,
                incl.labels = neighbours(p, j),
                col = "blue",
                size = 3)
  rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)
}


# demoNeighbours(1, 2)
# demoNeighbours(6, 2)
# demoNeighbours(16, 2)
# demoNeighbours(11, 2)
# demoNeighbours(172, 3)
