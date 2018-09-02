rcosmo::plotHPBoundaries(nside = 4)
rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 2, ordering = "nested")


# We start by finding the pixel h at resolution j-1
# that contains p, where p is specified at resolution j >= 1.
# Pixel h contains pixels h*4 - 1:4 + 1 at resolution j.
# Of these, t = p - p %% 4 + (p %% 4 != 0)*4 is the top pixel.
# Thus h = t/4.
#
# Parent gives the pixel h at resolution j - 1 that contains p,
# where p is specified at resoution j (notice it does not depend on j).
parent <- function(p)
{
  (p - p %% 4 + (p %% 4 != 0)*4)/4
}

children <- function(h)
{
  1:4 + (h-1)*4
}

# Hence, children(parent(p)) gives siblings of p.
# We thus define,
siblings <- function(p) {
  h <- (p - p %% 4 + (p %% 4 != 0)*4)/4
  1:4 + (h-1)*4
}

#' displayPixels
#'
#' Display the pixels spix at resolution j by colouring
#' in the grandchildren of spixat resolution plot.j
#'
#' @param j The resolution that spix are specified at.
#' @param boundary.j The resolution to display boundaries at. If
#' this is missing then boundaries will not be plotted.
#' @param plot.j The resolution to plot grandchildren at
#' @spix Integer vector. The pixel indices to display.
#' @incl.labels Integer vector of pixel indices to label at
#' resolution j.
#' @boundary.col The boundary colour.
#' @boundary.lwd The boundary line width.
#' @col The colour to make the grandchildren.
#' @size The size to make the grandchildren.
#'
#
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

# We have the siblings of p, but
# it remains to find all neighbours of p at resolution j

# We may need siblings at base pixel resolution
baseSiblings <- function(bp)
{
  # order: S,SE,E,NE,N,NW,W,SW
  # corners: S,E,N,W
  switch(bp,
         # north pole
         c(9 ,6,0,2,3,4,0,5),
         c(10,7,0,3,4,1,0,6),
         c(11,8,0,4,1,2,0,7),
         c(12,5,0,1,2,3,0,8),
         # equatorial
         c(0,9 ,6,1,0,4,8,12),
         c(0,10,7,2,0,1,5,9),
         c(0,11,8,3,0,2,6,10),
         c(0,12,5,4,0,3,7,11),
         # south pole
         c(11,10,0,6,1,5,0,12),
         c(12,11,0,7,2,6,0,9),
         c(9 ,12,0,8,3,7,0,10),
         c(10,9 ,0,5,4,8,0,11))
}

# The number of pixels at resolution j,
# within each base pixel, is 4^j

# Let b be the binary representation of the
# index of pixel p, within its base pixel,
# at resolution j. The decimal representation
# of b is given by
p2ibp <- function(p, j) #indexInBP
{
  (p-1) %% 4^j + 1
}
# p2ibp(2, 1) # index of pixel 2 at j = 1 (within base pixel 1)
# p2ibp(5, 1) # index of pixel 5 at j = 1 (within base pixel 2)

# The base pixel to which p belongs is
p2bp <- function(p, j)
{
  floor((p-1) / (4^j)) + 1
}
# p2bp(1,1)
# p2bp(4,1)
# p2bp(5,1)
# p2bp(5,2)

# Note that for all p we have
alwaysTrue <- function(p, j) {
  all((p2bp(p, j) - 1)*4^j + p2ibp(p, j) == p)
}
# Check:
# j <- 3
# all(alwaysTrue(1:(12*2^j^2), j))


bin2dec <- function(x, digits)
{
  pow <- 2^(0:31)[1:digits]
  sum(pow[as.logical(x)])
}

dec2bin = function(number, digits) {
  as.numeric(intToBits(number))[1:digits]
}

bin2f <- function(bin, j)
{
  even.bits <- bin[seq(2,2*j  , by = 2)]
  odd.bits  <- bin[seq(1,2*j-1, by = 2)]

  return(list(even = even.bits, odd = odd.bits))
}

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
#' @param p Integer vector. The pixel index (or indices) at resolution j .
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
onBPBoundary <- function(p,j)
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


## So far we can check if a pixel is on a boundary of a base pixel
# and we can find the siblings of a pixel.
# Now we want to find all neighbours of a pixel p that is not
# on the boundary of a base pixel.

# displaySiblings(26, 3, col = "red", size = 4)

## Convert p-1 to binary, separate even and odd digits,
# convert each to decimal, add {-1,0,1} to each,
# convert each to binary again, recombine to decimal

# Return neighbours of p assuming p does not border a BP boundary
# and only working with index in BP
neighboursInBP <- function(p, j)
{
  f <- bin2f(dec2bin(p-1, digits = 2*j), j = j)

  # increment/decrement each binary representation
  even.dec <- bin2dec(f$even, digits = j)
  odd.dec <- bin2dec(f$odd, digits = j)
  nbrs.f <- expand.grid(even = list(e = f$even,
                                    eu = dec2bin(even.dec + 1, digits = j),
                                    ed = dec2bin(even.dec - 1, digits = j)),
                        odd = list(o = f$odd,
                                   ou = dec2bin(odd.dec + 1, digits = j),
                                   od = dec2bin(odd.dec - 1, digits = j)))

  # convert the cartesian product of the list of evens/odds back to decimal
  unlist(
  mapply(FUN = function(even,odd) {
           bin <- f2bin(list(even = even, odd = odd), j = j)
           p <- bin2dec(bin, digits = 2*j) + 1
           return(p)
         },
         nbrs.f$even, nbrs.f$odd, SIMPLIFY = FALSE, USE.NAMES = FALSE))
}

j <- 2
p <- 6
displayPixels(boundary.j = j, j = j, plot.j = 5,
              spix = neighbours(p, j),
              boundary.col = "gray",
              boundary.lwd = 1,
              incl.labels = neighbours(p, j),
              col = "blue",
              size = 3)
rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)



# Find the pixel index p of a given pixel at ibp in bp
ibp2p <- function(ibp, bp, j)
{
  (bp - 1)*4^j + ibp
}


neighbours <- function(p, j)
{
  # Get the index in BP and the BP
  ibp <- p2ibp(p, j)
  bp <- p2bp

  # Get even and odd binary digits
  f <- bin2f(dec2bin(ibp-1, digits = 2*j), j = j)

  # separately increment/decrement the decimal representation of evens & odds
  even.dec <- bin2dec(f$even, digits = j)
  odd.dec <- bin2dec(f$odd, digits = j)

  #        S, SE, E, NE, N, NW, W, SW,
  ei <- c(-1, -1,-1,  0, 1,  1, 1,  0)
  oi <- c(-1,  0, 1,  1, 1,  0,-1, -1)

  nbrs <- data.frame(even = rep(even.dec,8) + ei,
                     odd  = rep(odd.dec, 8) + oi)

  # This df has all boundary crossing info in anti-clockwise order.
  # TRUE in 2 columns implies a corner pixel. TRUE in either column
  # implies the corresponding row number will give the correct
  # base pixel change when index from output of baseSiblings
  border <- nbrs < 0 | nbrs >= 2^j
  bcross <- apply(border, any, MARGIN = 1)

  # Base pixel targets after crossing border (zero means drop the pixel)
  bpt <- baseSiblings(bp)[bcross]


}
dec2bin(bin2dec(bin2f(dec2bin(p-1,2*j), j)$odd, digits = j) + 1, j)
