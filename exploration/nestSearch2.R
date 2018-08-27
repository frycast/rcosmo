rcosmo::plotHPBoundaries(nside = 4)
rcosmo::plotHPBoundaries(nside = 2, col = "blue", lwd = 2)


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
displaySiblings(4, 1)

# We have the siblings of p, but
# it remains to find all neighbours of p at resolution j

# We may need siblings at base pixel resolution
baseSiblings <- function(p)
{
  switch(p,
         c(2,3,4,5,6,9),
         c(1,3,4,6,7,10),
         c(1,2,4,7,8,11),
         c(1,2,3,5,8,12),
         c(1,4,6,8,9,12),
         c(1,2,5,7,9,10),
         c(2,3,6,8,10,11),
         c(3,4,5,7,9,10),
         c(1,5,6,10,11,12),
         c(2,6,7,9,11,12),
         c(3,7,8,9,10,12),
         c(4,5,8,9,10,11))
}

# The number of pixels at resolution j,
# within each base pixel, is 4^j

# Let b be the binary representation of the
# index of pixel p, within its base pixel,
# at resolution j. The decimal representation
# of b is given by
indexInBP <- function(p, j)
{
  (p-1) %% 4^j + 1
}
# indexInBP(2, 1) # index of pixel 2 at j = 1 (within base pixel 1)
# indexInBP(5, 1) # index of pixel 5 at j = 1 (within base pixel 2)

# The base pixel to which p belongs is
basePixel <- function(p, j)
{
  floor((p-1) / (4^j)) + 1
}
# basePixel(1,1)
# basePixel(4,1)
# basePixel(5,1)
# basePixel(5,2)

# Note that for all p we have
alwaysTrue <- function(p, j) {
  all((basePixel(p, j) - 1)*4^j + indexInBP(p, j) == p)
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
#' note \code{nside = 2^j}.
#'
#' Pixel p is on the boundary of a base pixel if and only if
#' the binary representation of \eqn{p - 1} has either
#' its even bits or its odd bits (or both) all zeros or all ones.
#'
#' @param p Integer vector. The pixel index (or indices) at resolution j .
#' @param j The resolution within which to specify p.
#'
#' @examples
#' j <- 3
#' spix <- 17:32 #try 1:16
#' on.b <- onBPBoundary(spix,j)
#' rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)
#' displayPixels(boundary.j = j, j = j, spix = spix[on.b],
#'               boundary.col = "gray", incl.labels = spix)
#'
#' @return A logical vector.
#'
onBPBoundary <- function(p,j)
{
  result <- logical(1)
  # The number of bits in the binary representation is 2*j (giving 4^j pixels).
  vapply(p, function(x) {
    bin <- dec2bin(indexInBP(x, j) - 1, digits = 2*j)
    f <- bin2f(bin, j = j)
    return(sum(f$even) %in% c(0,j)
           || sum(f$odd) %in% c(0,j))
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
neighbours <- function(p, j)
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
p <- 11
displayPixels(boundary.j = j, j = j, plot.j = 5,
              spix = neighbours(p, j),
              boundary.col = "gray",
              boundary.lwd = 1,
              incl.labels = neighbours(p, j),
              col = "blue",
              size = 3)
rcosmo::plotHPBoundaries(nside = 1, col = "blue", lwd = 3)
