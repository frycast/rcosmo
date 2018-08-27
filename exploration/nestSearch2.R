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


# We test siblings by colouring the siblings of some pixels at level j
displaySiblings <- function(p, j, plot.j = 5)
{
  rcosmo:::plotHPBoundaries(nside = 2^j, ordering = "nested", nums.col = "red")

  # We do this by plotting great-great-grandchildren of the siblings
  hp <- rcosmo:::HPDataFrame(nside = 2^plot.j,
                    spix = rcosmo:::pixelWindow(j1 = j, j2 = plot.ns,
                                               pix.j1 = siblings(p)))
  plot(hp, size = 3, add = TRUE)
}

# Looks good:
displaySiblings(4, 1)
displaySiblings(4, 2)
displaySiblings(27, 1)
displaySiblings(28, 1)
displaySiblings(28, 2)
displaySiblings(127, 3)


# We have the siblings of p, but
# it remains to find all neighbours of p at resolution j
rcosmo:::plotHPBoundaries(nside = 1, ordering = "nested")

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
indexInBP(2, 1) # index of pixel 2 at j = 1 (within base pixel 1)
indexInBP(5, 1) # index of pixel 5 at j = 1 (within base pixel 2)

# The base pixel to which p belongs is
basePixel <- function(p, j)
{
  floor((p-1) / (4^j)) + 1
}
basePixel(1,1)
basePixel(4,1)
basePixel(5,1)
basePixel(5,2)

# Note that for all p we have
alwaysTrue <- function(p, j) {
  all((basePixel(p, j) - 1)*4^j + indexInBP(p, j) == p)
}
# Check:
j <- 3
all(alwaysTrue(1:(12*2^j^2), j))


# TO DO: Find the binary representation of indexInBP(p).
# this is all 1s or all zeros (except the rightmost digit which is always 1)
# if and only if p borders a base pixel.


