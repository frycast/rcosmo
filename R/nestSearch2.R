


# We specify p at resolution j >= 1 so nside = 2^j >= 2.
# We want all neighbours n of p at resolution j.
p <- 26

# We start by finding the pixel h at resolution j-1
# that contains p.
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


siblings(4)

# plot.ns <- 5
# hp <- HPDataFrame(nside = plot.ns,
#                   spix = pixelWindow(j1 = 2, j2 = plot.ns,
#                                      pix.j1 = siblings(4)))
