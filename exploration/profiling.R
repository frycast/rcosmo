
# This script is sourced by profilingCall.R
library(rcosmo)
sky <- CMBDataFrame(nside = 32, coords = "cartesian", ordering = "nested")
sky.s <- CMBDataFrame(sky, sample.size = 1000)
hpdf <- HPDataFrame(sky.s, auto.spix = TRUE)


nside <- 32
pix1 <- pix <- apply(sky.s[,c("x","y","z")], MARGIN = 1, nestSearch,
      nside = nside, index.only = TRUE)


pix2 <- nestSearch(sky.s[,c("x","y","z")], nside = nside,
           index.only = TRUE)






