# test covCMB

source("covCMB.R")

require("rcosmo")
require("pracma")
require("spark")
source("ca2sph.R")
source("pix2vec.R")
source("nestSearch.R")
source("pointOnCircle.R")
source("pix2vec.R")
source("nest_search.R")
source("nest2ring.R")


## load CMB data in at HEALPix points with nested order
fi_fits  <- "CMB_map_smica1024.fits"

df <- CMBDataFrame(CMBData = fi_fits, coords = "HEALPix", ordering = "nested")

Trdf <- covCMB(df = df, rmin = 0, rmax = 0, Nr = 1, Nside = 1024, N_x_vec = 10)

print(" ****************************************** ")
print(" Plot Covariance for CMB of radius")
print(" ****************************************** ")
r <- Trdf$r
Tcov <- Trdf$Tcov
X11()
plot(r,Tcov,type = "l", col = "blue", lwd = 3,
     main = "Covariance for CMB")
points(r,Tcov,type = "p", col = "red", pch = 20, cex = 1.25)
