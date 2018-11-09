library(rcosmo)

### Looking at the pointOnCircle and nestSearch part
y <- pointOnCircle()
pix_y <- apply(y, MARGIN = 1, nestSearch, Nside = 16, index_only = TRUE)
## ??? pix_y AREN'T EVEN HEALPix INTEGERS ??? ###





###--------------------------------------------------------------###
### MICROBENCHMARK COMPARING nest_search OLD AND NEW ALGORITHMS  ###
## TEST VALUES SET UP ##
target <- c(0,0,1); j2 <- 3; j1 <- 0; pix.j1 <- 0
# Nside at level j2
nside.j2 <- 2^j2
if ( j1 <= 0 || pix.j1 == 0 ) {
  pix.j2 <- 1:(12*nside.j2^2)
} else {
  # number of points at level j2 in each box at level j1
  lev.diff <- 2^((j2-j1)*2)
  # pix indices at level j2
  pix.j2 <- (lev.diff*(pix.j1-1)+1):(lev.diff*pix.j1)
}
# Convert pix.j2 to Cartesian coordinates
xyz.j2 <- pix2vec(Nside = nside.j2, order_ring = FALSE, Pix = pix.j2)

### NEW ALGORITHM ###
t1 <- function()
{
  dists <-  apply(xyz.j2, MARGIN = 1,
                function(point) {
                  sum((point - target)^2)
                } )
  index.min <- which.min(dists)
  pnt0 <- as.numeric(xyz.j2[index.min,])
}

### OLD ALGORITHM ###
t2 <- function()
{
  index.min <- 1
  pnt0 <- as.numeric(xyz.j2[1,])
  for ( i in 1:nrow(xyz.j2) ) {

    pnt2 <- as.numeric(xyz.j2[i,])

    if (sum((pnt2-target)^2) < sum((pnt0-target)^2) ) {
      index.min <- i
      pnt0 <- as.numeric(xyz.j2[index.min,])
    }
  }
}

### BENCHMARKING ###
library(microbenchmark)
microbenchmark(
  t1(), # t1() is about 15 times faster
  t2()
)
###--------------------------------------------------------------###






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

Trdf <- covCMB(df = df, rmin = 0, rmax = 0.1, Nr = 10, Nside = 1024, N_x_vec = 10)

print(" ****************************************** ")
print(" Plot Covariance for CMB of radius")
print(" ****************************************** ")
r <- Trdf$r
Tcov <- Trdf$Tcov
X11()
plot(r,Tcov,type = "l", col = "blue", lwd = 3,
     main = "Covariance for CMB")
points(r,Tcov,type = "p", col = "red", pch = 20, cex = 1.25)
