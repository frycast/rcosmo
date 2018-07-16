######################################
######### INTRODUCING HEALPIX ########
######################################
library(rcosmo)

## Visualise base resolution pixals
cmbdf <- CMBDataFrame(nside = 64, ordering = "nested")
w1 <- window(cmbdf, in.pixels = c(1,9))
w2 <- window(cmbdf, in.pixels = c(2,10))
w3 <- window(cmbdf, in.pixels = c(4,12))
plot(w1, col = "blue", back.col = "white", size = 2)
plot(w2, col = "green", add = TRUE, size = 2)
plot(w3, col = "orange", add = TRUE, size = 2)
plotHPBoundaries(nside = 1, ordering = "nested",
                 incl.labels = 1:12, col = "red")

## Add more boundaries
plotHPBoundaries(nside = 16, ordering = "nested",
                 incl.labels = 1:12, col = "red")










######################################
######### IMPORT HEALPIX #############
######################################

#### PREVIOUS WAY TO IMPORT DATA ####
library(FITSio)

## DON'T RUN THIS, IT TAKES OVER 40 MINUTES, AND CRASHES!
data <- readFITS("../CMB_map_smica2048.fits")


#----------------------------------------------------------------

## Use memory mapping to connect to the file
map <- CMBReadFITS("../CMB_map_smica2048.fits", mmap = TRUE)

## Read the full map
sky <- CMBDataFrame(map)
sky

## Read a uniform sample
set.seed(1)
sky.sample <- CMBDataFrame(map, sample.size = 1e6)
plot(sky.sample, size = 2, back.col = "white")




######################################
######### CONVERSIONS    #############
######################################


## Add coordinates
coords(sky.sample) <- "spherical"
sky.sample


coords(sky.sample) <- "cartesian"
sky.sample


## Switch between ring and nested
ordering(sky.sample) <- "ring"







######################################
###### SUBSETTING AND COMBINING ######
######################################
## Extract an annulus from map
# Exterior of a disc of radius 0.5
dext <- CMBWindow(theta = pi/2, phi = 0, r = 0.5, set.minus = TRUE)
# Interior of a disc of radius 1
disc <- CMBWindow(theta = pi/2, phi = 0, r = 1)
# List of the two discs
annulus <- window(sky.sample, new.window = list(disc, dext))
## Extract a non-convex polygon
theta <-  rep(c(0.1,0.5), 4)
d <- pi/4; phi <- (0:7)*d
swin <- CMBWindow(theta = theta, phi = phi)
star <- window(sky.sample, new.window = swin)
plot(annulus, back.col = "black", size = 3)
plot(dext, lwd = 5); plot(disc, lwd = 5)
plot(sky.small, add = TRUE)
plot(star, add = TRUE, size = 3); plot(swin, lwd = 5)







######################################
######### STATISTICAL METHODS ########
######################################


## Calculate covariance
sky.small <- CMBDataFrame(map, sample.size = 1e5)
Cov <- covCMB(sky.small, max.dist = 0.03, num.bins = 100)
X11()
plot(Cov$dist, Cov$cov/Cov$cov[1],
     main = "Empirical correlation",
     xlab = "Distance",
     ylab = "Correlation")


######################################
###### SUMMARY, GEOMETRY ETC    ######
######################################

summary(annulus)

# area differences

######################################
############ UNIT TESTING ############
######################################

# CTRL + SHIFT + T
























