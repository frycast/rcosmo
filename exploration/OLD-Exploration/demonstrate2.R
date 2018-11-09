## Clean environment and detach all packages
rm(list = ls())
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),
       detach,character.only=TRUE,unload=TRUE)


df <- CMBDataFrame(nside = 1, I = 1:12)
plotHPBoundaries(nside = 1, ordering = "nested")
plotHPBoundaries(nside = 2, ordering = "nested", col = "blue")
plotHPBoundaries(nside = 4, ordering = "nested", col = "red")

##############################################################
######### Using the window generic ###########################
##############################################################

hpdf <- HPDataFrame(nside = 16, I = 1:(12*16^2))
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
plot(hpdf); plot(win1)
hpdf.win <- window(hpdf, new.window = win1)
plot(hpdf.win, col = "yellow", size = 4, add = TRUE)
attributes(hpdf.win)

##############################################################
######### Object of class CMBDat #############################
##############################################################

cmbdat <- CMBDat("../CMB_map_smica1024.fits", mmap = TRUE)
class(cmbdat)
is.CMBDat(cmbdat)
str(cmbdat)

##############################################################
#### investigating minDist, maxDist, geoDist #################
##############################################################

#### minDist was a cpp function but I converted internal
## and used an R wrapper.
## Gives the shortest distance from point to a data.frame
## or a CMBDataFrame (very simple)

cmbdf <- CMBDataFrame(nside = 1, spix = c(1,5,12), ordering = "ring")
plot(cmbdf, hp.boundaries = 1, col = "blue", size = 5)
p <- c(0,0,1)
minDist(cmbdf, p) # no need to have coordinates

hp <- HPDataFrame(nside = 1, I = rep(0,3), spix = c(1,5,12) )
minDist(hp, p) # notice no need to have coordinates

coords(hp) <- "cartesian"
df <- data.frame(x = hp$x, y = hp$y, z = hp$z)
minDist(df, p)

coords(hp) <- "spherical"
df <- data.frame(theta = hp$theta, phi = hp$phi)
minDist(df, p)


## maxDist is a generic function
## maxDist.CMBWindow
## maxDist.CMBDataFrame

win <- CMBWindow(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
maxDist(win)

#### New and improved geoDist
## geoDist is in Geometry.R and can take two data.frames
p1 <- data.frame(x = c(1, 0, 0), y = c(0, 1, 0), z = c(0, 0, 1))
p2 <- data.frame(theta = c(pi/4, pi/4), phi = c(0,pi/2))
geoDist(p1,p2, include.names = TRUE)



##############################################################
## Simple plot that isn't working on Andriy's PC #############
##############################################################

library(rcosmo)
path <- "C:/Users/danie/Downloads/CMB_maps/"
filename <- "COM_CMB_IQU-commander_1024_R2.02_full.fits"
map <- CMBDat(paste0(path,filename), mmap = TRUE)
sky.sample <- CMBDataFrame(map, sample.size = 100000)
plot(sky.sample, back.col = "white")




##############################################################
##### HPDataFrame Examples ###################################
##############################################################

## Specify locations as vectors and use auto.spix
hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
                   nside = 1, auto.spix = TRUE)
class(hp1)
pix(hp1)
plot(hp1, size = 5, hp.boundaries = 1)
plotHPBoundaries(nside = 1, ordering = "nested", col = "gray")
hp1

## Specify locations in data.frame and specify spix manually
d <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
hp2 <- HPDataFrame(d, nside = 2, spix = c(1,2,1))
plot(hp2, size = 5, hp.boundaries = 2)

## Force the above to be HEALPix center coordinates
class(hp2)
ordering(hp2)
nside(hp2)
pix(hp2)
hp2 <- coords(hp2, new.coords = "cartesian", healpix.only = TRUE)
class(hp2)
ordering(hp2)
nside(hp2)
pix(hp2)
plot(hp2, size = 5, hp.boundaries = 2)
hp2 <- ordering(hp2, new.ordering = "ring")
pix(hp2)
plot(hp2, hp.boundaries = 2, size = 5)
plotHPBoundaries(nside = 2, ordering = "ring", col = "gray")

## Do not specify locations (get all pixels at nside)
hp3 <- HPDataFrame(I = rep(0,12), nside = 1)
hp3
plot(hp3, size = 5, hp.boundaries = 1)
print(hp3)

## At low resolution, a few data points can
## occupy a large pixel area, e.g.:
hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
                   nside = 1, auto.spix = TRUE)
pix(hp1)
geoArea(hp1) #1/4 of the surface area
plot(hp1, size = 5, hp.boundaries = 1)


hp1 <- CMBDataFrame(nside = 1, spix = c(1,2,3))
pix(hp1)
geoArea(hp1) # pi = 1/4*(surface area of unit sphere)
plot(hp1, size = 5, hp.boundaries = 1)



