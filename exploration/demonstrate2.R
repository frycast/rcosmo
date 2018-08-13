## Clearn environment and detach all packages
rm(list = ls())
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),
       detach,character.only=TRUE,unload=TRUE)

## Simple plot
library(rcosmo)
path <- "C:/Users/danie/Downloads/CMB_maps/"
filename <- "COM_CMB_IQU-commander_1024_R2.02_full.fits"
map <- CMBReadFITS(paste0(path,filename), mmap = TRUE)
sky.sample <- CMBDataFrame(map, sample.size = 100000)
plot(sky.sample, back.col = "white")





##############################################################
##### HPDataFrame Examples ###################################
##############################################################

## Specify locations as vectors and use auto.spix
hp1 <- HPDataFrame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1),
                   nside = 1, auto.spix = TRUE)
class(hp)
pix(hp)
plot(hp, size = 5, hp.boundaries = 1)
plotHPBoundaries(nside = 1, ordering = "nested", col = "gray")

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
plot(hp3, size = 5, hp.boundaries = 1)



