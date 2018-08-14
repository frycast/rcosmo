# rm(list=ls())
# library(microbenchmark)
# library(sphereplot)
# library(pryr)
# library(Rcpp)
# library(Rcpp)
# library(tidyverse)
library(rcosmo)
# Install mmap from tarball
devtools::install_local("../mmap", force = TRUE)
library(mmap)


### GENERATE DOCUMENTATION
pack <- "rcosmo"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),
             "CMD", "Rd2pdf", shQuote(path)))






# Make a CMBDat object and take a sample data.frame
map <- CMBReadFITS("../CMB_map_smica1024.fits", mmap = TRUE)
s <- sample(1:(12*1024^2), 1000000)
dat.sample <- map$data[sort(s)]

# Get windows to make CMBDataFrames
win <- CMBWindow(x = 1, y = 0, z = 0, r = 0.1)
win2 <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                     theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))
cap <- window(map, new.window = win)
polygon <- window(map, new.window = win2)
plot(cap)
plot(polygon)

sky <- CMBDataFrame(map, sample.size = 1000000)
plot(sky)
plot(cap, add = TRUE)
plot(polygon, add = TRUE)




mmap::munmap(map$data)






## Test a non-CMB map
#http://cade.irap.omp.eu/dokuwiki/doku.php?id=chipass
sky.t <- CMBDataFrame("C:\Users\dfryer\Downloads\CHIPASS_1_1024.fits",
                      coords = "cartesian")




######## Experimenting code for thesis ################
sky <- CMBDataFrame("../CMB_map_smica1024.fits", include.masks = TRUE)
sky.sample <- CMBDataFrame(sky, sample.size = 100000, coords = "cartesian")
win <- CMBWindow(x = 1, y = 0, z = 0, r = 0.1)
cap <- window(sky, new.window = win)
plot(cap, back.col = "white")

disc <- CMBWindow(theta = pi/2, phi = -pi/8, r = 0.2)
polygon <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                     theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))

plot(sky.sample, back.col = "white")
plot(disc, lwd = 3, col = "blue"); plot(polygon, lwd = 3, col = "blue")
summary(polygon)
summary(disc)


d.exterior <- CMBWindow(theta = pi/2, phi = 0, r = 0.5,
                        set.minus = TRUE)

wins <- list(d.exterior, CMBWindow(theta = pi/2, phi = 0, r = 1))

sky.annulus <- window(sky.sample, new.window = wins)
plot(wins[[1]], col = "blue", lwd = 5)
plot(wins[[2]], col = "blue", lwd = 5)
plot(sky.annulus, back.col = "white", add = TRUE)
plot(wins[[1]], col = "blue", size = 5)
plot(wins[[2]], col = "blue", size = 5)

sky.cov <- covCMB(sky.sample, max.dist = 0.03, num.bins = 200)

# Apply the temperature mask
sky.masked <- sky[as.logical(sky$TMASK),]

# Visualise the temperature mask
plot(sky.sample, col = sky.sample$TMASK + 1, size = 3, back.col = "white")

# Calculate var and mean on unmasked sky
mean(sky$I)
var(sky$I)


# Calculate var and mean on masked sky
mean(sky.masked$I)
var(sky.masked$I)


alpha <- -10
A <- pixelArea(sky)
A*sum(sky.masked$I < alpha)
##########################################################

header(sky1)


win <- CMBWindow(x = 1, y = 0, z = 0, r = 0.05)
win2 <- CMBWindow(theta = pi/2, phi = 0, r = 0.05)
win3 <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                  theta = c(pi/2, pi/2, pi/10, pi/10))
plot(win3)

skywin <- window(sky, new.window = list(win,win2))

summary(skywin)
summary(win)
summary(win2)
summary(win3)

a <- summary(sky)
class(a)
a

library(sp)

#################################################################
######   DEMONSTRATE APPLYING MASKS #############################
#################################################################
library(gsl)
N <- 1000

sky <- CMBDataFrame("../CMB_map_smica1024.fits", include.masks = TRUE)
sky.masked <- sky[as.logical(sky$TMASK),]

sky.masked <- CMBDataFrame(sky.masked, coords = "cartesian")
sky.unmasked <- CMBDataFrame(sky, coords = "cartesian")

plot(sky.masked)
cov.masked <- covCMB(sky.masked, max.dist = 0.03, sample.size = 100000, num.bins = N)
cov.unmasked <- covCMB(sky.unmasked, max.dist = 0.03, sample.size = 100000, num.bins = N)

X11(); plot(cov.masked$dist, cov.masked$cov, main = "Masked")
X11(); plot(cov.unmasked$dist, cov.unmasked$cov, main = "Unmasked")

#We use rev to order cosines of distances from -1 to 1
C_est <- rev(cov.masked$cov)*1e9 # micro kelvin
s <- rev(cos(cov.masked$dist))
plot(s, C_est)

#Components of angular power spectrum
n<-1000
## See: http://www.gnu.org/software/gsl/doc/html/specfunc.html#legendre-functions-and-spherical-harmonics
#fl_func <- function(i) {4*pi*sum(C_est*legendre_Plm(i, 0, s, give=FALSE, strict=TRUE))/N} # This was the wrong function for Pl
fl_func <- function(l) {4*pi*sum(C_est/cov.masked$cov[1]*legendre_Pl(l, s))/N}
fl<- sapply(1:n, fl_func)

#Plot angular power spectrum
X11()
plot(fl*(1:n)*(1:n+1)/(2*pi), type="l", col=grey(.5))

plot(200:length(fl),fl[200:length(fl)], type="l", col=grey(.5))

#################################################################
######## DEMONSTRATE USE OF equiareal PAREMETER IN covCMB ######
################################################################

# Using data from above

# This was really slow for some reason or maybe it just froze
cov.ea <- covCMB(sky.unmasked, equiareal= TRUE, sample.size = 100000, num.bins = 1000)
cov.not.ea <- covCMB(sky.unmasked, equiareal = FALSE, sample.size = 100000, num.bins = 1000)



#################################################################
########## DEMONSTRATE nside = 2048 coordinates workaround ######
#################################################################

## When using CMB data with nside = 2048 many machines will not
## have enough RAM to assign cartesian coordinates, so the
## window function will crash since it implicitly assigns
## cartesian coordinates. The work-around is to use the
## in.pixels parameter of the window function after visually
## determining which pixels the window intersects, as in the
## example below.


#################################################################
########## DEMONSTRATE nest_search ##############################
#################################################################

# Closest point in level j2, when searching within
# pixel pix.j1 at level j1.
j1 <- 0
j2 <- 1
next.pix.j1 <- nestSearch_step(c(0,0,1), j2 = j2, j1 = j1,
                               pix.j1 = 1, demo.plot = TRUE)
plotHPBoundaries(2^j1, ordering = "nested")

# Deeper...
j1 <- 1
j2 <- 2
next.pix.j1 <- nestSearch_step(c(0,0,1), j2 = j2, j1 = j1,
                               pix.j1 = next.pix.j1$pix,
                               demo.plot = TRUE)

# Deeper...
j1 <- 2
j2 <- 3
nestSearch_step(c(0,0,1), j2 = j2, j1 = j1,
                pix.j1 = next.pix.j1$pix,
                demo.plot = TRUE)


# Deepest...
nestSearch(c(0,0,1), 32, demo.plot = TRUE)


# Also works for points outside the unit sphere and any nside
nestSearch(c(0,2,2), nside = 16224)


library(microbenchmark)
microbenchmark(
  nestSearch(c(0,0,1), 2048, j = c(1,2,3,4,5,6,7,8,9,10,11)),
  nestSearch(c(0,0,1), 2048, j = c(3,7,11)) )


# Problem with nested search j_0 = 0
label <- function(x, bpix)
{
  dists <-  apply(bpix, MARGIN = 1,
                  function(bp) {
                    sum((bp - x)^2)
                  } )
  return(which.min(dists))
}

## Fail example for j_0 = 0 (colour boundaries at nside = 1)
base1 <- CMBDataFrame(nside = 1, coords = "cartesian")[,c("x","y","z")]
cmbdf1 <- window(CMBDataFrame(nside = 128, coords = "cartesian"),
                 in.pixels = c(1,4,5), in.pixels.res = 0)
cols1 <- apply(cmbdf1[,c("x","y","z")], MARGIN = 1, label, bpix = base1)
cmbdf1[,c("x","y","z")] <- cmbdf1[,c("x","y","z")]*0.99

cols1[cols1 == 1] <- 6
plot(cmbdf1, col = cols1, size = 5, back.col = "white")
plotHPBoundaries(nside = 1, lwd = 7, col = "black", ordering = "nested",
                 nums.size = 1.7)

theta <- 0.95 #phi = 0
nestSearch(c(sin(theta), 0, cos(theta)), 32, j = 0:log2(32), demo.plot = TRUE)

## Fail example for j_0 = 3 (colour boundaries at nside = 8)
base2 <- CMBDataFrame(nside = 8, coords = "cartesian")[,c("x","y","z")]
cmbdf2 <- window(CMBDataFrame(nside = 1024, coords = "cartesian"), in.pixels = c(59,60,246,248), in.pixels.res = 3)
cols2 <- apply(cmbdf2[,c("x","y","z")], MARGIN = 1, label, bpix = base2)
cmbdf2[,c("x","y","z")] <- cmbdf2[,c("x","y","z")]*0.99

plot(cmbdf2, col = cols2, size = 5, back.col = "white")
plotHPBoundaries(nside = 8, lwd = 5, col = "black",
                 incl.labels = c(59,60,246,248), ordering = "nested")

theta <- 0.35
nestSearch(c(sin(theta), 0, cos(theta)), 64, j = 3:log2(64), demo.plot = TRUE)


#################################################################
########## DEMONSTRATE TMASK ####################################
#################################################################

## Test include.masks = TRUE
sky.m <- CMBDataFrame("../CMB_map_smica1024.fits", include.masks = TRUE,
                      coords = "cartesian")
plot(sky.m, col = sky.m$TMASK + 1, sample.size = 50000)



#################################################################
##### Experiment with functions for cheat sheet #################
#################################################################


sky <- CMBDataFrame("C:/Users/dfryer/Downloads/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits")
ff <- CMBReadFITS("C:/Users/dfryer/Downloads/COM_CMB_IQU-smica-field-Int_2048_R2.01_full.fits")

#plot full sky after random sample
plot(sky, sample.size = 100000)

#add detailed window and display clear window boundaries
plot(win2, add = TRUE); plot(poly); plot(disc)

#visualise healpix boundaries and pixel numbers
plotHPBoundaries(nside = 1, col = "blue", ordering = "nested")

poly <- CMBWindow(theta = c(1,2,2), phi = c(0,0,1))
disc <- CMBWindow(x = 0, y = 0, z = 1, r = 0.5)
win1 <- CMBDataFrame(sky, win = poly)
win2 <- CMBDataFrame(sky, win = list(poly, disc))

win2 <- CMBDataFrame(win2, yourVar = rep(1, nrow(win2)))

m1 <- matrix(c(0,1,0,1,0,0), nrow = 2, byrow = TRUE)
m2 <- data.frame(theta = c(pi/2,pi/2), phi = c(0,pi/4))
geoDist(m1,m2)


#################################################################
##### Change RGL viewpoint and create rotating view option ######
#################################################################

library(rcosmo)
library(rgl)
a <- CMBDataFrame(nside = 32, ordering = "nested", coords = "spherical")
plot(a)
rgl.viewpoint(theta = 0, phi = 0)


##################################################################
################ TESTING `[` and coords() ########################
##################################################################


# NESTED ORDERING
a1 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "spherical")
b1 <- cbind(a1, m = rep(1, 12))
c1 <- coords(a1, new.coords = "cartesian")
d1 <- coords(a1, new.coords = "spherical")
e1 <- coords(b1, new.coords = "cartesian")
f1 <- coords(b1, new.coords = "spherical")

a2 <- CMBDataFrame(nside = 1, ordering = "nested", coords = "cartesian")
b2 <- cbind(a2, m = rep(1, 12))
c2 <- coords(a2, new.coords = "cartesian")
d2 <- coords(a2, new.coords = "spherical")
e2 <- coords(b2, new.coords = "cartesian")
f2 <- coords(b2, new.coords = "spherical")

a3 <- CMBDataFrame(nside = 1, ordering = "nested")
b3 <- cbind(a3, m = rep(1, 12))
c3 <- coords(a3, new.coords = "cartesian")
d3 <- coords(a3, new.coords = "spherical")
e3 <- coords(b3, new.coords = "cartesian")
f3 <- coords(b3, new.coords = "spherical")

a4 <- CMBDataFrame(nside = 1, ordering = "nested")[,-1]
b4 <- cbind(a4, m = rep(1, 12))
c4 <- coords(a4, new.coords = "cartesian")
d4 <- coords(a4, new.coords = "spherical")
e4 <- coords(b4, new.coords = "cartesian")
f4 <- coords(b4, new.coords = "spherical")

a <- CMBDataFrame(nside = 1, ordering = "ring")
str(a)
str(a[,1])
str(a[1,])
str(a[1])

b <- CMBDataFrame(nside = 1, ordering = "ring", coords = "cartesian")
pix(b) <- sample(1:12, 12)
str(b)
str(b[,1])
str(b[1,])
str(b[1])
pix(b[1])
pix(b[1,])
pix(b[,1])
pix(b[1:2])
pix(b[4:7,])
pix(b[,2:3])




##################################################################
################ TESTING rbind & cbind ###########################
##################################################################

a <- CMBDataFrame(nside = 1, ordering = "nested")
w1 <- CMBWindow(theta = 0, phi = 0, r = 1)
w2 <- CMBWindow(theta = pi, phi = 0, r = 0.1)
a.w1 <- window(a, new.window = w1)
a.w2 <- window(a, new.window = w2)

a <- CMBDataFrame(nside = 16, ordering = "nested", coords = "cartesian")
w1 <- CMBWindow(theta = 0, phi = 0, r = 0.1)
w2 <- CMBWindow(theta = pi, phi = 0, r = 0.1)
w3 <- CMBWindow(theta = 0, phi = 0, r = 1)
a.w1 <- window(a, new.window = w1)
a.w2 <- window(a, new.window = w2)
a.w3 <- window(a, new.window = w3)

a.wins <- rbind(a.w1, a.w2)
all(pix(a.wins) == c(pix(a.w1), pix(a.w2)))
rbind(a.w1, a.w3) #correctly causes an error

cbind(a.w1, a.w2) #correctly causes an error

cbind(a.w1, m = 1:4)
cbind(m = 1:4, a.w1)




##################################################################
################ DEMONSTRATE rcosmo ##############################
##################################################################

?CMBDataFrame
cmbdf <- CMBDataFrame(nside = 32, ordering = "ring", coords = "spherical")
plot(cmbdf)


ns <- 2 # try: 1, 2, 16
plot(cmbdf, back.col = "black")
plotHPBoundaries(nside = ns, col = "red")
base <- CMBDataFrame(nside = ns, ordering = "ring", coords = "spherical")
plot(base, add = TRUE, size = 5, col = "yellow")




sky <- CMBDataFrame("../CMB_map_smica1024.fits", coords = "spherical")

sky

attributes(sky)



?CMBWindow
win.p <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                   theta = c(pi/2, pi/2, pi/10, pi/10))
plot(win.p, size = 2, back.col = "black")

attributes(win.p)
area(win.p)
maxDist(win.p)




win.d <- CMBWindow(theta = 0, phi = 0, r = 0.1)
plot(win.p, size = 2, back.col = "black")
plot(win.d, size = 2)




sky.excl <- window(sky, list(win.p, win.d))

area(sky.excl)
area(win.d) + area(win.p)

plot(sky.excl, add = TRUE, back.col = "black")





## SPIRAL
dt <- pi/5; dp <- pi/16; end <- pi + pi/20
theta <- rev(c(2*dt, 4*dt, 4*dt, 2*dt, dt, dt, 2*dt, 3*dt,
               3*dt, dt, 0, dt, 2*dt, end))
phi <- rev(c(0, 6*dp, 12*dp, 14*dp, 12*dp, 7*dp, 6*dp, 9*dp,
             5*dp, 4*dp, 8*dp, 15*dp, end, 14*dp))

spiral <- CMBWindow(theta = theta, phi = phi)
sky.spiral <- window(cmbdf, new.window = spiral)

plot(sky.spiral, back.col = "black")
plot(spiral, col = "red", add = TRUE, size = 1.5)
a <- lapply(triangulate(spiral), plot, col = "yellow")






?covCMB

sky.disc <- window(sky, new.window = CMBWindow(theta = 0, phi = 0, r = 0.05))
plot(sky.disc, back.col = "black")

cv <- covCMB(sky.disc, num.bins = 10)

plot(cv[,1], cv[,2], type = "b", xlab = "Distance",
     ylab = "Covariance", main = "Empirical Covariance of Disk")

cv













#####################################################################
######## DEMONSTRATE nest2ring and ring2nest ########################
#####################################################################

## This is a good test too
pix <- 1:48
all.equal((ring2nest(2, pix))[nest2ring(2,pix)],
          (nest2ring(2, pix))[ring2nest(2,pix)],
          pix)

pix <- sample(1:48, 5)
all.equal((ring2nest(2, pix))[nest2ring(2,pix)],
          (nest2ring(2, pix))[ring2nest(2,pix)],
          pix)

ring2nest(2, pix)
nest2ring(2, pix)



# Demonstrate ring order

cmbdf <- CMBDataFrame(nside = 8, ordering = "ring")
plot(cmbdf, type = 'l', col = "black", back.col = "white")
tolabel <- c(1,100:107,768)
plot(cmbdf[tolabel,], labels = tolabel, col = "red", add = TRUE)

#####################################################################
######### DEMONSTRATE plot HealpixBoundaries ########################
#####################################################################

# With dots
ns <- 1
cmbdf <- CMBDataFrame(nside = ns, ordering = "nested",
                      coords = "spherical")
plot(cmbdf, back.col = "black", col = "blue", size = 6)
plotHPBoundaries(ns, col = "red")
plotHPBoundaries(2, col = "green")


# With pixel indices NESTED
plotHPBoundaries(2, col = "red", ordering = "nested")

# With pixel indices RING
plotHPBoundaries(2, col = "red", ordering = "ring")


## Stage 1 (nside 1)
# Plot individual pixels in different colours and only label a few
cmbdf2 <- CMBDataFrame(nside = 64, ordering = "nested")
plot(window(cmbdf2,in.pixels = c(1,9)), col = "blue", back.col = "white")
plot(window(cmbdf2,in.pixels = c(2,10)), col = "green", add = TRUE)
#plot(window(cmbdf2,in.pixels = 3), col = "yellow", add = TRUE)
plot(window(cmbdf2,in.pixels = c(4,12)), col = "orange", add = TRUE)
#plot(window(cmbdf2,in.pixels = 5), col = "orange", add = TRUE)
#plot(window(cmbdf2,in.pixels = 6), col = "red", add = TRUE)
#plot(window(cmbdf2,in.pixels = 7), col = "white", add = TRUE)
#plot(window(cmbdf2,in.pixels = 8), col = "gray", add = TRUE)
#plot(window(cmbdf2,in.pixels = 9), col = "blue", add = TRUE)
#plot(window(cmbdf2,in.pixels = 10), col = "green", add = TRUE)
#plot(window(cmbdf2,in.pixels = 11), col = "yellow", add = TRUE)
#plot(window(cmbdf2,in.pixels = 12), col = "orange", add = TRUE)
plotHPBoundaries(nside = 1, ordering = "nested",
                 incl.labels = c(1,2,3,4,5,6,9), col = "red")




## Stage 2 (nside 2)
# Scale the locations so the pixels don't overshadow the lines
cmbdf3 <- CMBDataFrame(nside = 256, ordering = "nested",
                       coords = "cartesian")
cmbdf3.w <- window(cmbdf3,in.pixels = 1)
cmbdf3.w[,c("x","y","z")] <- cmbdf3.w[,c("x","y","z")]*0.99

# Colour base pixel 1 strongly
plot(cmbdf3.w, col = "light blue",
     back.col = "white", add = TRUE, size = 1.2)
plotHPBoundaries(nside = 2, col = "black", lwd = 1,
                 ordering = "nested",
                 nums.col = "red", incl.labels = c(1,2,3,4))

# Add center labels for the pixels at nside 2 in base pixel 1
# centers <- rcosmo::CMBDataFrame(nside = 2, ordering = "nested",
#                                 coords = "cartesian")
# rcosmo:::plot.CMBDataFrame(centers[c(1,2,3,4),], add = TRUE, col = "red",
#                            labels = c(1,2,3,4), font = 2)

# Plot boundaries and add pixel colours
plotHPBoundaries(nside = 2, col = "red")
plot(window(cmbdf2,in.pixels = 2), col = "green", add = TRUE)
plot(window(cmbdf2,in.pixels = 4), col = "purple", add = TRUE)
plot(window(cmbdf2,in.pixels = 5), col = "orange", add = TRUE)
plot(window(cmbdf2,in.pixels = 6), col = "red", add = TRUE)



##Stage 1 (nside 1) version 2
plot(window(cmbdf2,in.pixels = 1), col = "blue", back.col = "white")
plot(window(cmbdf2,in.pixels = 2), col = "green", add = TRUE)
#plot(window(cmbdf2,in.pixels = 3), col = "yellow", add = TRUE)
plot(window(cmbdf2,in.pixels = 4), col = "orange", add = TRUE)
#plot(window(cmbdf2,in.pixels = 5), col = "orange", add = TRUE)
#plot(window(cmbdf2,in.pixels = 6), col = "red", add = TRUE)
#plot(window(cmbdf2,in.pixels = 7), col = "white", add = TRUE)
#plot(window(cmbdf2,in.pixels = 8), col = "gray", add = TRUE)
plot(window(cmbdf2,in.pixels = 9), col = "blue", add = TRUE)
plot(window(cmbdf2,in.pixels = 10), col = "green", add = TRUE)
#plot(window(cmbdf2,in.pixels = 11), col = "yellow", add = TRUE)
plot(window(cmbdf2,in.pixels = 12), col = "orange", add = TRUE)
plotHPBoundaries(nside = 1, ordering = "nested",
                 incl.labels = c(1,2,4,5,6,9,10,12), col = "black",
                 nums.col = "red")

plotHPBoundaries(nside = 1, ordering = "nested",
                 incl.labels = c(3,7,8,11), col = "black",
                 nums.size = 0.7, nums.col = "red")

# This function labels pixels according to the
# way that nestSearch sees them (by distance)
label <- function(x, bpix)
{
  dists <-  apply(bpix, MARGIN = 1,
                  function(bp) {
                    sum((bp - x)^2)
                  } )
  return(which.min(dists))
}

# Generate base pixels
base <- CMBDataFrame(nside = 1,
                     coords = "cartesian")[,c("x","y","z")]

# Generate pixels at high resolution within base pixels 1, 4, 5
cmbdf <- window(CMBDataFrame(nside = 128, coords = "cartesian"),
                in.pixels = c(1,4,5), in.pixels.res = 0)

# Apply the label function to colour pixels the way that nestSearch sees them
cols <- apply(cmbdf[,c("x","y","z")], MARGIN = 1,
              label, bpix = base)
cols[cols == 1] <- 8 # remove black colour
cmbdf[,c("x","y","z")] <- cmbdf[,c("x","y","z")]*0.99

# Plot all pixels at target resolution
plot(cmbdf, col = cols, size = 5, back.col = "white")

# Plot base pixel boundaries
plotHPBoundaries(nside = 1, lwd = 7, col = "black",
                 ordering = "nested", nums.size = 1.7)

# The problematic pixel is located at (theta, phi) = (0.95, 0)
theta <- 0.95

# Conduct nested search for the nearest neighbour to
# the problematic pixel, starting at base pixel resolution
nestSearch(c(sin(theta), 0, cos(theta)), 32, j = 0:log2(32),
           demo.plot = TRUE)

#####################################################################
######### TEST window() for cmbdf ###################################
#####################################################################


# window on cmbdf
a <- CMBDataFrame(nside = 1, ordering = "nested")
w11 <- CMBWindow(x = 0, y = 0, z = 1, r = 1)
w12 <- CMBWindow(theta = 0, phi = 0, r = 1)
w13 <- CMBWindow(theta= pi/2, phi = 0, r = 1)
w21 <- CMBWindow(x = 0, y = 0, z = -1, r = 1)
w22 <- CMBWindow(theta = pi, phi = 0, r = 1)

a.w11 <- window(a, new.window = w11)
a.w12 <- window(a, new.window = w12)
a.w13 <- window(a, new.window = w13)
a.w21 <- window(a, new.window = w21)
a.w22 <- window(a, new.window = w22)

all(pix(a.w11) == pix(a.w12) & pix(a.w12) == c(1,2,3,4))
all(pix(a.w21) == pix(a.w22) & pix(a.w22) == c(9, 10, 11, 12))

b.w11w12 <- window(a, new.window = list(w11, w12))
all.equal(a.w11, b.w11w12, check.attributes = FALSE)

b.w11w21 <- window(a, new.window = list(w11, w21))
c.w11w21 <- rbind(a.w11,a.w21)
all.equal(c.w11w21, b.w11w21, check.attributes = FALSE)
all.equal(window(c.w11w21), window(b.w11w21))

b.w11w21w13 <- window(a, new.window = list(w11, w21, w13))
c.w11w21w13 <- rbind(a.w11, a.w21, a.w13)
c.w11w21w13 <- c.w11w21w13[order(pix(c.w11w21w13)),] # SHOW THIS PARTICULAR TRICK IN DOCUMENTATION
all.equal(c.w11w21w13, b.w11w21w13, check.attributes = FALSE)
all.equal(window(c.w11w21w13), window(b.w11w21w13))

#####################################################################
######### DEMONSTRATE window ########################################
#####################################################################


## Use the in.pixels argument of window.
## We start by using plotHPBoundaries to determine which pixel the
## window lies in
sky <- readRDS("C:/Users/danie/Downloads/CMB_map_smica_1024.Rda")
fullsky <- readRDS("C:/Users/danie/Downloads/CMB_map_smica_2048.Rda")

size <- 0.5 # The horizon at recombination is today ~1degree on the sky = ~0.02 radians
win <- CMBWindow(theta = c(pi/2, pi/2, pi/2-size, pi/2-size),
                 phi = c(0,size,size,0))
plotHPBoundaries(nside = 1, ordering = "nested")
plot(win)
pixelWin <- pixelWindow(j1 = 0, j2 = log2(1024), pix.j1 = c(1,5))
plot(sky[pixelWin,])
skywin <- window(skyfull, new.window = win, in.pixels = c(1,5))
plot(skywin)





# It is much quicker if you already have cartesian coords
# because no call to pix2coords_internal
cmb1 <- CMBDataFrame(nside = 16, ordering = "nested")
cmb2 <- CMBDataFrame(nside = 16, ordering = "nested", coords = "cartesian")
w <- CMBWindow(x = 0, y = 0, z = 1, r = 1)
library(microbenchmark)
microbenchmark(
  window(cmb1, new.window = w),
  window(cmb2, new.window = w)
)



cmbdf <- CMBDataFrame(nside = 128, ordering = "nested",
                      coords = "cartesian")
plot(cmbdf, back.col = "black")


## LIST OF DISC WITH POLYGON
win.p <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                  theta = c(pi/2, pi/2, pi/10, pi/10))
plot(win.p)
win.d <- CMBWindow(theta = 0, phi = 0, r = 0.1)
plot(win.d)
excl <- window(cmbdf, list(win.p, win.d))
plot(excl, size = 1.2, col = "red")


## DISK WINDOW
win1 <- CMBWindow(x = 1, y = 0, z = 0, r = 0.5, set.minus = FALSE)
disc <- window(cmbdf, win1)
plot(cmbdf, back.col = "black")
plot(disc, col = "red", size = 1.2, add = TRUE)
plot(win1, size = 1.2, col = "yellow")

## MINUS.DISK WINDOW
win2 <- CMBWindow(x = 1, y = 0, z = 0, r = 1, set.minus = TRUE)
disc <- window(cmbdf, win2)
plot(cmbdf, back.col = "black")
plot(disc, col = "red", size = 1.2, add = TRUE)
plot(win2, size = 1.2, col = "yellow")

## LIST OF DISC WITH MINUS.DISC (UNION)
discs <- window(cmbdf, list(win1, win2), intersect = FALSE) # use win1, win2 above
plot(cmbdf, back.col = "black")
plot(discs, col = "red", size = 1.2, add = TRUE)
plot(win1, size = 1.2, col = "yellow")
plot(win2, size = 1.2, col = "yellow")


######### SLOW!!!!
## SMALL WINDOW ON SPARSE SKY
sky <- CMBDataFrame(nside = 1024, ordering = "nested", coords = "cartesian")
plot(sky, sample.size = 20000, back.col = "black")
win <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                 theta = c(pi/2, pi/2, pi/10, pi/10))
detailed <- window(sky, new.window = win)
plot(detailed, add = TRUE, size = 1)
plot(win)



## list of overlapping polygons
win1 <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                  theta = c(pi/2, pi/2, pi/10, pi/10))
win2 <- CMBWindow(phi = c(0, pi/5, pi/5, 0),
                  theta = c(pi/4, pi/4, pi/5, pi/5))
cmbdf.win <- window(cmbdf, new.window = list(win1, win2))
plot(cmbdf.win)
plot(win1, col = "yellow")
plot(win2, col = "yellow")

### list of overlapping minus.polygons
win1 <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                  theta = c(pi/2, pi/2, pi/10, pi/10),
                  set.minus = TRUE)
win2 <- CMBWindow(phi = c(0, pi/5, pi/5, 0),
                  theta = c(pi/4, pi/4, pi/5, pi/5),
                  set.minus = TRUE)
cmbdf.win <- window(cmbdf, new.window = list(win1, win2))
plot(cmbdf.win, add = TRUE, col = "red", size = 1.2)
plot(win1, col = "yellow")
plot(win2, col = "yellow")

### list of overlapping polygons and minus.polygons (INTERSECT)
win1 <- CMBWindow(phi = c(0, pi/10, pi/10, 0),
                  theta = c(pi/2, pi/2, pi/10, pi/10),
                  set.minus = FALSE)
win2 <- CMBWindow(phi = c(0, pi/5, pi/5, 0),
                  theta = c(pi/4, pi/4, pi/5, pi/5),
                  set.minus = TRUE)
plot(cmbdf)
cmbdf.win <- window(cmbdf, new.window = list(win2, win1))
plot(cmbdf.win, add = TRUE, col = "red", size = 1.2)
plot(win1, col = "yellow")
plot(win2, col = "yellow")


### list of disjoint polygons and minus.polygons (UNION)
win1 <- CMBWindow(phi = c(0, pi/2, pi/2),
                  theta = c(pi/2, pi/2, 0),
                  set.minus = TRUE)
win2 <- CMBWindow(phi = c(2*pi/10, 2*pi/10, 3*pi/10, 3*pi/10),
                  theta = c(2*pi/10, 3*pi/10, 3*pi/10, 2*pi/10),
                  set.minus = FALSE)
plot(cmbdf)
cmbdf.win <- window(cmbdf, new.window = list(win1,win2), FALSE)
plot(cmbdf.win, col = "green", size = 1.4, add = TRUE)
plot(win1, col = "yellow")
plot(win2, col = "yellow")




## STAR
df <- data.frame(theta = c(pi-pi/40, pi/2+pi/10, pi/2,
                            pi/2-pi/10, 0, pi/2-pi/10, pi/2, pi/2+pi/10),
                 phi = c(3*pi/2, 3*pi/2+pi/10, 0, 3*pi/2+pi/10,
                          0, 3*pi/2-pi/10, pi+pi/40, 3*pi/2-pi/10))
star <- CMBWindow(df)
cmbdf.st <- window(cmbdf, star)
plot(cmbdf, back.col = "black")
plot(cmbdf.st, col = "green", size = 1.5, add = TRUE)
a <- lapply(triangulate(star), plot)


## SPIRAL
dt <- pi/5
dp <- pi/16
end <- pi + pi/20
theta <- rev(c(2*dt, 4*dt, 4*dt, 2*dt, dt, dt, 2*dt, 3*dt,
           3*dt, dt, 0, dt, 2*dt, end))
phi <- rev(c(0, 6*dp, 12*dp, 14*dp, 12*dp, 7*dp, 6*dp, 9*dp,
         5*dp, 4*dp, 8*dp, 15*dp, end, 14*dp))
df <- data.frame(theta = theta, phi = phi)
spiral <- CMBWindow(df)
cmbdf.sp <- window(cmbdf, new.window = spiral)
plot(cmbdf, back.col = "black")
plot(cmbdf.sp, col = "red", add = TRUE, size = 1.5)
a <- lapply(triangulate(spiral), plot, col = "yellow")

## Have a look at the locations and order of the vertices
## to check that they are anti-clockwise
df.xyz <- rcosmo::sph2car(df)
rgl::plot3d(df.xyz, add = TRUE, col = "yellow", size = 10)
rgl::plot3d(df.xyz[1,], add = TRUE, col = "yellow", size = 10)
rgl::plot3d(df.xyz[2,], add = TRUE, col = "yellow", size = 10)
rgl::plot3d(df.xyz[7,], add = TRUE, col = "green", size = 15)



### List of disjoint polygons
mouth <- CMBWindow(theta = c(pi/2,pi/3,pi/2-pi/20,pi/3),
                 phi = c(3*pi/2,3*pi/2 + pi/4, 3*pi/2, 3*pi/2 - pi/4))
eye1 <- CMBWindow(theta = c(pi/4+0.15, pi/4+0.15, pi/4-0.1, pi/4-0.1),
                 phi = c(3*pi/2+pi/8-0.1, 3*pi/2+pi/8+0.1, 3*pi/2+pi/8+0.1,
                         3*pi/2+pi/8-0.1))
eye2 <- CMBWindow(theta = c(pi/4+0.15, pi/4+0.15, pi/4-0.1, pi/4-0.1),
                 phi = c(3*pi/2-pi/8-0.1, 3*pi/2-pi/8+0.1, 3*pi/2-pi/8+0.1,
                         3*pi/2-pi/8-0.1))
face <- list(mouth, eye1, eye2)
cmbdf.face <- window(cmbdf, new.window = face)
plot(cmbdf, back.col = "black")
plot(cmbdf.face, col = "red", add = TRUE, size = 1.5)
a <- lapply(face, plot, col = "yellow")

### List of disjoint minus polygons
mouth <- CMBWindow(theta = c(pi/2,pi/3,pi/2-pi/20,pi/3),
                   phi = c(3*pi/2,3*pi/2 + pi/4, 3*pi/2, 3*pi/2 - pi/4),
                   set.minus = TRUE)
eye1 <- CMBWindow(theta = c(pi/4+0.15, pi/4+0.15, pi/4-0.1, pi/4-0.1),
                  phi = c(3*pi/2+pi/8-0.1, 3*pi/2+pi/8+0.1, 3*pi/2+pi/8+0.1,
                          3*pi/2+pi/8-0.1),
                  set.minus = TRUE)
eye2 <- CMBWindow(theta = c(pi/4+0.15, pi/4+0.15, pi/4-0.1, pi/4-0.1),
                  phi = c(3*pi/2-pi/8-0.1, 3*pi/2-pi/8+0.1, 3*pi/2-pi/8+0.1,
                          3*pi/2-pi/8-0.1),
                  set.minus = TRUE)
face <- list(mouth, eye1, eye2)
cmbdf.face <- window(cmbdf, new.window = face)
plot(cmbdf, back.col = "black")
plot(cmbdf.face, col = "red", add = TRUE, size = 1.5)
a <- lapply(face, plot, col = "yellow")



## non-convex polygon using triangulate implicitly
win <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                 theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))
plot(win)
cmbdf.wins <- window(cmbdf, new.window = win)
plot(cmbdf.wins, col = "red", add = TRUE, size = 1)
winType(window(cmbdf.wins))
a <- lapply(window(cmbdf.wins), plot)

## non-convex polygon using triangulate explicitly
win <- CMBWindow(phi = c(0, pi/4, pi/4, pi/5),
                 theta = c(pi/2, pi/2, pi/4, pi/2 - pi/20))
wins <- triangulate(win)
cmbdf.wins <- window(cmbdf, new.window = wins)
plot(cmbdf.wins, add = TRUE, size = 2)
winType(window(cmbdf.wins)) # extra polygon because triangulate
a <- lapply(window(cmbdf.wins), plot)


## gnomonic projection of non-convex polygon onto plane z = 1

      ### haven't done this yet


#### BEACH BALL WITH WINDOW BOUNDARIES

## Up-quad1 (counter-clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,0,pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, size = 1)
plot(win1)

## Up-quad2 (counter-clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,pi/2,pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "yellow")
plot(win1)

## Up-quad3 (counter-clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,pi,3*pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "orange")
plot(win1)

## Up-quad4 (counter-clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,3*pi/2, 2*pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "green")
plot(win1)

## Down-quad1 (counter-clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(0,0,pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "green")
plot(win1)

## Down-quad2 (counter-clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(pi/2,pi/2, pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "orange")
plot(win1)

## Down-quad3 (counter-clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(pi,pi,3*pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "yellow")
plot(win1)

## Down-quad4 (counter-clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi,pi/2), phi = c(3*pi/2,3*pi/2,2*pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, col = "blue")
plot(win1)





##Example1 up-quad1 (clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,pi/2,0))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE, add = TRUE, size = 1.5)
plot(win1, axes = TRUE)

##Example1 up-quad2 (clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,pi, pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

##Example1 up-quad3 (clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,3*pi/2, pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

##Example1 up-quad4 (clockwise):
win1 <- CMBWindow(theta = c(0,pi/2,pi/2), phi = c(0,2*pi, 3*pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

## Example1 down-quad1 (clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi/2,pi), phi = c(0,pi/2,0))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

## Example1 down-quad2 (clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi/2,pi), phi = c(pi/2, pi, pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

## Example1 down-quad3 (clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi/2,pi), phi = c(pi,3*pi/2,pi))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)

## Example1 down-quad4 (clockwise):
win1 <- CMBWindow(theta = c(pi/2,pi/2,pi), phi = c(3*pi/2,2*pi,3*pi/2))
cmbdf.win1 <- window(cmbdf, new.window = win1)
plot(cmbdf.win1, axes = TRUE)
plot(win1, axes = TRUE)





## Example2 (anticlockwise): north cap
win2 <- CMBWindow(theta = rep(pi/3, 10),
                  phi = seq(0,9*2*pi/10, length.out = 10))
cmbdf.win2 <- window(cmbdf, new.window = win2)
plot(cmbdf.win2)
plot(win2)

### cLOCKWISE WINDOW EXAMPLE SHOWS THAT IT'S A GOOD IDEA TO PLOT WIN TO CHECK
## Example2 (clockwise): Note the use of 'rev' function here:
win2 <- CMBWindow(theta = rep(pi/3, 10),
                  phi = rev(seq(0,9*2*pi/10, length.out = 10)))
cmbdf.win2 <- window(cmbdf, new.window = win2)
plot(cmbdf.win2)
plot(win2)



## Example3: another non-convex window
win3 <- CMBWindow( theta = c(pi/2, pi/2, pi/4, pi/4, 0 ),
                   phi   = c(   0,    1,    1,  0.5, 0 ) )
cmbdf.win3 <- window(cmbdf, new.window = win3)
plot(cmbdf.win3)
plot(win3)



#####################################################################
######### DEMONSTRATE CMBWindow attributes and helpers ##############
#####################################################################

# disc
win <- CMBWindow(x = 0, y = 0, z = 1, r = pi/2)
winType(win)
area(win)
all.equal(area(win), 4*pi/2)

# minus.disc
win <- CMBWindow(x = 0, y = 0, z = 1, r = 0.5, set.minus = TRUE)
winType(win)

# minus.polygon
win <- CMBWindow(theta = c(0,pi/2,pi/2),
                 phi   = c(0,0   ,pi/2),
                 set.minus = TRUE)
winType(win)


# anticlockwise
win <- CMBWindow(theta = c(0,pi/2,pi/2),
                 phi   = c(0,0   ,pi/2))
winType(win)
area(win)

# clockwise
win2 <- CMBWindow(theta = c(0,pi/2, pi/2),
                  phi   = c(0,pi/2, 0  ))
area(win2)

# clockwise + anticlockwise = area of unit sphere
area(win) + area(win2) == 4*pi

# Note that maxDist does not restrict to travel within window
maxDist(win)
maxDist(win2)


coords(win)
coords(win) <- "cartesian"
win

win <- CMBWindow(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
coords(win) <- "spherical"
win


# Disc window
dwin <- CMBWindow(theta = 0, phi = 0, r = 1)
winType(dwin)
dwin


#####################################################################
######### DEMONSTRATE car2sph and sph2car ###########################
#####################################################################

df <- data.frame(theta = c(1,2,3), phi = c(0,0,0))
df
df.xyz <- sph2car(df)
df <- car2sph(df.xyz)
df
df.xyz <- data.frame(x = c(1,0,0), y = c(0,1,0), z = c(0,0,1))
df <- car2sph(df.xyz)



#####################################################################
######### DEMONSTRATE covCMB ########################################
#####################################################################

### Using simulated data
## Set parameters
ns <- 16
num.bins <- 10
npix <- 12*ns^2
sim.I <- rnorm(npix)

cmbdf.fake <- CMBDataFrame(nside = ns, coords = "cartesian",
                      ordering = "nested",
                      intensities = sim.I)
# Full sky fake data
C <- covCMB(cmbdf.fake)



### Using real data ###########################################################
#sky <- CMBDataFrame("../CMB_map_smica1024.fits", coords = "cartesian")
#plot(sky, sample.size = 500000)
sky <- readRDS("../fullsky.Rds")

## Window
size <- 0.2
win <- CMBWindow(theta = c(pi/2, pi/2, pi/2-size, pi/2-size), phi = c(0,size,size,0))
cmbdf.win <- window(sky, new.window = win)
# plot(cmbdf.win, add = TRUE)
# plot(win)
C <- covCMB(cmbdf.win, max.dist = size, num.bins = 20)
plot(C$dist, C$cov/C$cov[1])

## Full sky sample
sky.s <- sampleCMB(sky, sample.size = 50000)
C.sky <- covCMB(sky.s, max.dist = pi, num.bins = 200)
plot(C.sky$dist, C.sky$cov/C$cov[1])



######### BENCHMARK covCMB_internal1 vs covCMB_internal2 ############

library(microbenchmark)

breaks1 <- seq(0, pi, length.out = 10+1)[-1]
breaks2 <- seq(cos(pi), 1, length.out = 10+1)[-(10+1)]
sky.fake <- CMBDataFrame(nside = 16, coords = "cartesian", ordering = "nested", intensities = rnorm(12*16^2))


mb <- microbenchmark(covCMB_internal1(sky.fake, breaks1),
                     covCMB_internal2(sky.fake, breaks2))
summary(mb)$mean[1]/summary(mb)$mean[2] # Internal2 is about 5 times faster


######### DEMONSTRATION OF SPEED INCREASE FROM Rcpp #################

# R version (using list)
f <- function(cmbdf)
{
  dists <- list()
  n <- nrow(cmbdf)

  for ( i in 1:(n-1) ) {
    #cat("i: ", i, "\n")

    dists[[i]] <- vector(mode = "numeric", length = n-i)

    for ( j in (i+1):n ) {

      # Distances d are stored as d(xi, xj) = dists[[i]][j-i] (where i < j)
      dists[[i]][j-i] <- geoDist(cmbdf[i, c("x","y","z")],
                                 cmbdf[j, c("x","y","z")])

    }
  }

  return(dists)
}

# Test
library(microbenchmark)
cmbdf <- CMBDataFrame(nside = 1, coords = "cartesian",
                      ordering = "nested")
mb <- microbenchmark(f(cmbdf),           # R version
                     geoDistList(cmbdf)) # C++ rcosmo version
# C++ function is about 2500 times faster
summary(mb)$mean[1]/summary(mb)$mean[2]









########### Calculation of geoDistList output size in GB ##############

## Specify n:
n <- 50000
## Else specify nside:
# nside <- 256
# n <- 12*nside^2
## Most of the space is consumed by the sum_{k=1}^{n-1} k doubles
## which take up 8 bytes of space each. When integers are used as
## with distBinList these are 4 bytes so space is halved.
num.doubles <- n*(n-1)/2
rough.size <- (num.doubles*8 + (n-1)*8)/1000000000
rough.size # n = 50,000 gives ~10GB list!















###############################################################################################
######################### CMBDataFrame DEMONSTRATION ##########################################
###############################################################################################

library(rcosmo)
library(magrittr)
library(sphereplot)
library(rcosmo)

# Example 1
df <- CMBDataFrame(nside = 16, ordering = "ring")
class(df)
plot(df)

# Example 2
df <- CMBDataFrame(nside = 16, ordering = "nested",
                   intensities = rnorm(12*16^2),
                   coords = "cartesian")
plot(df)

# Example 3
sky <- CMBDataFrame(CMBData = "../CMB_map_smica1024.fits")
sky.sample <- sampleCMB(sky, sample.size = 1000000)
plot(sky.sample)

# Example 4 (same effect as example 3)
sky.sample2 <- CMBDataFrame(df, sample.size = 1000000, ordering = "ring",
                            coords = "cartesian")
plot(sky.sample2)


sky
summary(df)
summary(sky)

ordering(df)
ordering(df) <- "ring"
ordering(df, newOrdering = "ring")

coords(df)
coords(df) <- "cartesian"
coords(df, newCoords = "cartesian")

head(pix(df))
#pix(df) <- some_vector
nside(df)


ls("package:rcosmo")
lsf.str("package:rcosmo")








########################################################################
############### OLD ATTEMPTED POINT IN POLYGON DEMO ####################
########################################################################


# BOUNDARY POINT FINDER FUNCTION
# Works by rotating the two vertices of each line to the equator before sampling, then rotating back
polygonBoundary <- function( vertices.xyz, eps = 0.001 )
{
  rows <- nrow(vertices.xyz)
  if(rows < 3){
    stop("Polygon must have at least 3 vertices.");
  }
  #vertices_xyz_cycled <- rbind(vertices_xyz[2:nrow(vertices_xyz),], vertices_xyz[1,])

  boundary <- data.frame()
  for ( row in 1:rows )
  {
    #Geodesic is the shortest distance between two_points
    V1 <- vertices.xyz[row,]
    V2 <- vertices.xyz[1 + (row %% rows),]

    normal <- vector_cross(V1, V2)
    rotated <- rodrigues(as.matrix(normal)[1,], c(0,0,1), rbind(V1,V2))

    two.longitudes <- atan2(rotated[,2], rotated[,1]) # latitudes are both pi/2

    line.longitudes <- seq(two.longitudes[1], two.longitudes[2], by = eps)
    line <- data.frame(phi = line.longitudes, theta = rep(pi/2,length(line.longitudes)))
    line <- sph2car( line )

    line.rotated <-  as.data.frame(rodrigues(c(0,0,1), as.matrix(normal)[1,], line))
    names(line.rotated) <- c("x","y","z")
    boundary <- rbind(boundary, line.rotated)
  }
  boundary
}

## TEST BOUNDARY POINT FINDER FUNCTION
cmbdf <- CMBDataFrame(nside = 16, ordering = "ring",
                      coords = "cartesian")
plot(cmbdf)

vertices.sph <- data.frame( theta = c(pi/2, pi/2, pi/4, pi/4, 0 ),
                            phi   = c(   0,    1,    1,  0.5, 0 ) )
vertices.xyz <- sph2car(vertices.sph)
rgl::plot3d(vertices.xyz, col = 'red', type = 'p', size = 12, pch = 2, add = TRUE)

boundary.xyz <- polygonBoundary( vertices.xyz, 0.001 )
rgl::plot3d( boundary.xyz, col = 'red', type = 'p', size = 3.2, pch = 3, add = TRUE)
boundary <- car2sph(boundary.xyz)



######## THE REST OF THIS CODE HASN'T BEEN UPDATED AND PROBABLY DOESN'T WORK

# Get ALL HEALPix points in the square
df_square <- df[  df$theta  <= max(vertices_sph$theta)  & df$theta  >= min(vertices_sph$theta)
                  & df$phi <= max(vertices_sph$phi) & df$phi >= min(vertices_sph$phi),]
df_square_xyz <- as.data.frame(sph2car(df_square$phi, df_square$theta, deg = FALSE))
plot3d(df_square_xyz, col = 'blue', type = 'p', size = 1.6, pch = 3, add = TRUE)


head(polygon_boundary_points)
df_square


### THIS SECTION ONLY USEFUL IF THE SAMPLE SIZE IS SMALL I.E. NOT DENSE IN SKY ###
# Round latitudes to some tolerance based on eps then split into lists by latitude
eps <- 0.001
digits <- ceiling(log10(1/eps) + 1)
polygon_boundary_points$theta <- round(polygon_boundary_points$theta, digits)
polygon_boundary_points_split <- split(polygon_boundary_points, polygon_boundary_points$theta)
polygon_boundary_lats <- as.numeric(names(polygon_boundary_points_split))
# Get isolatitude values from the square
isolats_in_square <- unique(df_square$theta)
# For each latitude in the polygon, find the nearest isolat in the square
for ( polylat in polygon_boundary_lats )
{
  polygon_boundary_points_split[[as.character(polylat)]]$theta <-
    isolats_in_square[ which.min(abs(isolats_in_square - polylat)) ]
}
# Bind the list pack together
library(data.table)
polygon_boundary_points <- rbindlist(polygon_boundary_points_split) #faster than do.call("rbind",...)
# Plot result to check sanity
polygon_boundary_points_xyz <- as.data.frame(sph2car(
  polygon_boundary_points$phi,
  polygon_boundary_points$theta,
  deg = FALSE))
plot3d(polygon_boundary_points_xyz, col = 'yellow', type = 'p', size = 5, pch = 3, add = TRUE)
#-----------------------------------------------------------------------------------





###############################################################################################
####################################### OLD DEMONSTRATION CODE ################################
###############################################################################################


rm(list=ls())
library(rcosmo)
library(Rcpp)
library(rgl)
library(sphereplot) # for sph2car()
#help(rgl)
library(FITSio)
library(R.matlab) # for importing the colour map
library(geosphere)
sourceCpp("pix2ang.cpp")
sourceCpp("distGeo.cpp")
sourceCpp("nestIndex.cpp")
source("readFITScmb.R")
source("CMBDataFrame.R")
# Function takes longitude in [0,2pi] and transforms it
# to longitude in [-pi,pi]
lonWrap <- function(lon) {
  (lon + 180) %% 360 - 180
}












##### PLAN (TO DO) #####
# Select only the points on the square that lie on the planes that
# define the geodesics, since they are already on the sphere they
# must lie on the geodesics. Also, they already must lie on the
# correct isolatitude rings. Then figure out how to cut off the excess.
vertices_sph <- data.frame( phi = c(0, 1, 1, 0.5, 1, 0), theta = c(0, 0, 0.2, 0.2, 0.5, 1) )
vertices_xyz <- as.data.frame(sph2car(vertices_sph$phi, vertices_sph$theta, deg = FALSE))
df_square <- df[  df$theta  <= max(vertices_sph$theta)  & df$theta  >= min(vertices_sph$theta)
                  & df$phi <= max(vertices_sph$phi) & df$phi >= min(vertices_sph$phi),]
df_square_xyz <- as.data.frame(sph2car(df_square$phi, df_square$theta, deg = FALSE))
plot3d(df_square_xyz, col = 'blue', type = 'p', size = 3.3, pch = 3, add = TRUE)

vertices <- nrow(vertices_xyz)
boundary_points <- data.frame()
for ( point in 1:nrow(df_square_xyz) )
{
  p <- df_square_xyz[point,]
  for ( vertex in 1:vertices )
  {

    # two_points <- rbind(vertices_xyz[vertex,],
    #                     vertices_xyz[1 + (vertex %% vertices),])

    # if (u x v) dot p = 0
    # if ( isTRUE(all.equal(sum(vector_cross(two_points[1,], two_points[2,])*p), 0)) )
    # {
    #   boundary_points <- rbind(boundary_points, p)
    #   break
    # }
  }
}
names(boundary_points) <- c("x","y","z")
plot3d(boundary_points, col = 'green', type = 'p', size = 7, pch = 3, add = TRUE)



# POINT IN POLYGON --------------------------------------------------------

# Import data
df <- CMBDataFrame("../CMB_map_smica1024.fits", sampleSize = 100000)
plotCMB(df) # Currently if the file is in nested and the user specifies ordering = ring
            # then CMBDataFrame doesnt automatically
            # convert to ring, so we should add that (or at least throw a warning)

# TEST LATEST -------------------------------------------------------------
# Trying to find polygon edge pixels
source("exploration/rodrigues.R") # Rodrigues and vector_cross

# BOUNDARY POINT FINDER FUNCTION
# Works by rotating the two vertices of each line to the equator before sampling, then rotating back
polygon_boundary <- function( vertices_xyz, eps = 0.001 )
{
  rows <- nrow(vertices_xyz)
  if(rows < 3){
    stop("Polygon must have at least 3 vertices.");
  }
  #vertices_xyz_cycled <- rbind(vertices_xyz[2:nrow(vertices_xyz),], vertices_xyz[1,])

  boundary_points <- data.frame()
  for ( row in 1:rows )
  {
    #Geodesic is the shortest distance between two_points
    two_points <- rbind(vertices_xyz[row,], vertices_xyz[1 + (row %% rows),])
    normal_vector <- vector_cross(two_points[1,], two_points[2,])
    two_points_rotated <- rodrigues(as.matrix(normal_vector)[1,], c(0,0,1), two_points)
    two_longitudes <- atan2(two_points_rotated[,2], two_points_rotated[,1]) # latitudes are both 0

    line_longitudes <- seq(two_longitudes[1], two_longitudes[2], by = eps)
    line <- as.data.frame( sph2car( phi = line_longitudes,
                                    theta = rep(0,length(line_longitudes)),
                                    deg = FALSE ) )

    line_rotated_back <-  as.data.frame(rodrigues(c(0,0,1), as.matrix(normal_vector)[1,], line))
    names(line_rotated_back) <- c("x","y","z")
    boundary_points <- rbind(boundary_points, line_rotated_back)
  }
  boundary_points
}

# Use boundary point finder function to find some boundary points for a polygon
vertices_sph <- data.frame( phi = c(0, 1, 1, 0.5, 1, 0), theta = c(0, 0, 0.2, 0.2, 0.5, 1) )
vertices_xyz <- as.data.frame(sph2car(vertices_sph$phi, vertices_sph$theta, deg = FALSE))
plot3d(vertices_xyz, col = 'red', type = 'p', size = 15, pch = 3, add = TRUE)
polygon_boundary_points_xyz <- polygon_boundary( vertices_xyz, 0.001 )
plot3d(polygon_boundary_points_xyz, col = 'red', type = 'p', size = 3.2, pch = 3, add = TRUE)
polygon_boundary_points <- as.data.frame( car2sph(polygon_boundary_points_xyz$x,
                                                  polygon_boundary_points_xyz$y,
                                                  polygon_boundary_points_xyz$z,
                                                  deg = FALSE)[,c("phi","theta")] )




# Get a plane through two of the points for comparison geodesic
two_points <- sph2car(c(1,0), c(0.5,1), deg = FALSE)
normal_vector <- vector_cross(two_points[1,], two_points[2,])
planes3d(normal_vector[1], normal_vector[2], normal_vector[3]) # Geodesic should follow this plane


# Select square with lowest to highest both theta and long, plot it
df_square <- df[  df$theta  <= max(vertices_sph$theta)  & df$theta  >= min(vertices_sph$theta)
                & df$phi <= max(vertices_sph$phi) & df$phi >= min(vertices_sph$phi),]
df_square_xyz <- as.data.frame(sph2car(df_square$phi, df_square$theta, deg = FALSE))
row.names(df_square_xyz) <- row.names(df_square) # This is important! sph2car ditches the row names
plot3d(df_square_xyz, col = 'blue', type = 'p', size = 3.3, pch = 3, add = TRUE)


head(polygon_boundary_points)
head(df_square)


#############################################################################################
### THIS SECTION ONLY USEFUL IF THE SAMPLE SIZE IS SMALL I.E. NOT DENSE IN SKY ###
# Round latitudes to some tolerance based on eps then split into lists by latitude
eps <- 0.001
digits <- ceiling(log10(1/eps) + 1)
polygon_boundary_points$theta <- round(polygon_boundary_points$theta, digits)
polygon_boundary_points_split <- split(polygon_boundary_points, polygon_boundary_points$theta)
polygon_boundary_lats <- as.numeric(names(polygon_boundary_points_split))
# Get isolatitude values from the square
isolats_in_square <- unique(df_square$theta)
# For each latitude in the polygon, find the nearest isolat in the square
for ( polylat in polygon_boundary_lats )
{
  polygon_boundary_points_split[[as.character(polylat)]]$theta <-
    isolats_in_square[ which.min(abs(isolats_in_square - polylat)) ]
}
# Bind the list pack together
library(data.table)
polygon_boundary_points <- rbindlist(polygon_boundary_points_split) #faster than do.call("rbind",...)
# Plot result to check sanity
polygon_boundary_points_xyz <- as.data.frame(sph2car(
  polygon_boundary_points$phi,
  polygon_boundary_points$theta,
  deg = FALSE))
plot3d(polygon_boundary_points_xyz, col = 'yellow', type = 'p', size = 5, pch = 3, add = TRUE)
################################################################################################

#loop through all points in the square using a raycasting algorithm
rows <- nrow(df_square)
for ( row in 1:rows )
{
  inside <- "unknown"
  # If we are at the start of an isolat ring
  if ( row > 1
     & df_square$theta[row] != df_square$theta[row-1] )
  {
    # Check if we start inside a polygon

  }
}


rows <- nrow(df_square)
thislat <- df_square$theta[1]
inside <- "unknown"
count <- 0
for ( row in 1:rows )
{
  prevlat <- thislat
  thislat <- df_square$theta[row]

  # We are at the start of an isolat ring if...
  if ( thislat != prevlat )
  {
    inside <- "unknown"
  }

  if ( thislat  )

}






# OLD / WORK IN PROGRESS --------------------------------------------------------


# Parametric form for great circle joining u,v is ucos(t) + wsin(t)
# where w = (u x v) x u, since then u, w are orthonormal vectors in the plane of the circle
vertices_sph <- data.frame( phi = c(0, 1, 1, 0.5, 1, 0), theta = c(0, 0, 0.2, 0.2, 0.5, 1) )
vertices_xyz <- as.data.frame(sph2car(vertices_sph$phi, vertices_sph$theta, deg = FALSE))
u <- vertices_xyz[1,]
v <- vertices_xyz[2,]
w <- vector_cross(vector_cross(u,v),u)
w <- w/sqrt(sum(w^2))
t <- seq(0,pi/2,0.01)
u <- matrix(as.numeric(u), byrow = TRUE, ncol = 3, nrow = length(t))
w <- matrix(as.numeric(w), byrow = TRUE, ncol = 3, nrow = length(t))
great_circ_uv <- u*cos(t) + w*sin(t)
plot3d(great_circ_uv, col = 'blue', type = 'p', size = 3.2, pch = 3, add = TRUE)


# Specify the polygon in spherical, convert to cartesian, plot it
polygon_sph <- data.frame( phi = c(0.2, 0.4, 0.4, 0.2, 0.2), theta = c(0.2,0.2,0.4,0.4, 0.2) )
polygon_xyz <- as.data.frame(sph2car(polygon_sph$phi, polygon_sph$theta, deg = FALSE))
plot3d(polygon_xyz, col = 'white', type = 'l', size = 20, cex = 30, pch = 3, add = TRUE)

# These are the high/low-est lat/phi of the polygon vertices
highest_lat <- max(polygon_sph$theta)
lowest_lat <- min(polygon_sph$theta)
highest_long <- max(polygon_sph$phi)
lowest_long <- min(polygon_sph$phi)

# Select the strip from lowest to highest lat, plot it
df_strip <- df[df$theta < highest_lat & df$theta > lowest_lat, ]
df_strip_xyz <- as.data.frame(sph2car(df_strip$phi, df_strip$theta, deg = FALSE))
plot3d(df_strip_xyz, col = 'green', type = 'p', size = 7, cex = 30, pch = 3, add = TRUE)

# Select strip with lowest to highest both lat and phi, plot it
df_square <- df[df$theta <= highest_lat & df$theta >= lowest_lat & df$phi <= highest_long & df$phi >= lowest_long,]
df_square_xyz <- as.data.frame(sph2car(df_square$phi, df_square$theta, deg = FALSE))
plot3d(df_square_xyz, col = 'blue', type = 'p', size = 8, pch = 3, add = TRUE)

# Stereographic projection of the spherical polygon vecrtices and contents
df_plane <- data.frame(x = df_square_xyz$x/(1 - df_square_xyz$z),
                       y = df_square_xyz$y/(1 - df_square_xyz$z),
                       z = rep(1, dim(df_square_xyz)[1]) )
polygon_plane <- data.frame(x = polygon_xyz$x/(1-polygon_xyz$z),
                            y = polygon_xyz$y/(1-polygon_xyz$z),
                            z = rep(1, dim(polygon_xyz)[1]) )
plot3d(df_plane, col = 'blue', type = 'p', size = 8, pch = 3)
plot3d(polygon_plane, col = 'yellow', type = 'l', size = 8, pch = 3, add = TRUE)

# Trying to find polygon edge pixels
library(Rcpp)
sourceCpp("exploration/distGeo.cpp") # For segment distances to determine tolerance

# Just draw one line for practice
distGeo(as.matrix(data.frame(x = 1, y = 1)), as.matrix(data.frame(x = 1, y = 1)))
eps <- 0.001
dist <- distGeo(as.matrix(data.frame(phi = 0.2, theta = 0.2)),
                as.matrix(data.frame(phi = 0.4, theta = 0.4)) )
n <- dist/eps

line <- data.frame(phi = c(seq(0.2,0.4, length.out = n)),
                   theta = c(seq(0.2,0.4, length.out = n)) )
line_xyz <- as.data.frame(sph2car(line$phi, line$theta, deg = FALSE))
plot3d(line_xyz, col = 'green', type = 'p', size = 10, pch = 3, add = TRUE)

## Sample simple polygon boundary
polygon_sph <- polygon_sph[1:(nrow(polygon_sph)-1),] # We no longer need the extra (0.2,0.2) row
# BOUNDARY POINT FINDER FUNCTION
polygon_boundary <- function( vertices_sph, eps = 0.001 )
{
  vertices_sph_rot <- rbind(vertices_sph[2:nrow(vertices_sph),], vertices_sph[1,])
  dist <- distGeo(as.matrix(vertices_sph), as.matrix(vertices_sph_rot))

  boundary_points <- data.frame()
  for ( row in 1:nrow(vertices_sph) )
  {
    n <- dist[row]/eps
    line <- data.frame(phi = seq(vertices_sph$phi[row], vertices_sph_rot$phi[row], length.out = n),
                       theta = seq(vertices_sph$theta[row], vertices_sph_rot$theta[row], length.out = n) )
    boundary_points <- rbind(boundary_points, line)
  }

  boundary_points
}
# Plot now
polygon_points <- polygon_boundary( polygon_sph )
polygon_points_xyz <- as.data.frame(sph2car(polygon_points$phi, polygon_points$theta, deg = FALSE))
plot3d(polygon_points_xyz, col = 'red', type = 'p', size = 10, pch = 3, add = TRUE)

#Try it with a more complicated polygon
polygon_sph2 <- data.frame( phi = c(0, 1, 1, 0.5, 1, 0), theta = c(0, 0, 0.2, 0.2, 0.5, 1) )
polygon_points2 <- polygon_boundary( polygon_sph2 )
polygon_points_xyz2 <- as.data.frame(sph2car(polygon_points2$phi, polygon_points2$theta, deg = FALSE))
plot3d(polygon_points_xyz2, col = 'red', type = 'p', size = 5, pch = 3, add = TRUE)

## Check that my lines are geodesics
# FUNCTION FOR CROSS PRODUCT AND RODRIGUES FORMULA
source("exploration/rodrigues.R")
# Test: Check the Rodrigues function translates polygon to north pole
poly_rot_check <- rodrigues(a = sph2car(0,1, deg = FALSE),
                            b = c(0,0,1),
                            polygon_points_xyz2)
plot3d(poly_rot_check, col = 'red', type = 'p', size = 7, pch = 3, add = TRUE)

# Get points grid on geodesic between two points (phi, theta) = (1,0.5) and (0,1)
two_points <- sph2car(c(1,0), c(0.5,1), deg = FALSE)
plot3d(two_points, col = 'yellow', type = 'p', size = 20, pch = 3, add = TRUE)
normal_vector <- vector_cross(two_points[1,], two_points[2,])
planes3d(normal_vector[1], normal_vector[2], normal_vector[3]) # Geodesic should follow this plane
two_points_rot <- rodrigues(normal_vector, c(0,0,1), two_points)
two_points_rot_long <- atan2(two_points_rot[,2],two_points_rot[,1])

# Distance between two points should be preserved during rotation
p1 <- two_points[1,]; p2 <- two_points[2,]
theta <- acos( sum(p1*p2) / ( sqrt(sum(p1 * p1)) * sqrt(sum(p2 * p2)) ) )
all.equal(two_points_rot_long[2] - two_points_rot_long[1], theta) # Checks out
# Distance between two points should agree with distGeo
distGeo(matrix(c(1,0.5), ncol = 2), matrix(c(0,1), ncol = 2)) # DOES NOT AGREE!! ERROR IN DIST GEO!!
# Luckily distGeo is not needed anymore:
# NEW BOUNDARY POINT FINDER FUNCTION
# BOUNDARY POINT FINDER FUNCTION
# Works by rotating the two vertices of each line to the equator before sampling, then rotating back
polygon_boundary <- function( vertices_xyz, eps = 0.001 )
{
  rows <- nrow(vertices_xyz)
  if(rows < 3){
    stop("Polygon must have at least 3 vertices.");
  }
  #vertices_xyz_cycled <- rbind(vertices_xyz[2:nrow(vertices_xyz),], vertices_xyz[1,])

  boundary_points <- data.frame()
  for ( row in 1:rows )
  {
    cat(paste("row = ", row))
    two_points <- rbind(vertices_xyz[row,], vertices_xyz[1 + (row %% rows),])
    normal_vector <- vector_cross(two_points[1,], two_points[2,])
    two_points_rotated <- rodrigues(as.matrix(normal_vector)[1,], c(0,0,1), two_points)
    two_longitudes <- atan2(two_points_rotated[,2], two_points_rotated[,1]) # latitudes are both 0

    line_longitudes <- seq(two_longitudes[1], two_longitudes[2], by = eps)
    line <- as.data.frame( sph2car( long = line_longitudes,
                                    theta = rep(0,length(line_longitudes)),
                                    deg = FALSE ) )

    line_rotated_back <-  as.data.frame(rodrigues(c(0,0,1), as.matrix(normal_vector)[1,], line))
    names(line_rotated_back) <- c("x","y","z")
    boundary_points <- rbind(boundary_points, line_rotated_back)
  }

  boundary_points
}

boundary_points <- data.frame()

vertices_xyz2 <- as.data.frame(sph2car(polygon_sph2$long, polygon_sph2$theta, deg = FALSE))
plot3d(vertices_xyz2, col = 'red', type = 'p', size = 20, pch = 3, add = TRUE)
polygon_boundary_points <- polygon_boundary( vertices_xyz2 )
plot3d(polygon_boundary_points, col = 'red', type = 'p', size = 5, pch = 3, add = TRUE)

rows <- nrow(vertices_xyz2)
two_points <- rbind(vertices_xyz2[2,], vertices_xyz2[1 + (2 %% rows),])
normal_vector <- vector_cross(two_points[1,], two_points[2,])
two_points_rotated <- rodrigues(as.matrix(normal_vector)[1,], c(0,0,1), two_points)
two_longitudes <- atan2(two_points_rotated[,2], two_points_rotated[,1]) # latitudes are both 0
line_longitudes <- seq(two_longitudes[1], two_longitudes[2], by = eps)

line <- sph2car( long = line_longitudes,
                        theta = rep(0,length(line_longitudes)),
                        deg = FALSE )
line_rotated_back <-  as.data.frame(rodrigues(c(0,0,1), as.matrix(normal_vector)[1,], line))
names(line_rotated_back) <- c("x","y","z")

boundary_points <- rbind(boundary_points, line_rotated_back)

# a <- normal_vector
# b <- c(0,0,1)
# norm_a <- sqrt(sum(a * a))
# norm_b <- sqrt(sum(b * b))
# !isTRUE( all.equal(a/norm_a, b/norm_b, check.attributes = FALSE) )
# {
#   k <- vector_cross(a,b)
#   k <- k/sqrt(sum(k^2)) # normalised k
#   K <- matrix( c( 0   , -k[3], k[2],
#                   k[3], 0    , -k[1],
#                   -k[2], k[1] , 0    ), nrow = 3, byrow = TRUE)
#   theta <- acos( sum(a*b) / ( norm_a*norm_b ) ) # The angle between a and b
#   I <- diag(c(1,1,1))
#   R <- I + sin(theta)*K + (1-cos(theta))*K%*%K #Rodrigues' Formula.
#   p_xyz <- t(R%*%t(p_xyz))
# }

p_xyz

two_points_rotated <- rodrigues(normal_vector, c(0,0,1), two_points)





# TRY THE CMBDataFrame CONSTRUCTOR and plotCMB function ------------------------------------

# Method 1: Read the data from FITS and construct the CMBDataFrame
df <- CMBDataFrame("../CMB_map_smica1024.fits")

# Method 2: Read the data first, then construct the CMBDataFrame
cmbdat <- readFITScmb("../CMB_map_smica1024.fits")
# We'll just take the sample pixels 1,2 and 3
df <- CMBDataFrame(CMBData = cmbdat, spix = c(2,4,6))
# Alternatively specify a sample size for a random sample
df2 <- CMBDataFrame(CMBData = cmbdat, sampleSize = 800000)
plotCMB(df2)

#Plotting with sample pixels
Nside <- 1024     # specify the Nside parameter of the CMB map
N <- 12*256^2     # specify the number of sample pixels
spix <- sample(seq(1,12*Nside^2),N)
df <- CMBDataFrame(CMBData = cmbdat, spix = spix)
#Now make the plot:
sm <- matrix(c(df$long, df$theta, rep(1,N)), nrow = N)
smx <- sph2car(sm, deg = FALSE)
mat <- readMat("cmbmap.mat")
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])
cols <- colmap[cut(df$I,length(colmap))]
open3d()
bg3d("black")
plot3d(smx, col = cols, type = "p", cex = 5, pch = 3, add = TRUE)






# TRY THE GEOSPHERE PACKAGE  ----------------------------------------------
LA <- c(-118.40, 33.95)
NY <- c(-73.78, 40.63)
MELB <- c(144.96,-37.81)
data(wrld)
plot(wrld, type='l')
gc <- greatCircle(LA, NY)
gc2 <- greatCircle(MELB,NY)
lines(gc, lwd=2, col='blue')
lines(gc2, lwd=2, col='blue')
gci <- gcIntermediate(LA, NY)
gci2 <- gcIntermediate(MELB, NY, breakAtDateLine = TRUE)
lines(gci, lwd=4, col='green')
lines(gci2[[1]], lwd=4, col='red')
lines(gci2[[2]], lwd=4, col='red')
points(rbind(LA, NY), col='red', pch=20, cex=2)

mp <- midPoint(MELB, NY)
points(mp, col = 'blue', pch=20, cex=2)














# TRY THE distGeo FUNCTION ------------------------------------------------

Nside <- 256
sph <- pix2angC(Nside)
lon <- sph[,2]
lat <- pi/2 - sph[,1]
dists <- distGeo(matrix(c(0,0), nrow = 1), cbind(lon, lat))


LA <- c(-118.40, 33.95)
NY <- c(-73.78, 40.63)
onCirc <- sph[which(onGreatCircle(LA,NY, sph[,1:2])),][,1:2]

Npix <- length(onCirc)/2
m <- matrix(c(onCirc[,2], pi/2 - onCirc[,1], rep(1,Npix)), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
plot3d(mx, col = "blue", type = 'p', cex = 2, add = TRUE)



# Run tests
# library(testthat)
# testthat::context("Convert HEALPix to spherical coordinates and (j,i) indices")
#
# testthat::test_that("Output matrix is as expected for RING ordering with small Nside", {
#   testthat::expect_equal_to_reference(pix2angC(2,FALSE), "../tests/testthat/references/pix2angRING_02.rds")
#   testthat::expect_equal_to_reference(pix2angC(4,FALSE), "../tests/testthat/references/pix2angRING_04.rds")
#   testthat::expect_equal_to_reference(pix2angC(8,FALSE), "../tests/testthat/references/pix2angRING_08.rds")
#   testthat::expect_equal_to_reference(pix2angC(16,FALSE), "../tests/testthat/references/pix2angRING_16.rds")
#   testthat::expect_equal_to_reference(pix2angC(32,FALSE), "../tests/testthat/references/pix2angRING_32.rds")
# })
#
# testthat::test_that("Output matrix is as expected for NEST ordering with small Nside", {
#   testthat::expect_equal_to_reference(pix2angC(2,TRUE), "../tests/testthat/references/pix2angNEST_02.rds")
#   testthat::expect_equal_to_reference(pix2angC(4,TRUE), "../tests/testthat/references/pix2angNEST_04.rds")
#   testthat::expect_equal_to_reference(pix2angC(8,TRUE), "../tests/testthat/references/pix2angNEST_08.rds")
#   testthat::expect_equal_to_reference(pix2angC(16,TRUE), "../tests/testthat/references/pix2angNEST_16.rds")
#   testthat::expect_equal_to_reference(pix2angC(32,TRUE), "../tests/testthat/references/pix2angNEST_32.rds")
# })
#
# testthat::test_that("Sample pixel results agree with full pixel results with small Nside", {
#   for (Ns in c(2,4,8,16,32)){
#     for (i in c(1,20,12*Ns^2)){
#       eval(bquote(testthat::expect_equal(pix2angC(.(Ns),TRUE,.(i))[1,],pix2angC(.(Ns),TRUE)[.(i),])))
#       eval(bquote(testthat::expect_equal(pix2angC(.(Ns),FALSE,.(i))[1,],pix2angC(.(Ns),FALSE)[.(i),])))
#     }
#   }
# })


# TRY OUT THE NEW readFITScmb FUNCTION ------------------------------------
# timeT <- system.time(
# cmbdatT <- readFITS("CMB_map_smica1024.fits")
# )
time <- system.time(
  cmbdat <- readFITScmb("../CMB_map_smica1024.fits")
)
# Check first 10 intensity values
# head(cmbdat$col[[1]], n = 10)
# Check last 10 intensity values
# tail(cmbdat$col[[1]], n = 10)
# TRY OUT THE SHRUNKEN CMB TEST MAP
cmbTestMap <- readFITScmb("CMB_testmap_1024_10cols.fits")
# Verify shrunken map matches full map
# for (i in 1:5){
#   cat(all.equal(head(cmbdat$col[[i]], n = 10), cmbTestMap$col[[i]]))
#   cat(" ")
# }

# TESTING RING ORDERING -------------------------------------------
Nside <- 2
Npix <- 12*Nside^2
angSphereC <- pix2angC(Nside, Nest = FALSE)

# Plots of result:
m <- matrix(c(angSphereC[,2], pi/2 - angSphereC[,1], rep(1,Npix)), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
plot3d(mx, col = "blue", type = 'l', cex = 1)
plot3d(mx, col = "blue", type = 'p', cex = 1)





# TESTING NESTED ORDERING -------------------------------------------------
Nside <- 2
Npix <- 12*Nside^2
sph <- pix2angC(Nside)

#spho <- sph[order(sph[,3],sph[,4]),]
spho <- sph[order(sph[,3],sph[,4]),]

m <- matrix(c(long = spho[,2],
              lat = pi/2 - spho[,1],
              r = rep(1,Npix)), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
plot3d(mx, col = "blue", type = 'l', cex = 1)
plot3d(mx, col = "red", cex = 30, pch = 3)






# TEST NEST PLOT WITH COLOUR -----------------------------------------------
#library(plotrix)
tNside <- 32
tNpix <- 12*tNside^2
tI <- seq(-1e5,1e5,length.out = tNpix) #rnorm(Npix,0,1e5)
#cols <- color.scale(I, cs1=seq(1,0,0), cs2=seq(1,1,0), cs3=c(0,1,0)) #cs1=1, cs2=1, cs3=c(0,1,0)
tpalette <- colorRampPalette(c("blue","red"))(tNpix)
tcols <- tpalette[cut(tI, tNpix)]
plot(tI, col = tcols)

tsph <- pix2angC(tNside, TRUE)
#spho <- sph[order(sph[,3],sph[,4]),]
tm <- matrix(c(long = tsph[,2], lat = pi/2 - tsph[,1], r = rep(1,tNpix)), nrow = tNpix)
tmx <- sph2car(tm, deg = FALSE)
plot3d(tmx, col = tcols, cex = 30, pch = 3)
plot3d(tmx, col = "blue", type = 'l', cex = 1)
plot3d(tmx, col = tcols, type = 'p', cex = 1)
plot3d(tmx, col = tcols, type = 'h', cex = 1)
plot3d(tmx, col = tcols, type = 's', cex = 1)




# PLOT WITH FULL IMPORTED DATA -------------------------------------
cmbdat <- readFITScmb("CMB_map_smica1024.fits")
Nside <- as.numeric(cmbdat$hdr[81])
Npix <- 12*Nside^2
sph <- pix2angC(Nside)

# Data frame with longitude, latitude and intensity
CMBI <- data.frame(long = sph[,2], lat = pi/2 - sph[,1], I = cmbdat$col[[1]])

# #Full plot, very expensive
# palette <- colorRampPalette(c("blue","red"))(Npix)
# I <- CMBI$I
# cols <- palette[cut(I, Npix)]
# m <- matrix(c(CMBI$long, CMBI$lat, rep(1,Npix)), nrow = Npix)
# mx <- sph2car(m, deg = FALSE)
# plot3d(mx, col = cols, cex = 30, pch = 3, add = TRUE)

# SMALLER PLOT USING RANDOM SAMPLE (QUICKER BUT CLEAR)
# Take a sample from CMBI
sNside <- 256
sNpix <- 12*sNside^2
sample <- sample(seq(1,Npix),sNpix)
sCMBI <- CMBI[sample,]  # can rm(CMBI) now
sm <- matrix(c(sCMBI$long, sCMBI$lat, rep(1,sNpix)), nrow = sNpix)
smx <- sph2car(sm, deg = FALSE)

# Get colour map
mat <- readMat("cmbmap.mat")
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])
# palette <- colorRampPalette(c("blue","green"))(sNpix)
# sCols <- palette[cut(sCMBI$I, sNpix)]

#Plot
cols <- colmap[cut(sCMBI$I,length(colmap))] # partition the colour map
open3d()
bg3d("black")
plot3d(smx, col = cols, type = "p", cex = 5, pch = 3, add = TRUE)
# plot3d(smx, col = sCols, type = "s", cex = 5, pch = 3, add = TRUE)

# PLOT JUST ONE BASE PIXEL
Nbase <- Npix/12
bCMBI <- CMBI[seq(1,Nbase),]
bm <- matrix(c(bCMBI$long, bCMBI$lat, rep(1,Nbase)), nrow = Nbase)
bmx <- sph2car(bm, deg = FALSE)

#Plot
cols <- colmap[cut(bCMBI$I,length(colmap))]
open3d()
bg3d("black")
plot3d(bmx, col = cols, type = "p", cex = 5, pch = 3, add = TRUE)





# TRY AN ALTERNATIVE COLOUR MAP -------------------------------------------
library(colorspace)
#altColmap <- diverge_hcl(256)
#altColmap <- sequential_hcl(256, c = 0, power = 2.2)
#altColmap <- heat_hcl(256, c = c(80, 30), l = c(30, 90), power = c(1/5, 1.3))
altColmap <- rainbow_hcl(256, start = 1, end = 240)

#Plot
altCols <- altColmap[cut(sCMBI$I,length(altColmap))] # partition the colour map
open3d()
bg3d("black")
plot3d(smx, col = altCols, type = "p", cex = 5, pch = 3, add = TRUE)
# plot3d(smx, col = sCols, type = "s", cex = 5, pch = 3, add = TRUE)




# IMPORT/PLOT THE SIMPLE RANDOM SAMPLE CMB MAP ----------------------------
# IN FUTURE WE COULD PASS THE INDICES INTO A SAMPLE CMB MAP COLUMN IN PYTHON.
# FOR NOW WE IMPORT THE INDICES FROM A .CSV
sCMB <- readFITScmb("CMB_testmap_1024_256sample.fits")
spix <- read.table("exploration/CMB_testmap_1024_256sampleIndices.csv", sep = ",")[,1]
Nside <- as.numeric(sCMB$hdr[51])
sph <- pix2angC(Nside, spix = spix)
# We have to reorder the data to match the sample, which is annoying
sph <- sph[match(spix, sph[,5]),]
# Take a look at the sample pixels with spherical coords
head(sph, n = 10)
tail(sph, n = 10)
# Plotting related
sCMBI <- data.frame(long = sph[,2], lat = pi/2 - sph[,1], I = sCMB$col[[1]])
sm <- matrix(c(sCMBI$long, sCMBI$lat, rep(1,sNpix)), nrow = sNpix)
smx <- sph2car(sm, deg = FALSE)
mat <- readMat("exploration/cmbmap.mat")
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])
cols <- colmap[cut(sCMBI$I,length(colmap))]
open3d()
bg3d("black")
plot3d(smx, col = cols, type = "p", cex = 5, pch = 3, add = TRUE)
