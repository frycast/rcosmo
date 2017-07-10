#rm(list=ls())
library(Rcpp)
library(rgl)
library(sphereplot) # for sph2car()
#help(rgl)
library(FITSio)
library(R.matlab) # for importing the colour map
sourceCpp("exploration/pix2ang.cpp")
source("exploration/readFITScmb.R")

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
angSphereC <- pix2angC(Nside, FALSE)

# Plots of result:
m <- matrix(c(angSphereC[,2], pi/2 - angSphereC[,1], rep(1,Npix)), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
plot3d(mx, col = "blue", type = 'l', cex = 1)
plot3d(mx, col = "blue", type = 'p', cex = 1)





# TESTING NESTED ORDERING -------------------------------------------------
Nside <- 2
Npix <- 12*Nside^2
sph <- pix2angC(Nside, TRUE)

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
sph <- pix2angC(Nside, TRUE)

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
sCMB <- readFITScmb("exploration/CMB_testmap_1024_256sample.fits")
spix <- read.csv("exploration/CMB_testmap_1024_256sampleIndices.csv")[,1]
sNside <- 256
sNpix <- 12*sNside^2
sph <- pix2angC(sNside, TRUE, spix)
CMBI <- data.frame(long = sph[,2], lat = pi/2 - sph[,1], I = cmbdat$col[[1]])
sm <- matrix(c(sCMBI$long, sCMBI$lat, rep(1,sNpix)), nrow = sNpix)
smx <- sph2car(sm, deg = FALSE)
mat <- readMat("cmbmap.mat")
colmap <- rgb(mat$map[,1], mat$map[,2], mat$map[,3])
cols <- colmap[cut(sCMBI$I,length(colmap))]
open3d()
bg3d("black")
plot3d(smx, col = cols, type = "p", cex = 5, pch = 3, add = TRUE)
