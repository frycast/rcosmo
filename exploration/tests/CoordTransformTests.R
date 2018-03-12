rm(list=ls())
library(Rcpp)
library(rgl)
library(R.matlab)
library(microbenchmark)
library(rcosmo)


##################################################################
### RECENT AND UP-TO DATE TESTS      #############################
##################################################################



## check that rcosmo::pix2coords agrees with python results
res <- vector(mode = "logical", length = 4)

# Compare with Python's results in NESTED CARTESIAN
Dxn <- pix2coords(nside = Nside, nest = TRUE, cartesian = TRUE )[,1:3]
m <- paste0("exploration/tests/verify_py/nest_ns",toString(Nside),".mat")
Pxn <- readMat(m)
Pxn <- t(Pxn$hpt.nest)
res[1] <- all.equal(Pxn, Dxn, check.attributes = FALSE)

# Compare with Python's results in RING CARTESIAN
Dxr <- pix2coords(nside = Nside, nest = FALSE, cartesian = TRUE )[,1:3]
m <- paste0("exploration/tests/verify_py/ring_ns",toString(Nside),".mat")
Pxr <- readMat(m)
Pxr <- t(Pxr$hpt.ring)
res[2] <- all.equal(Pxr, Dxr, check.attributes = FALSE)

# Compare with Python's results in NESTED SPHERICAL
Dsn <- pix2coords(nside = Nside, nest = TRUE, cartesian = FALSE )[,1:2]
Dsn <- as.data.frame(Dsn)
names(Dsn) <- c("theta", "phi")
Dsn <- as.matrix(sph2car(Dsn))
res[3] <- all.equal(Pxn, Dsn, check.attributes = FALSE)

# Compare with Python's results in RING SPHERICAL
Dsr <- pix2coords(nside = Nside, nest = FALSE, cartesian = FALSE )[,1:2]
Dsr <- as.data.frame(Dsr)
names(Dsr) <- c("theta", "phi")
Dsr <- as.matrix(sph2car(Dsr))
res[4] <- all.equal(Pxr, Dsr, check.attributes = FALSE)

# results
res


#### TESTING nest2ring vs nest2ringR (by Yuguang from matlab)

pix <- seq(1,48)
all.equal(as.numeric(nest2ring(2, pix)),
          as.numeric(nest2ringR(2, pix)),
          check.attributes = FALSE)









#######################################################################
################### OLDER OUTDATED TESTS ##############################
#######################################################################

# source("pix2vec.R")
# sourceCpp("pix2ang.cpp")
# sourceCpp("pix2ang_2.cpp")
# sourceCpp("pix2coord.cpp")

### COMPARING WITH NESTED ORDERING -----------------------------------

# Danny's method
p2a <- pix2coords(Nside = Nside, Nest = TRUE, cartesian = TRUE)
D <- p2a[,1:3]

# Yu Guang's method
Y <- as.matrix(pix2vec(Nside, "nested"))

# Find rows where Yu Guang's Y disagrees with Danny's D
eq <- c()
for (i in 1:Npix)
{
  eq[i] <- isTRUE(all.equal(Y[i,],D[i,], check.attributes = FALSE))
}

## Number and ratio of equal and unequal elements
sum(eq)
sum(!eq)
sum(eq)/sum(!eq) # Always 5


# Inspect some values where the rows are not equal
head(Y[!eq,], n = 10) # Y consists of the same (repeated) row when Nside = 2, but not for higher Nside
head(D[!eq,], n = 10) # Suspect the z column of D always matches Y but other columns dont.

# Show that z column is all equal
all.equal(Y[,3],D[,3])    # TRUE

# Suspect that the problem is with the pixel-in-ring indices in pix2coords?

# These are the NESTED isolatitude-ring indices where it doesn't agree
p2a[,4][!eq]

# These are the NESTED pixel in ring indices where it doesn't agree
p2a[,5][!eq]



# COMPARISON REPEATED WITH RING ORDERING (ALSO DOESNT AGREE) ---------------------------------


# Danny's method
p2ar <- pix2coords(Nside = Nside, Nest = TRUE, cartesian = TRUE)
Dr <- p2ar[,1:3]

# Yu Guang's method
Yr <- as.matrix(pix2vec(Nside,"ring"))

# Find rows where Yu Guang's Y1 disagrees with Danny's D
eqr <- c()
for (i in 1:Npix)
{
  eqr[i] <- isTRUE(all.equal(Yr[i,],Dr[i,], check.attributes = FALSE))
}

## Number and ratio of equal and unequal elements
sum(eqr)
sum(!eqr)
sum(eqr)/sum(!eqr) # Many more are unequal

# Inspect some values where the rows are not equal
head(Yr[!eqr,], n = 10)
head(Dr[!eqr,], n = 10)

# This time the z columns are not equal!
all.equal(Yr[,3],Dr[,3])    # FALSE

# These are the pixel in ring indices where it doesn't agree
sort(unique(p2ar[,5][!eqr]))

# These are the RING isolatitude-ring indices where it doesn't agree
sort(unique(p2ar[,4][!eqr]))





# NOW LOOK AT PLOTS FOR COMPARISON ----------------------------------------

# NESTED (positions disagree in polar caps)
open3d()
bg3d("black")
plot3d(Y, type = "p", col = "white", cex = 5, pch = 3, add = TRUE)
plot3d(D, type = "p", col = "red", cex = 5, pch = 3, add = TRUE)

# RING (positions agree everywhere but different order)
open3d()
bg3d("black")
plot3d(Yr, type = "p", col = "white", cex = 5, pch = 3, add = TRUE)
plot3d(Dr, type = "p", col = "red", cex = 5, pch = 3, add = TRUE)




## COMPARE WITH PYTHON HEALPy ----------------------------

# Compare with Python's results in NESTED
rdm <- paste0("verify_py/nest_ns",toString(Nside),".mat")
py <- readMat(rdm)
py <- t(py$hpt.nest)

all.equal(py,Y, check.attributes = FALSE)



# Compare with Python's results in RING
rdm <- paste0("verify_py/ring_ns",toString(Nside),".mat")
pyr <- readMat(rdm)
pyr <- t(pyr$hpt.ring)

all.equal(pyr, Yr, check.attributes = FALSE)




### DOES NEW FUNCTION pix2ang_2 AGREE WITH HEALPy? ---------------------

# NESTED NO!
D2 <- pix2angC_2(4, Nest = TRUE, cartesian = TRUE)[,1:3]
Y <- as.matrix(pix2vec(4,"nested"))

all.equal(D2,Y, check.attributes = FALSE)

# RING YES!
D2r <- pix2angC_2(2, Nest = FALSE, cartesian = TRUE)[,1:3]
Yr <- as.matrix(pix2vec(2,"ring"))

all.equal(D2r,Yr, check.attributes = FALSE)



## DOES NEW FUNCTION pix2coordC AGREE WITH HEALPy? -------------------

sourceCpp("pix2coord.cpp")

# NESTED YES!
D2 <- pix2coordC(16, Nest = TRUE, cartesian = TRUE)[,1:3]
Y <- as.matrix(pix2vec(16,"nested"))

all.equal(D2,Y, check.attributes = FALSE)


# RING YES!
D2r <- pix2coordC(16, Nest = FALSE, cartesian = TRUE)[,1:3]
Yr <- as.matrix(pix2vec(16,"ring"))

all.equal(D2r,Yr, check.attributes = FALSE)

# pix2coordC WINS!!



# TESTING THE NEW nest2ringC -------------------------------

pix <- seq(1,48)
all.equal(as.numeric(nest2ring(2, pix)),
          nest2ringC(2, pix),
          check.attributes = FALSE)





### BENCHMARKS ------------------------------------------------

# pix2coords vs pix2vec NESTED
bn <- microbenchmark(
  pix2coords(Nside = 2, Nest = TRUE, cartesian = TRUE),
  pix2vec(2,"nested")
)
sn <- summary(bn)$mean[2]/summary(bn)$mean[1]
cat(paste("pix2coords is about", round(sn,0),
          "times faster than pix2vec for NESTED\n"))

# pix2coords vs pix2vec RING
br <- microbenchmark(
  pix2coords(Nside = 2,Nest = FALSE, cartesian = TRUE),
  pix2vec(2,"ring")
)
sr <- summary(br)$mean[2]/summary(br)$mean[1]
cat(paste("pix2coords is about", round(sr,0),
          "times faster than pix2vec for RING\n"))


# NEW pix2ang_2 VERSUS pix2coords NESTED
bn <- microbenchmark(
  pix2angC_2(Nside = 4, Nest = TRUE, cartesian = TRUE),
  pix2coords(Nside = 4, Nest = TRUE, cartesian = TRUE)
)
sn <- summary(bn)$mean[2]/summary(bn)$mean[1]
cat(paste("pix2angC_2 is about", round(sn,0),
          "times the speed of pix2coords for NESTED\n"))


# NEW pix2ang_2 VERSUS pix2coords RING
bn <- microbenchmark(
  pix2angC_2(Nside = 4, Nest = FALSE, cartesian = TRUE),
  pix2coords(Nside = 4, Nest = FALSE, cartesian = TRUE)
)
sn <- summary(bn)$mean[2]/summary(bn)$mean[1]
cat(paste("pix2angC_2 is about", round(sn,0),"times the speed of pix2coords for RING\n"))


# nest2ring versus nest2ringC
bp <- microbenchmark(
  nest2ringC(8, pix),
  nest2ring(8, pix)
)
sp <- summary(bp)$mean[2]/summary(bp)$mean[1]
cat(paste("nest2ringC is about", round(sp,0),"times the speed of nest2ring\n"))


# NEW pix2ang_2 VERSUS NEW pix2coordC RING
bn <- microbenchmark(
  pix2angC_2(Nside = 4, Nest = FALSE, cartesian = TRUE),
  pix2coordC(Nside = 4, Nest = FALSE, cartesian = TRUE)
)
sn <- summary(bn)$mean[2]/summary(bn)$mean[1]
cat(paste("pix2angC_2 is about", round(sn,0),"times the speed of pix2coord for RING\n"))


# NEW pix2ang_2 VERSUS NEW pix2coordC NESTED
bn <- microbenchmark(
  pix2coordC(Nside = 4, Nest = TRUE, cartesian = TRUE),
  pix2angC_2(Nside = 4, Nest = TRUE, cartesian = TRUE)
)
sn <- summary(bn)$mean[2]/summary(bn)$mean[1]
cat(paste("pix2coordC is about", round(sn,0),"times the speed of pix2angC_2 for NESTED\n"))
