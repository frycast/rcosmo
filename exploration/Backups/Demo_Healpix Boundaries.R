rm(list = ls())
# Used libraries and source functions
library(rgl)
library(sphereplot)
source("exploration/HealpixBoundaries.R")
# Parameter setting
Nside <- 4
interval <- pi/90

# Part I#
for (k in seq(1,Nside,1)){
  start_phi<- pi*k/(2*Nside)
  end_phi<- pi/2
  S<-PixBoundariesPolarCap(Nside, k, interval, start_phi, end_phi, 0)
  for (m in seq(0,3,1)){
  #Northern Hemisphere
  C<-sph2car(S[,2]+pi*m/2,pi/2-S[, 1], deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l", add=TRUE)

  #Southern Hemisphere
  C<-sph2car(S[,2]+pi*m/2,S[, 1]-pi/2, deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l", add=TRUE)

  }

  start_phi<- 0
  end_phi<- pi/2-pi*k/(2*Nside)
  S<-PixBoundariesPolarCap(Nside, k, interval, start_phi, end_phi, 1)
  for (m in seq(0,3,1)){
    #Northern Hemisphere
    C<-sph2car(S[, 2]+pi*m/2,pi/2-S[, 1], deg = FALSE)
    plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)

    #Southern Hemisphere
    C<-sph2car(S[,2]+pi*m/2, S[, 1]-pi/2,deg = FALSE)
    plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)
  }
}


# Part II#
for (k in seq(0,3,1)){
  #Northern Hemisphere
  PHI<- seq(0,acos(2/3),interval)
  if (PHI[length(PHI)]!= acos(2/3))  PHI = c(PHI, acos(2/3))
  S<-cbind(PHI, rep(1,length(PHI))*k*pi/2)
  C<-sph2car(S[,2],pi/2-S[,1], deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)

  #Southern Hemisphere
  PHI<- seq(acos(-2/3),pi,interval)
  if (PHI[length(PHI)]!= pi)  PHI = c(PHI, pi)
  S<-cbind(PHI, rep(1,length(PHI))*k*pi/2)
  C<- sph2car(S[,2],pi/2-S[,1], deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)

}


## Part III: Equatorial Belt Area
for (k in seq(-3*Nside,Nside-1,1)){
  start_theta<-acos(-2/3)
  start_phi<-(-4/3+4*k/(3*Nside))*3*pi/8
  end_phi<- 4*k/(3*Nside)*3*pi/8
  S<-PixBoundariesEquatorialBelt(Nside, k, interval, start_phi, end_phi, start_theta)
  C<-sph2car(S[,2],pi/2-S[,1], deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)


  start_theta<- acos(2/3)
  temp<- -start_phi
  start_phi<- -end_phi
  end_phi<- temp
  aa=k+13
  S<-PixBoundariesEquatorialBelt(Nside, k, interval, start_phi, end_phi, start_theta)
  C<-sph2car(S[,2],pi/2-S[,1], deg = FALSE)
  plot3d(C[,1],C[,2],C[,3], col = "black", type = "l",add=TRUE)
}

# #######To reflect the boundaries with the healpix points
# Npix <- 12*Nside^2
# angSphereC <- pix2angC(Nside, FALSE)
#
# m <- matrix(c(angSphereC[,2], pi/2 - angSphereC[,1], rep(1,Npix)), nrow = Npix)
# mx <- sph2car(m, deg = FALSE)
# plot3d(mx, col = "blue", type = 'p', pch=16, add = TRUE)
# spheres3d(0,0,0,lit=FALSE,color="green")
