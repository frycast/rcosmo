rm(list = ls())
# Used libraries and source functions
library(rgl)
library(sphereplot)
source("exploration/HealpixBoundaries.R")

# Parameter setting
Nside <- 8
interval <- pi/90

## HELPER FUNCTION 1
plotPixel <- function(S)
{
  C <- rcosmo::sph2car(S)
  rgl::plot3d(C[,"x"],C[,"y"],C[,"z"], col = "black", type = "l", add=TRUE)
}




#### Part I
for (k in seq(1,Nside,1)) {

  start_phi <- pi*k/(2*Nside)
  end_phi <- pi/2
  S <- pbPolarCap(Nside, k, interval, start_phi, end_phi, 0)
  for ( m in seq(0,3,1) ){
    #Northern Hemisphere
    plotPixel( data.frame(theta = S[,1], phi = S[,2] + pi*m/2) )
    #Southern Hemisphere
    plotPixel( data.frame(theta = S[,1] - pi, phi = S[,2] + pi*m/2) )
  }


  start_phi<- 0
  end_phi <- pi/2-pi*k/(2*Nside)
  S <- pbPolarCap(Nside, k, interval, start_phi, end_phi, 1)
  for ( m in seq(0,3,1) ){
    #Northern Hemisphere
    plotPixel( data.frame(theta = S[,1], phi = S[,2] + pi*m/2) )
    #Southern Hemisphere
    plotPixel( data.frame(theta = S[,1] - pi, phi = S[,2] + pi*m/2) )
  }
}


#### Part II#
for (k in seq(0,3,1)){
  #Northern Hemisphere
  PHI<- seq(0,acos(2/3),interval)
  if (PHI[length(PHI)]!= acos(2/3))  PHI = c(PHI, acos(2/3))
  S<-cbind(PHI, rep(1,length(PHI))*k*pi/2)
  plotPixel( data.frame(theta = S[,1], phi = S[,2]) )

  #Southern Hemisphere
  PHI<- seq(acos(-2/3),pi,interval)
  if (PHI[length(PHI)]!= pi)  PHI = c(PHI, pi)
  S<-cbind(PHI, rep(1,length(PHI))*k*pi/2)
  plotPixel( data.frame(theta = S[,1], phi = S[,2]) )

}


#### Part III: Equatorial Belt Area
for (k in seq(-3*Nside,Nside-1,1)){
  start_theta<-acos(-2/3)
  start_phi<-(-4/3+4*k/(3*Nside))*3*pi/8
  end_phi<- 4*k/(3*Nside)*3*pi/8
  S<-pbEqBelt(Nside, k, interval, start_phi, end_phi, start_theta)
  plotPixel( data.frame(theta = S[,1], phi = S[,2]) )

  start_theta<- acos(2/3)
  temp<- -start_phi
  start_phi<- -end_phi
  end_phi<- temp
  aa=k+13
  S<-pbEqBelt(Nside, k, interval, start_phi, end_phi, start_theta)
  plotPixel( data.frame(theta = S[,1], phi = S[,2]) )
}

