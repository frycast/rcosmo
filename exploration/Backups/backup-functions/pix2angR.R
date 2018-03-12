rm(list=ls())
Nside <- 10
Nrings <- 4*Nside - 1                    # Number of isolatitude rings.
Npix <- 12*Nside^2                         # Total number of pixels.
NpixPC <- sum(4*seq(1,Nside-1))       # Numer of pixels in each polar cap region
NringsPC <- Nside - 1                    # Number of rings in each polar cap region
pixArea <- 4*pi/Npix                       # Area of each pixel (unit sphere).

# There are 2Nside - 1 rings in the equatorial zone.
# There are 4Nside - 1 rings in total.
# Hence there are Nside rings in each polar zone.
# Eq-zone rings have 4Nside pixels each.
# Num polar zone rings = Nside closest to eq-zone and decreases by 4 pix per ring.
# Hence:
pixPerRing <- c(4*seq(1,Nside),          # Stores number of pixels in ring i.
                rep(4*Nside,2*Nside - 1),  
                4*rev(seq(1,Nside)))     # Check: length(pixPerRing) == Nrings
cPixPerRing <- cumsum(pixPerRing)        # Cumulative sum of pix per ring.

# Isolatitude-ring indices for North Polar, North Equatorial, South Equatorial and South Polar regions.
NPindex <- c(1,Nside - 1)
NEindex <- c(Nside,2*Nside)
SEindex <- c(2*Nside+1,3*Nside)
SPindex <- c(3*Nside+1,4*Nside - 1)

# Pixel indices for regional boundaries with RING ordering scheme
bpiRingNP <- c(0,cPixPerRing[NPindex[2]]-1)
bpiRingNE <- c(cPixPerRing[NPindex[2]],cPixPerRing[NEindex[2]] - 1) 
bpiRingSE <- c(cPixPerRing[NEindex[2]],cPixPerRing[SEindex[2]] - 1)
bpiRingSP <- c(cPixPerRing[SEindex[2]],cPixPerRing[SPindex[2]] - 1)

##CHECK: e.g. for Nside = 3:
#        bpiRingNP  # = 0 3
#        bpiRingNE  # = 12 59
#        bpiRingSE  # = 60 95
#        bpiRingSP  # = 96 107








# HELPER FUNCTIONS -------------------------------------------------------------------------


######## NORTH POLAR CAP  (1 <= i < Nside) 

# i = ring index (north->south), 
# p = pixel index (starts at 0).
p2iNP <- function(p){
  ph <- (p+1)/2
  floor(sqrt(ph-sqrt(floor(ph)))) + 1
}

# j = pixel-in-ring index.
p2jNP <- function(p){
  i <- p2iNP(p)
  p + 1 - 2*i*(i-1)
}

# z = cos(theta), theta = colatitude in [0,pi] (north-south)
p2zNP <- function(p){
  i <- p2iNP(p)
  1 - i^2/(3*Nside^2) 
}

# phi = eastward longitude in [0,2pi]
p2phiNP <- function(p){
  i <- p2iNP(p)
  j <- p2jNP(p)
  pi/(2*i)*(j - 1/2) 
}

##I DONT THINK THESE CHECKS ARE WORKING:
##CHECK: p2jNP(sum(pixPerRing[1:Nside]) - 1) = 4*Nside
##CHECK: p2jNP(sum(pixPerRing[1:Nside-1]) - 1) = 4*Nside - 4
##CHECK: dp <- p2iNP(p + pixPerRing[p2iNP(p)] - p2jNP(p) + 1) 
#        # should return i + 1 where p lies on ring i, that is:
#        p2iNP(dp) - p2iNP(p) = 1.

##CHECK:  print ring and pixel in ring indices for whole region
#         for (p in 0:(sum(pixPerRing[1:(Nside-1)])-1) ){
#           cat(p2iNP(p))
#           cat(" ")
#           cat(p2jNP(p))
#           cat("\n")
#         }




######## NORTH EQUATORIAL BELT (Nside <= i <= 2Nside)

p2iNE <- function(p){
  pp <- p - 2*Nside*(Nside - 1)
  floor(pp/(4*Nside)) + Nside
}

p2jNE <- function(p){
  pp <- p - 2*Nside*(Nside - 1)
  pp %% (4*Nside) + 1
}

p2zNE <- function(p){
  i <- p2iNE(p)
  4/3 - 2*i/(3*Nside)
}

p2phiNE <- function(p){
  i <- p2iNE(p)
  j <- p2jNE(p)
  s <- (i - Nside + 1) %% 2
  pi/(2*Nside)*(j - s/2)
}

##CHECK: print ring and pixel in ring indices for whole region
#         for (p in (sum(pixPerRing[1:(Nside-1)]))
#                  :(sum(pixPerRing[1:(2*Nside)])-1) ){
#          cat(p2iNE(p))
#          cat(" ")
#          cat(p2jNE(p))
#          cat("\n")
#         }

##CHECK: The area of a grid element, for p in equatorial zone:
#         p <- sum(pixPerRing[1:Nside])                  # gives first pixel in ring i = Nside + 1 
#         dp <- p + pixPerRing[p2iNE(p)] - p2jNE(p) + 1  # gives first pixel on ring i + 1
#         dz <- p2zNE(dp) - p2zNE(p)                     # increment z by 1 grid unit
#         dphi <- p2phiNE(p + 1) - p2phiNE(p)            # increment phi by 1 grid unit
#         area <- abs(dphi*dz)
#         all.equal(area, pixArea)






######## SOUTH EQUATORIAL BELT (2Nside < i <= 3Nside)

# p2i, p2j, p2phi, p2z should be unchanged from NORTH EQUATORIAL BELT

p2iSE <- p2iNE
p2jSE <- p2jNE
p2phiSE <- p2phiNE
p2zSE <- p2zNE

##CHECK: print ring and pixel in ring indices for whole region
#         for (p in (sum(pixPerRing[1:(2*Nside)]))
#                  :(sum(pixPerRing[1:(3*Nside)])-1) ){
#                cat(p2iNE(p))
#                cat(" ")
#                cat(p2jNE(p))
#                cat("\n")
# }

##CHECK: 
#         p <- sum(pixPerRing[1:(2*Nside)])                  # gives first pixel in ring i = 2*Nside + 1 
#         dp <- p + pixPerRing[p2iNE(p)] - p2jNE(p) + 1      # gives first pixel on ring i + 1
#         dz <- p2zSE(dp) - p2zSE(p)                         # increment z by 1 grid unit
#         dphi <- p2phiNE(p + 1) - p2phiNE(p)                # increment phi by 1 grid unit
#         area <- abs(dphi*dz)
#         all.equal(area, pixArea)






######## SOUTH POLAR CAP (3Nside < i <= 4Nside-1)

# CURRENTLY THESE ARE STUFFED UP:
p2iSP <- function(p){
  ps <- Npix - p - 1
  4*Nside - p2iNP(ps)
}
p2jSP <- function(p){
  ps <- Npix - p - 1
  i <- p2iSP(p)
  pixPerRing[i] - p2jNP(ps) + 1
}

##CHECK:  e.g. for Nside = 3 we have NpixPC = 108
#          look at cumsum(pixPerRing) and note that 96 is
#          the index of the first south polar cap pixel:
#             p2iSP(96)        # should give i = 10.
#             p2iSP(104)       # should give i = 11.
#             p2jSP(96)        # should give j = 1.
#             p2jSP(104)       # should give j = 1.
#          In general see:
#          for (p in sum(pixPerRing[1:(3*Nside)]):(Npix-1)){
#             cat(p2iSP(p))
#             cat(" ")
#             cat(p2jSP(p))
#             cat("\n")
#          }



p2zSP <- function(p)
{
  i <- p2iSP(p)
  is <- Nrings - p2iSP(p) + 1
  is^2/(3*Nside^2) - 1 
}


p2phiSP <- function(p)
{
  i <- p2iSP(p)
  is <- Nrings - p2iSP(p) + 1
  j <- p2jSP(p)
  pi/(2*is)*(j - 1/2) 
}








# MAIN FUNCTIONS ----------------------------------------------------------

pix2ang <- function(p){
  if (bpiRingNP[1] <= p & p <= bpiRingNP[2]) {
    
    return(c(acos(p2zNP(p)), p2phiNP(p)))
    
  } else if (bpiRingNE[1] <= p & p <= bpiRingNE[2]) {
    
    return(c(acos(p2zNE(p)), p2phiNE(p)))
    
  } else if (bpiRingSE[1] <= p & p <= bpiRingSE[2]) {
    
    return(c(acos(p2zSE(p)), p2phiSE(p)))
    
  } else {
    
    return(c(acos(p2zSP(p)), p2phiSP(p)))
    
  }
}

pix2angz <- function(p){
  if (bpiRingNP[1] <= p & p <= bpiRingNP[2]) {
    
    return(c(p2zNP(p), p2phiNP(p)))
    
  } else if (bpiRingNE[1] <= p & p <= bpiRingNE[2]) {
    
    return(c(p2zNE(p), p2phiNE(p)))
    
  } else if (bpiRingSE[1] <= p & p <= bpiRingSE[2]) {
    
    return(c(p2zSE(p), p2phiSE(p)))
    
  } else {
    
    return(c(p2zSP(p), p2phiSP(p)))
    
  }
}


# CURRENTLY USELESS
ang2xyz <- function(theta, phi){
  list(x = cos(phi)*sin(theta), 
       y = sin(theta)*sin(phi), 
       z = cos(theta))
}


# TESTING -----------------------------------------------------------------

testSphere <- seq(0,Npix - 1)
angSphere <- t(sapply(testSphere, pix2ang))


library(sphereplot)
m <- matrix(c(angSphere[,2], pi/2 - angSphere[,1], rep(1,Npix)), nrow = Npix)
mx <- sph2car(m, deg = FALSE)
plot3d(mx, col = "blue", type = 'l', cex = 1)
plot3d(mx, col = "red", cex = 30, pch = 3, add = TRUE)

