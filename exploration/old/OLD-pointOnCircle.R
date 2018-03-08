# Inputs: C - circle center (theta1, phi1).
#         r - any radius r in (0,pi)
#         Nside - e.g., 1024, 2048
#
# Properties of a circle on a unit sphere -----------------------
# Suppose circle has center (theta1, phi1)
# and suppose a point (theta2, phi2) lies 
# on the circle. Then,
# Radius of Associated Planar Circle = sin(r) (proof: translate to north pole, then r = theta).
# Circumpherence = 2*pi*sin(r)
# Radius (r) = arccos( cos(theta1)cos(theta2) + sin(theta1)sin(theta2)cos(phi1 - phi2) )
# 
# ---------------------------------------------------------------
#
## ALGORITHM
# 1. Determine a minimum sampling distance, eps, in terms of Nside.
# 2. Generate 2*pi*sin(r)/eps points p on the circle:
#           - Generate them with center at North Pole first
#           - Then rotate center to C using Rodrigues' Formula Rotation Matrix.
# 3. Find the closest healpix point h to each p.
# 4. Return the h in data frame.



library(sphereplot)
library(Rcpp)
source("pix2vec.R")


pointOnCircle <- function(Nside = 16, radius = 0.2, center = c(1,0.7))
{

C_xyz <- c(cos(center[2])*sin(center[1]), # Cartesian
           sin(center[2])*sin(center[1]),
           cos(center[1]))

# DEFINE MINIMUM SAMPLE DISTANCE
N <- 5*Nside      # COME UP WITH GOOD FORMULA FOR N (NOT THIS SAMPLE ONE)           
eps <- 2*pi/N

# MAKE THE SAMPLE POINTS IN SPHERICAL THEN CONVERT TO CARTESIAN
# AT NORTH POLE (unrotated)
p_sph <- matrix(c( rep(radius,N),
                   seq(0,2*pi-eps,eps)
), nrow = N, ncol = 2)
p_xyz <- matrix(c( cos(p_sph[,2])*sin(p_sph[,1]), # Cartesian.
                   sin(p_sph[,2])*sin(p_sph[,1]),
                   rep(cos(p_sph[1,1]),N)
), nrow = N, ncol = 3)

# JUST FOR VISUALISATION
# Healpix grid for reference
Npix <- 12*Nside^2
mx <- pix2vec(Nside, order_ring = TRUE)


# TRANSLATE THE SAMPLE POINTS
# Using Rodrigues' Rotation Formula
I <- diag(c(1,1,1))
# 1. Rotate about y axis by theta according to Right Hand Rule.
Ky <- matrix( c(0,0,1,
                0,0,0,
                -1,0,0), nrow = 3, byrow = TRUE)
Ry <- I + sin(C[1])*Ky + (1-cos(C[1]))*Ky%*%Ky #Rodrigues' Formula.
py <- t(Ry%*%t(p_xyz))

# 2. Rotate about z axis by phi according to Right Hand Rule.
Kz <- matrix( c(0,-1,0,
                1, 0,0,
                0, 0,0), nrow = 3, byrow = TRUE)
Rz <- I + sin(C[2])*Kz + (1-cos(C[2]))*Kz%*%Kz #Rodrigues' Formula.
pyz <- t(Rz%*%t(py))


# TO DO: NEST SEARCH TO FIND CLOSEST HEALPIX
pyz
}
