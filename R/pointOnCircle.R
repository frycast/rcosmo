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



pointOnCircle <- function(Nside = 16, radius = 0.2, center = c(1,0.7))
{
  # DEFINE MINIMUM SAMPLE DISTANCE
  N <- 5*Nside      # COME UP WITH GOOD FORMULA FOR N (NOT THIS SAMPLE ONE)
  eps <- 2*pi/N

  # MAKE THE SAMPLE POINTS IN SPHERICAL THEN CONVERT TO CARTESIAN
  # AT NORTH POLE (unrotated)
  north.sph <- matrix(c( rep(radius,N),
                         seq(0,2*pi-eps,eps)
                       ),
                      nrow = N, ncol = 2)

  north.xyz <- matrix(c( cos(north.sph[,2])*sin(north.sph[,1]), # Cartesian.
                         sin(north.sph[,2])*sin(north.sph[,1]),
                         rep(cos(north.sph[1,1]),N)
                       ),
                      nrow = N, ncol = 3)

  # TRANSLATE THE SAMPLE POINTS
  # Using Rodrigues' Rotation Formula
  I <- diag(c(1,1,1))
  # 1. Rotate about y axis by theta according to Right Hand Rule.
  Ky <- matrix( c(0,0,1,
                  0,0,0,
                  -1,0,0), nrow = 3, byrow = TRUE)
  Ry <- I + sin(center[1])*Ky + (1-cos(center[1]))*Ky%*%Ky #Rodrigues' Formula.
  py <- t(Ry%*%t(north.xyz))

  # 2. Rotate about z axis by phi according to Right Hand Rule.
  Kz <- matrix( c(0,-1,0,
                  1, 0,0,
                  0, 0,0), nrow = 3, byrow = TRUE)
  Rz <- I + sin(center[2])*Kz + (1-cos(center[2]))*Kz%*%Kz #Rodrigues' Formula.
  pyz <- t(Rz%*%t(py))

  ###-------- JUST FOR VISUALISATION ----------###
  ## Healpix grid for reference
  # library(rgl)
  # Npix <- 12*Nside^2
  # mx <- pix2vec(Nside, order_ring = TRUE)
  # center.xyz <- c(cos(center[2])*sin(center[1]), # Cartesian
  #                 sin(center[2])*sin(center[1]),
  #                 cos(center[1]))
  ## plot
  # rgl::open3d()
  # rgl::bg3d("black")
  # rgl::plot3d(mx, col = "gray", type = 'p', pch = "e", cex = 1)
  # rgl::plot3d(north.xyz, col = "yellow", type = 'p', cex = 1, add = TRUE)
  # rgl::plot3d(center.xyz[1],center.xyz[2],center.xyz[3], col = "red", cex = 1, add = TRUE)
  # rgl::plot3d(-center.xyz[1],-center.xyz[2],-center.xyz[3], col = "red", cex = 1, add = TRUE)
  # C_line <- matrix(c(center.xyz,-center.xyz), nrow = 2, ncol = 3, byrow = TRUE)
  # rgl::plot3d(C_line, type = 'l', col = "red", cex = 1, add = TRUE)
  # rgl::plot3d(py, col = "green", add = TRUE)
  # rgl::plot3d(pyz, col = "purple", add = TRUE)
  ###------------------------------------------###

  # TO DO: NEST SEARCH TO FIND CLOSEST HEALPIX
  return(pyz)
}
