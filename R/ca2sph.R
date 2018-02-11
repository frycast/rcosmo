# compute the spherical coordinates (theta,phi) of xyz which is in Cartesian coordinates
# dim(xyz) = [N 3]

#'@export
ca2sph <- function(xyz = matrix(c(0,0,1),ncol=3,byrow=TRUE)) {
  N_xyz <- dim(xyz)[1] # number of points
  theta <- acos(xyz[,3])
  x1 <- xyz[,1]
  x2 <- xyz[,2]
  x3 <- xyz[,3]
  phi <- c(rep(0,N_xyz))
  t <- phi
  # x3==1
  logic_x3 <- x3==1 | x3==-1
  phi[logic_x3] <- 0
  # x2>=0 & -1<x3<1
  logic_x2_1 <- x2>=0 & x3<1 & x3>-1
  t[logic_x2_1] <- x1[logic_x2_1]/sqrt(1-(x3[logic_x2_1])^2)
  t[t>1] <- 1
  t[t<(-1)] <- -1
  phi[logic_x2_1] <- acos(t[logic_x2_1])
  # x2<0 & x3<1
  logic_x2_2 <- x2<0 & x3<1 & x3>-1
  t[logic_x2_2] <- x1[logic_x2_2]/sqrt(1-(x3[logic_x2_2])^2)
  t[t>1] <- 1
  t[t<(-1)] <- -1
  phi[logic_x2_2] <- 2*pi-acos(t[logic_x2_2])

  ca <- data.frame(theta = theta, phi = phi)
  return(ca)
}
