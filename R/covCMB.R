#' Covariance for CMB
#'
#' The function \code{covCMB} computes the covariance for CMB.
#'
#' @param rmin,rmax are the minimum and maximum of radii.
#' @param Nr is the number of radii in between rmin and rmax for which covariance of CMB is evalueated.
#' @param Nside is the Nside for which the HEALPix points are used.
#' @param N_x_vec is the number of points x for each radius, the number of points y for each
#' is 2^(ceil(log2((sqrt(N_x_vec))))) which is equivalent to sqrt(N_x_vec).
#'
#' @return the output is the data frame of radius r and the covariance Tcov.
#'
#' @examples
#' # compute the covariance of CMB at Nside = 1024 and radii between 10^(-6) and pi-10^(-6) with 10 radii
#' covCMB(rmin = 10^(-6), rmax = pi-10^(-6), Nr = 10, Nside = 1024, N_x_vec = 10)
#'
#' @export
covCMB <- function(df = CMBDataFrame(CMBData = "CMB_map_smica1024.fits", coords = "HEALPix", ordering = "nested"), rmin = 10^(-6), rmax = pi-10^(-6), Nr = 10, Nside = 1024, N_x_vec = 10){

# require("rcosmo")
# require("pracma")
# require("spark")
# source("ca2sph.R")
# source("pix2vec.R")
# source("nestSearch.R")
# source("pointOnCircle.R")
# source("pix2vec.R")
# source("nest_search.R")

# # function: linspace
#  linspace <- function(a,b,N){
#    # generate N points between a and b
#    h <- (b-a)/(N-1)
#    lN <- a + h*((1:N)-1)
#    return(lN)
#  }

print(" ****************************************** ")
print("  Compute Covariance Function for CMB ")

print("  *****************************************")
###
r_vec <- seq(rmin,rmax, length.out = Nr)
print(paste0("  ** Radius: min  <- ",rmin,", max  <- ",rmax,", Nr = ", Nr))
###
Nside_x <- Nside # Nside for target points
Pix_x <- 12*Nside_x^2 # number of HEALPix points
pix_x_vec <- floor(seq(1,Pix_x,length.out = N_x_vec))
print(paste0("  ** Number of x: N_x  = ",N_x_vec))
###
# N_y  <- Nside # number of points N_y = 2^n and N_y ~ sqrt(N_x_vec) ######
  N_y  <- 2^(ceil(log2((sqrt(N_x_vec)))))
 print(paste("  ** Number of points on each circle: N_y  =",N_y))
###
Nside_y  <- Nside # Nside for y
 print(paste0("  ** Nside for y: Nside_y  = ",Nside_y))
###
 print("  ****************************************")

T_I  <- df$I

## compute covariance
T_r  <- c(rep(0,Nr)) # covariance

for(r_i in 1:Nr){
  r  <- r_vec[r_i] # radius
  print(paste0("  - compute covariance for radius ",round(r,digits = 5)))

  if (r!=0) {
  T_yx  <- c(rep(0,N_x_vec))
  for(i_pix_x  in 1:N_x_vec) {
    pix_x  <- pix_x_vec[i_pix_x]
    x  <- pix2vec(Nside = Nside_x,order_ring=TRUE,Pix=pix_x)

    T_x  <- T_I[pix_x] # T at x

    x_sph_1 <- ca2sph(matrix(c(x$x,x$y,x$z), ncol = 3, byrow = TRUE))
    x_sph <- c(x_sph_1$theta,x_sph_1$phi)
    y <- pointOnCircle(Nside = N_y, radius = r, center = x_sph)
    # cat(paste("\ni_pix_x: ", i_pix_x, "\n"))
    # cat("---------------------------------------------------------\n")
    # cat(paste("\nOn circle    : ", round(y,10), "\nCircle center: ", round(x,10),"\n"))
    # cat("\n---------------------------------------------------------\n")

    N_y_1 <- dim(y)[1]
    yhp_y  <- c(rep(0,N_y_1))
    pix_y  <- c(rep(0,N_y_1))
# compute the closest HEALPix points to each points of x_2
    for(i in 1:N_y_1) {
      tp  <- y[i,]
      pix_y[i]  <- nestSearch(tp = tp, Nside = Nside_y, index_only = TRUE)
    }

# comput the sum of temperature intensity over the HEALPix on the circle
# with center y and radius r
    T_y  <- sum(T_I[pix_y])

    T_yx[i_pix_x]  <- T_y*T_x
  }

# compute the average over all points x
  T_r[r_i]  <- sum(T_yx)/N_x_vec
  }
  if (r==0) {
    T_x_1 <- T_I[pix_x_vec]
    T_r[r_i]  <- sum(T_x_1*T_x_1)/N_x_vec
  }

  print(paste0("  - Tcov: ",T_r[r_i]))

# print(paste0("  == Running time: #.4f",t))
  print("  ****************************************")
}

Trdf <- data.frame(r = r_vec, Tcov = T_r)
return(Trdf)
}
