#' (old/unworking) Covariance for CMB
#'
#' This function computes the covariance for CMB.
#' It does not place data into bins, but instead generates
#' points on a circle of radius r and then finds the closest
#' HEALPix point to each, using nestSearch.
#'
#' @param rmin,rmax are the minimum and maximum of radii.
#' @param Nr is the number of radii in between rmin and rmax for which covariance of CMB is evaluated.
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
covCMBold <- function(df = CMBDataFrame(CMBData = "CMB_map_smica1024.fits",
                                     coords = "HEALPix",
                                     ordering = "nested"),
                   rmin = 10^(-6), rmax = pi-10^(-6),
                   Nr = 10, Nside = 1024, N_x_vec = 10)
{

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
  Pix_x <- 12*Nside^2 # number of HEALPix points
  pix_x_vec <- floor(seq(1,Pix_x,length.out = N_x_vec))
  print(paste0("  ** Number of x: N_x  = ",N_x_vec))
  ###

  ### ??? WHY CHOOSE N_y TO BE THIS ???? (USED TO BE USED BELOW IN poinOnCircle(Nside = N_y))
  # N_y  <- Nside # number of points N_y = 2^n and N_y ~ sqrt(N_x_vec) ######
  N_y  <- 2^(ceiling(log2((sqrt(N_x_vec)))))
  print(paste("  ** Number of points on each circle: N_y  =", N_y))

  ###
  print(paste0("  ** Nside for y: Nside  = ", Nside))
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

      for(i_pix_x  in 1:N_x_vec)
      {
        pix_x  <- pix_x_vec[i_pix_x]
        x  <- pix2vec(Nside = Nside, order="ring",Pix=pix_x)

        T_x  <- T_I[pix_x] # T at x

        x_sph_1 <- ca2sph(matrix(c(x$x,x$y,x$z), ncol = 3, byrow = TRUE))
        x_sph <- c(x_sph_1$theta,x_sph_1$phi)
        y <- pointOnCircle(Nside = Nside, radius = r, center = x_sph)
        # cat(paste("\ni_pix_x: ", i_pix_x, "\n"))
        # cat("---------------------------------------------------------\n")
        # cat(paste("\nOn circle    : ", round(y,10), "\nCircle center: ", round(x,10),"\n"))
        # cat("\n---------------------------------------------------------\n")

        #     pix_y  <- vector(mode = "numeric", nrow(y))
        # # compute the closest HEALPix points to each points of x_2
        #     for(i in seq_along(pix_y)) {
        #       pix_y[i]  <- nestSearch(target = y[i,], Nside = Nside, index.only = TRUE)
        #     }
        # THIS LINE DOES ALL OF THE ABOVE IN ONE LINE
        pix_y <- apply(y, MARGIN = 1, nestSearch, Nside = Nside, index.only = TRUE)

        # comput the sum of temperature intensity over the HEALPix on the circle
        # with center y and radius r
        T_y  <- sum(T_I[pix_y])

        T_yx[i_pix_x]  <- T_y*T_x
      }

      # compute the average over all points x
      T_r[r_i]  <- sum(T_yx)/N_x_vec

    } else {
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
