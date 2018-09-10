
#' Plot HEALPix pixel boundaries
#'
#' Plot the HEALPix pixel boundaries at \code{nside}
#'
#' @param nside the HEALPix nside parameter (integer number \eqn{2^k})
#' @param eps controls the smoothness of the plot, smaller eps
#' implies more samples
#' @param col the colour of plotted boundary lines
#' @param lwd the thickness of the plotted boundary lines
#' @param ordering optionally specify an ordering scheme
#' from which to plot HEALPix pixel numbers. Can be
#' either "ring" or "nested"
#' @param incl.labels If \code{ordering} is specified then
#' this parameter sets the pixel indices that will be
#' displayed (default is all indices at \code{nside})
#' @param nums.col specifies the colour of pixel numbers
#' if \code{ordering} is specified
#' @param nums.size specifies the size of pixel numbers
#' if \code{ordering} is specified
#' @param font A numeric font number from 1 to 5,
#' used if \code{ordering} is specified
#' @param ... arguments passed to \code{rgl::plot3d}
#'
#' @return produces a plot of the HEALPix pixel boundaries
#'
#' @examples
#'
#' displayPixelBoundaries(1, eps = pi/90, col = "red")
#' displayPixelBoundaries(2, eps = pi/90, col = "green")
#'
#' @export
displayPixelBoundaries <- function(nside, eps = pi/90,
                             col = "gray",
                             lwd = 1, ordering,
                             incl.labels = 1:(12*nside^2),
                             nums.col = col, nums.size = 1,
                             font = 2, ...)
{

  ### Part I
  for (k in seq(1,nside,1)) {

    start_phi <- pi*k/(2*nside)
    end_phi <- pi/2
    S <- pbPolarCap(nside, k, eps, start_phi, end_phi, 0)
    for ( m in seq(0,3,1) ){
      #Northern Hemisphere
      plotPixel( data.frame(theta = S[,1], phi = S[,2] + pi*m/2),
                 col = col, lwd = lwd, ...)
      #Southern Hemisphere
      plotPixel( data.frame(theta = S[,1] - pi, phi = S[,2] + pi*m/2),
                 col = col, lwd = lwd, ...)
    }


    start_phi<- 0
    end_phi <- pi/2-pi*k/(2*nside)
    S <- pbPolarCap(nside, k, eps, start_phi, end_phi, 1)
    for ( m in seq(0,3,1) ){
      #Northern Hemisphere
      plotPixel( data.frame(theta = S[,1], phi = S[,2] + pi*m/2),
                 col = col, lwd = lwd, ...)
      #Southern Hemisphere
      plotPixel( data.frame(theta = S[,1] - pi, phi = S[,2] + pi*m/2),
                 col = col, lwd = lwd, ...)
    }
  }

  ### PART II
  for (k in seq(0,3,1)){
    #Northern Hemisphere
    PHI <- seq(0,acos(2/3),eps)
    if (PHI[length(PHI)] != acos(2/3))  PHI = c(PHI, acos(2/3))
    S <- cbind(PHI, rep(1,length(PHI))*k*pi/2)
    plotPixel( data.frame(theta = S[,1], phi = S[,2]),
               col = col, lwd = lwd, ...)

    #Southern Hemisphere
    PHI <- seq(acos(-2/3),pi,eps)
    if (PHI[length(PHI)] != pi)  PHI = c(PHI, pi)
    S <-cbind(PHI, rep(1,length(PHI))*k*pi/2)
    plotPixel( data.frame(theta = S[,1], phi = S[,2]),
               col = col, lwd = lwd, ...)

  }

  ### Part III: Equatorial Belt Area
  for (k in seq(-3*nside,nside-1,1)){
    start_theta <- acos(-2/3)
    start_phi <- (-4/3+4*k/(3*nside))*3*pi/8
    end_phi <- 4*k/(3*nside)*3*pi/8
    S <- pbEqBelt(nside, k, eps, start_phi, end_phi, start_theta)
    plotPixel( data.frame(theta = S[,1], phi = S[,2]),
               col = col, lwd = lwd, ...)

    start_theta <- acos(2/3)
    temp <- -start_phi
    start_phi <- -end_phi
    end_phi <- temp
    S <- pbEqBelt(nside, k, eps, start_phi, end_phi, start_theta)
    plotPixel( data.frame(theta = S[,1], phi = S[,2]),
               col = col, lwd = lwd, ...)
  }

  if ( !missing(ordering) )
  {
    if ( !identical(ordering, "ring") && !identical(ordering, "nested") )
    {
      stop("ordering, if specified, must be 'ring' or 'nested'")
    }

    centers <- rcosmo::CMBDataFrame(nside = nside, ordering = ordering,
                                    coords = "cartesian")[incl.labels,]

    plot.CMBDataFrame(centers, add = TRUE, col = nums.col,
                              cex = nums.size, labels = incl.labels,
                              font = font)
  }
}


## HELPER FUNCTION 1
plotPixel <- function(S, col = "black", lwd = 1, ...)
{
  C <- sph2car(S)
  rgl::plot3d(C[,"x"],C[,"y"],C[,"z"],
              type = "l", add = TRUE, col = col, lwd = lwd, ...)
}

pbEqBelt<-function(n, k, eps, start_phi,
                   end_phi, start_theta){

  delta_phi <- eps
  phi <- start_phi-delta_phi
  delta_theta <- eps
  theta <- start_theta-delta_theta
  P <- numeric(0)
  index <- 0
  while (phi < end_phi){
    tmp <- suppressWarnings(pbNextEqBelt(n, k, phi, theta,
                            delta_phi, delta_theta))
    index <- index+1
    phi <- tmp[[1]]
    theta <- tmp[[2]]
    dTheta_dPhi<-tmp[[3]]
    a_temp <- cbind(theta, phi)
    P <- rbind(P,a_temp)
    dL_dPhi <- sqrt(dTheta_dPhi^2+sin(theta)^2)
    delta_phi <- eps/dL_dPhi
    delta_theta <- dTheta_dPhi*delta_phi
    # Add the end edge point
    if (phi+delta_phi > end_phi) delta_phi <- end_phi-phi
  }
  return(P)
}


## HELPER FUNCTION 2
pbNextEqBelt <- function(n, k, prev_phi, prev_theta,
                         delta_phi, delta_theta){

  phi <- prev_phi+delta_phi
  a <- 2/3-4*k/(3*n)
  b <- 8/(3*pi)

  #With consideration of the complex values
  if ( (a+b*phi) >= -1 && (a+b*phi) <= 1 )
  {
    theta0 <- acos(a+b*phi)
    if ( (a-b*phi) >= -1 && (a-b*phi) <= 1 )
    {
      theta1 <- acos(a-b*phi)
    }
    else
    {
      theta1 <- Im(acos(a-b*phi))
    }
  }
  else
  {
    theta0 <- Im(acos( a + b*phi ))
    if ( (a - b*phi) >= -1 && (a - b*phi) <= 1 ){
      theta1 <- acos(a-b*phi)
    }else{
      theta1 <- Im(acos( a - b*phi ))
    }
  }

  exp_theta <- prev_theta + delta_theta
  if (abs(theta0-exp_theta) < abs( theta1 - exp_theta )) {
    theta <- Re(theta0)
    dTheta_dPhi <- -b/sin(theta)
  } else {
    theta <- Re(theta1)
    dTheta_dPhi <- b/sin(theta)
  }
  return( list(phi, theta, dTheta_dPhi) )
}


## HELPER FUNCTION 3
pbPolarCap <- function(n, k, eps,
                       start_phi, end_phi, flag){

  delta_phi <- eps
  phi <- start_phi - delta_phi
  P <- numeric(0)

  while (phi < end_phi){
    tmp <- pbNextPolarCap(n, k, phi, delta_phi, flag)
    phi <- tmp[[1]]
    theta <- tmp[[2]]
    dTheta_dPhi <- tmp[[3]]
    a_temp <- cbind(theta, phi)
    P <- rbind(P,a_temp)
    dL_dPhi <- sqrt( dTheta_dPhi^2 + sin(theta)^2 )
    delta_phi<- eps/dL_dPhi

    if (phi+delta_phi > end_phi) delta_phi <- end_phi-phi   #add the end edge point
  }
  return(P)
}


## HELPER FUNCTION 4
pbNextPolarCap <- function(n, k, prev_phi, delta_phi, flag){

  phi <- prev_phi + delta_phi
  b <- -k^2*pi^2/(12*n^2)

  if (flag==0){
    theta <- acos( 1 + b/phi^2 )           #Eq. 19, here cos(theta)=z;
    dTheta_dPhi <- 2*b/(phi^3*sin(theta))
  }else{
    theta <- acos( 1 + b/( phi-pi/2)^2 )   #Eq. 20, here cos(theta)=z;
    dTheta_dPhi <-  2*b/( (phi-pi/2)^3*sin(theta) )
  }

  return( list(phi, theta, dTheta_dPhi) )
}






#' Display the pixels and grandchildren
#'
#' Display the pixels spix at resolution j by colouring
#' in the grandchildren of spix at resolution plot.j
#'
#' @param j The resolution that spix are specified at.
#' @param boundary.j The resolution to display boundaries at. If
#' this is missing then boundaries will not be plotted.
#' @param plot.j The resolution to plot grandchildren at
#' @param spix Integer vector. The pixel indices to display.
#' These must be in nested order.
#' @param incl.labels Integer vector of pixel indices to label at
#' resolution j.
#' @param boundary.col The boundary colour.
#' @param boundary.lwd The boundary line width.
#' @param col The colour to make the grandchildren.
#' @param size The size to make the grandchildren.
#'
#'
#'@examples
#'
#'## Example 1
#'
#' ## Plot base pixels 1,2,3 by colouring their grandchildren at resolution
#' ## 5 (by default). No pixel boundaries.
#' displayPixels(j=0, spix=c(1,2,3))
#'
#' ## Plot base pixels 1,2,3 display and their boundaries (boundary.j=0)
#' displayPixels(0,0, spix=c(1,2,3))
#'
#' ## Plot base pixels 1,2,3 by colouring their grandchildren at resolution 2
#' displayPixels(0,0, plot.j = 2, spix=c(1,2,3))
#'
#'## Example 2
#'
#' demoNeighbours <- function(p,j) {
#'   neighbours(p, j)
#'   displayPixels(boundary.j = j, j = j, plot.j = 5,
#'                 spix = neighbours(p, j),
#'                 boundary.col = "gray",
#'                 boundary.lwd = 1,
#'                 incl.labels = neighbours(p, j),
#'                 col = "blue",
#'                 size = 3)
#'   rcosmo::displayPixelBoundaries(nside = 1, col = "blue", lwd = 3)
#' }
#'
#'demoNeighbours(1,2)
#'
#'@export
displayPixels <- function(boundary.j, j, plot.j = 5, spix,
                          boundary.col = "gray",
                          boundary.lwd = 1,
                          incl.labels = 1:(12*4^boundary.j),
                          col = "blue",
                          size = 3)
{
  if ( !missing(boundary.j) ) {
    rcosmo::displayPixelBoundaries(nside = 2^boundary.j,
                                    ordering = "nested",
                                    nums.col = "red",
                                    col = boundary.col,
                                    lwd = boundary.lwd,
                                    incl.labels = incl.labels)
  }

  # We do this by plotting grandchildren of the siblings
  gchild <- rcosmo::pixelWindow(j1 = j,
                                 j2 = plot.j,
                                 pix.j1 = spix)

  hp <- rcosmo::HPDataFrame(nside = 2^plot.j,
                             spix = gchild)
  plot.HPDataFrame(hp, add = TRUE, col = col, size = size)
}
