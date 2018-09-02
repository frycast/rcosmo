#' Compute spherical harmonic values at given points on the sphere.
#'
#' The function \code{sphericalharmonics} computes the spherical harmonic values for the given 3D Cartesian coordinates.
#' @param L  The degree of spherical harmonic (L=0,1,2,...)
#' @param m  The order number of the degree-L spherical harmonic (m=-L,-L+1,...,L-1,L)
#' @param xyz Dataframe for given points in 3D cartesian coordinates
#' @return values of spherical harmonics
#' @examples
#' ## Calculate spherical harmonic value at the point (0,1,0) with L=5, m=2
#' point<-data.frame(x=0,y=1,z=0)
#' sphericalharmonics(5,2,point)
#' ## Calculate spherical harmonic values at the point (1,0,0), (0,1,0), (0,0,1) with L=5, m=2
#' points<-data.frame(diag(3))
#' sphericalharmonics(5,2,points)
#'
#' @keywords spherical harmonic
#' @references  See also https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
#' It uses equation (7) in Hesse, K., Sloan, I. H., & Womersley, R. S. (2010).
#' Numerical integration on the sphere. In Handbook of Geomathematics (pp. 1185-1219).
#' Springer Berlin Heidelberg.
#' but instead of the order k=1,...,2L+1 in this book we use m=k-L-1.
#'
#' @export
sphericalharmonics <- function(L,m,xyz){
  if (L==0){
      Y <- matrix(1,dim(xyz)[1],1)
      return (Y)
  }

  if ((m<(-L))|| (m>L)){
    sprintf("Warning: m should be in [-%i,%i]", L)
  }else{
    m_abs<-abs(m)
    # normalization constant c_{L,m}
    n <- 1:m_abs
    c_ellm_2 <- prod(1+m_abs/(L-n+1))
    c_ellm <- (sqrt(2)/2^m_abs)*sqrt((2*L+1)/(4*pi)*c_ellm_2)
    # convert xyz to polar coordinates
    x <- xyz[,1]
    y <- xyz[,2]
    z <- xyz[,3]
    phi <- matrix(0,length(x),1)
    t <- phi
    # x3==1
    logic_x3<-(z==1|z==-1)
    phi[logic_x3]<-0
    # x2>=0 & -1<x3<1
    logic_x2_1 <- (y>=0&z<1&z>-1)
    t[logic_x2_1] <- x[logic_x2_1]/sqrt(1-(z[logic_x2_1])^2)
    t[t>1] <- 1
    t[t<(-1)] <- -1
    phi[logic_x2_1] = acos(t[logic_x2_1])
    # x2<0 & x3<1
    logic_x2_2 <- (y<0&z<1&z>-1)
    t[logic_x2_2] <- x[logic_x2_2]/sqrt(1-(z[logic_x2_2])^2)
    t[t>1] <- 1
    t[t<(-1)] <- -1
    phi[logic_x2_2] <- 2*pi-acos(t[logic_x2_2])
    if (m>0){
    # Compute Y_{L,k}
    Y <- c_ellm*(1-z^2)^(m_abs/2)*jacobiPol(m_abs,m_abs,L-m_abs,z)*cos(m_abs*phi)
    # The m-th derivative of the Legendre polynomial
    # P_L^(m)(z) = P_{L-m}^(m,m)(z)
    }else if (m==0){
      Y <- sqrt((2*L+1)/(4*pi))*jacobiPol(0,0,L,z)
    }else{
    Y <- c_ellm*(1-z^2)^(m_abs/2)*jacobiPol(m_abs,m_abs,L-m_abs,z)*sin(m_abs*phi)
    }
  }
  return (Y)
}
