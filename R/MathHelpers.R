
## HELPER FUNCTION, CROSS PRODUCT
vector_cross <- function(a, b) {
  if(length(a)!=3 || length(b)!=3){
    stop("Cross product is only defined for 3D vectors.");
  }
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (a[i1]*b[i2] - a[i2]*b[i1])
}







# HELPER FUNCTION: RODRIGUES ROTATION FORMULA SEE WIKIPEDIA PAGE:
# Rotation axis k is defined by being away from a and towards b.
# This function rotates a directly to b (rotating whole sphere points p_xyz).
rodrigues <- function(a,b,p.xyz)
{
  norm.a <- as.numeric(sqrt(sum(a * a)))
  norm.b <- as.numeric(sqrt(sum(b * b)))
  if ( !isTRUE( all.equal(a/norm.a, b/norm.b) ) )
  {
    k <- vector_cross(a,b)
    k <- k/sqrt(sum(k^2)) # normalised k
    K <- matrix( c(  0   , -k[3],  k[2],
                     k[3], 0    , -k[1],
                    -k[2], k[1] ,  0  ), nrow = 3, byrow = TRUE)
    theta <- acos( sum(a*b) / ( norm.a*norm.b ) ) # The angle between a and b
    I <- diag(c(1,1,1))
    R <- I + sin(theta)*K + (1-cos(theta))*K%*%K #Rodrigues' Formula.
    p.xyz <- t(R%*%t(p.xyz))
  }

  p.xyz
}





#' Calculate Jacobi polynomial values
#'
#' Calculate Jacobi polynomial values of degree L at given point T in [-1,1].
#'
#' @param a,b The parameters of Jacobi polynomial
#' @param L  The degree of Jacobi polynomial
#' @param t Given point in [-1,1].
#' @return Jacobi polynomial values
#' @examples
#'
#' jacobiPol(0,0,5,0)
#' jacobiPol(2,-5,2,-1)
#' jacobiPol(1,2,4,0.5)
#'
#' @keywords internal
#' @source \url{http://dlmf.nist.gov/18.9}
#' @export
jacobiPol <- function(a,b,L,t) {

  if (L == 0){
    YJ <- matrix(1,length(t),1)
  }else if (L == 1){
    YJ <- (a - b) / 2 + (a + b + 2) / 2 * t
  }else{
    pMisb1 <- matrix( 1, length(t), 1)
    pMi <- (a-b) / 2 + (a + b + 2) / 2 * t
    for (i in seq(2,L,1)){
      c <- 2 * i + a + b
      tmppMisb1 <- pMi
      pMi <- ((c - 1) * c * (c - 2) * t * pMi + (c - 1) *
                (a ^ 2 - b ^ 2) * pMi - 2 * (i + a - 1) *
                (i + b - 1) * c * pMisb1) /
                (2 * i * (i + a + b) * (c-2) )
      pMisb1 <- tmppMisb1
    }
    YJ <- pMi
  }
  return (YJ)
}






#' Compute spherical harmonic values at given points on the sphere.
#'
#' The function \code{sphericalHarmonics} computes
#' the spherical harmonic values
#' for the given 3D Cartesian coordinates.
#' @param L  The degree of spherical harmonic (L=0,1,2,...)
#' @param m  The order number of the degree-L spherical
#' harmonic (m=-L,-L+1,...,L-1,L)
#' @param xyz Dataframe for given points in 3D cartesian coordinates
#' @return values of spherical harmonics
#'
#' @examples
#' ## Calculate spherical harmonic value at
#' ## the point (0,1,0) with L=5, m=2
#' point<-data.frame(x=0,y=1,z=0)
#' sphericalHarmonics(5,2,point)
#'
#' ## Calculate spherical harmonic values at
#' ## the point (1,0,0), (0,1,0), (0,0,1) with L=5, m=2
#' points<-data.frame(diag(3))
#' sphericalHarmonics(5,2,points)
#'
#' @references  See
#' https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
#'
#' It uses equation (7) in Hesse, K., Sloan, I. H., &
#' Womersley, R. S. (2010).
#' Numerical integration on the sphere.
#' In Handbook of Geomathematics (pp. 1185-1219).
#' Springer Berlin Heidelberg,
#'
#' but instead of the order k=1,...,2L+1 in the book we use m=k-L-1.
#'
#' @export
sphericalHarmonics <- function(L,m,xyz){
  if (L == 0){
    Y <- matrix(1,dim(xyz)[1],1)/sqrt(4*pi)
    return (Y)
  }

  if ((m < -L) || (m > L)){
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
    logic_x2_2 <- (y < 0 & z < 1 & z > -1)
    t[logic_x2_2] <- x[logic_x2_2]/sqrt(1-(z[logic_x2_2])^2)
    t[t>1] <- 1
    t[t<(-1)] <- -1
    phi[logic_x2_2] <- 2*pi-acos(t[logic_x2_2])
    if (m>0){
      # Compute Y_{L,k}
      Y <- c_ellm*(1-z^2)^(m_abs/2)*
        jacobiPol(m_abs,m_abs,L-m_abs,z)*cos(m_abs*phi)
      # The m-th derivative of the Legendre polynomial
      # P_L^(m)(z) = P_{L-m}^(m,m)(z)
    }else if (m==0){
      Y <- sqrt((2*L+1)/(4*pi))*jacobiPol(0,0,L,z)
    }else{
      Y <- c_ellm*(1-z^2)^(m_abs/2)*
        jacobiPol(m_abs,m_abs,L-m_abs,z)*sin(m_abs*phi)
    }
  }
  return (Y)
}
