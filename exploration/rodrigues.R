## HELPER FUNCTION, CROSS PRODUCT
vector_cross <- function(a, b) {
  if(length(a)!=3 || length(b)!=3){
    stop("Cross product is only defined for 3D vectors.");
  }
  i1 <- c(2,3,1)
  i2 <- c(3,1,2)
  return (a[i1]*b[i2] - a[i2]*b[i1])
}

# RODRIGUES FUNCTION SEE WIKIPEDIA PAGE:
# Rotation axis k is defined by being away from a and towards b.
# This function rotates a directly to b (rotating whole sphere points p_xyz).
rodrigues <- function(a,b,p_xyz)
{
  norm_a <- sqrt(sum(a * a))
  norm_b <- sqrt(sum(b * b))
  if ( !isTRUE( all.equal(a/norm_a, b/norm_b, check.attributes = FALSE, use.names = FALSE) ) )
  {
    k <- vector_cross(a,b)
    k <- k/sqrt(sum(k^2)) # normalised k
    K <- matrix( c( 0   , -k[3], k[2],
                  k[3], 0    , -k[1],
                 -k[2], k[1] , 0    ), nrow = 3, byrow = TRUE)
    theta <- acos( sum(a*b) / ( norm_a*norm_b ) ) # The angle between a and b
    I <- diag(c(1,1,1))
    R <- I + sin(theta)*K + (1-cos(theta))*K%*%K #Rodrigues' Formula.
    p_xyz <- t(R%*%t(p_xyz))
  }

  p_xyz
}
