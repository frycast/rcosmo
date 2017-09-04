#' Calculate Jacobi polynomial values of degree L at given point T in [-1,1]. 
#' 
#' @param (a,b) The parameters of Jacobi polynomial  
#' @param L  The degree of Jacobi polynomial
#' @param T Given point in [-1,1]. 
#' @return Jacobi polynomial values 
#' @examples 
#' JacobiRecursive(0,0,5,0)
#' JacobiRecursive(1,2,4,0.5)
#' @keywords Jacobi,Orthogonal polynomials.
#' @source \url{http://dlmf.nist.gov/18.9}
#' @export
JacobiRecursive<-function (a,b,L,T){
  
if (L==0){
  YJ<-matrix(1,length(T),1)
}else if (L==1){
  YJ<-(a-b)/2+(a+b+2)/2*T
}else{
  pMisb1<-matrix(1,length(T),1)
  pMi<-(a-b)/2+(a+b+2)/2*T
  for (i in seq(2,L,1)){
      c<-2*i+a+b
      tmppMisb1<- pMi
      pMi<-((c-1)*c*(c-2)*T*pMi+(c-1)*(a^2-b^2)*pMi-2*(i+a-1)*(i+b-1)*c*pMisb1)/(2*i*(i+a+b)*(c-2))
      pMisb1<-tmppMisb1
  }
  YJ<- pMi
}
  return (YJ)
}
