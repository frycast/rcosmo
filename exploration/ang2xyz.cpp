//Includes/namespaces
#include <Rcpp.h>
#include <bitset>
using namespace Rcpp;

//'@title
//'ang2xyzC
//'@description
//'Converts spherical coordinates to cartesian coordinates.
//'
//'@param theta A vector of colattitude.
//'@param phi A vector of longitude.
//'
//'@details
//'This is a place holder
//'@name ang2xyzC
//'
//'@export
// [[Rcpp::export]]
NumericMatrix pix2angC(NumericVector theta, NumericVector phi){
  if (theta.length() != phi.length()) {
    throw std::invalid_argument("theta and phi lengths do not match");
  }
  NumericMatrix xyz(theta.length(), 3);
  
  // Finish this later
  
  return(xyz);
}