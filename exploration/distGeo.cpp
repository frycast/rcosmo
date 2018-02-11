// WARNING: THIS FUNCTION IS INACCURATE AND USING
// acos( dot(a,b)/(|a||b|) ) WORKS BETTER

//Includes/namespaces
#include <Rcpp.h>
#include <bitset>
using namespace Rcpp;

//'@title
//'distGeo
//'@description
//'Calculates the geodesic distance between points on a sphere
//'using Vincenty's algorithm.
//'
//'@param p1 matrix of longitude and lattitude coordinates (in that order).
//'@param p2 as above. Number of rows of p2 must equal number of rows of p1
//'unless p1 has only one row.
//'
//'@details
//'The radius of the sphere is assumed to be 1.
//'
//'@return
//' If p1 has just one row then the distance from that row to each row
//' in p2 is returned. If p1 has the same number of rows as p2 then the
//' distance from each row of p1 to the corresponding row of p2 is returned.
//'
//'@name distGeo
//'//' @export
// [[Rcpp::export]]
NumericVector distGeo(NumericMatrix p1, NumericMatrix p2) {

  NumericVector lat2 = p2(_,0);
  NumericVector lon2 = p2(_,1);
  int N = lat2.length();
  NumericVector lat1 = rep_len(p1(_,0),N);
  NumericVector lon1 = rep_len(p1(_,1),N);
  NumericVector out(lat2.length());

  if (p1(_,0).size() != 1 && p1(_,0).size() != N) {
    throw std::invalid_argument("p1 has invalid number of rows");
  }

  NumericVector x = sqrt(
      pow(cos(lat2) * sin(lon1 - lon2),2) +
      pow(cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2),2)
      );
  NumericVector y = sin(lat1) * sin(lat2) +
      cos(lat1) * cos(lat2) * cos(lon1 - lon2);

  for (int i = 0; i < lat2.length(); ++i) {
   out(i) = std::atan2(x(i),y(i));
  }

  return(out);
}
