#include <Rcpp.h>
using namespace Rcpp;



// CREATE A HELPER FUNCTION FOR geoDist,
// WHEN A MATRIX IS PASSED TO geoDist the
// NUMBER CRUNCHING CAN BE DONE HERE



// // '@title
// // 'minDist_internal1
// // '
// // '@param cmbdf a \code{data.frame}
// // '@param point a point on the unit sphere in cartesian coordinates
// // '
// // '@return the shortest distance from \code{point} to \code{cmbdf}
// // '
// // '@name minDist_internal1
// // '
// // '@export
// [[Rcpp::export]]
double minDist_internal1(Rcpp::DataFrame cmbdf, NumericVector point) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  double px = point[1];
  double py = point[2];
  double pz = point[3];

  // Find the minimum distance
  double mindot = -1;
  for ( int i = 0; i < n; i++ )
  {
    // Find d(xi, xj)
    double dot = x[i]*px + y[i]*py + z[i]*pz;

    if ( dot > mindot ) mindot = dot;
  }

  return acos(std::min<double>(std::max<double>( mindot, -1),1));
}





// // '@title
// // 'minDist_internal2
// // '
// // '@param cmbdf A \code{data.frame}.
// // '
// // '@return The shortest distance between any pair of
// // 'points in \code{cmbdf}
// // '
// // '@name minDist_internal2
// // '
// // '@export
// [[Rcpp::export]]
double minDist_internal2(Rcpp::DataFrame cmbdf) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];

  double mindot = -1;
  for ( int i = 0; i < n-1; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      // Find d(xi, xj)
      double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];

      if ( dot > mindot ) mindot = dot;
    }
  }

  return acos(std::min<double>(std::max<double>( mindot, -1),1));
}




// //'@title
// //'maxDist_internal2
// //'
// // '@param cmbdf a \code{data.frame}
// // '@param point a point on the unit sphere in cartesian coordinates
// // '
// // '@return the longest geodesic distance from \code{point} to \code{cmbdf}
// // '
// // '@name minDist_internal1
// //'
// //'@name maxDist_internal2
// //'
// [[Rcpp::export]]
double maxDist_internal1(Rcpp::DataFrame cmbdf, NumericVector point) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  double px = point[1];
  double py = point[2];
  double pz = point[3];

  // Find the minimum distance
  double maxdot = 1;
  for ( int i = 0; i < n; i++ )
  {
    // Find d(xi, xj)
    double dot = x[i]*px + y[i]*py + z[i]*pz;

    if ( dot < maxdot ) maxdot = dot;
  }

  return acos(std::min<double>(std::max<double>( maxdot, -1),1));
}




// //'@title
// //'maxDist_internal
// //'
// // '@param cmbdf A \code{data.frame}.
// // '
// // '@return The longest geodesic distance between any pair of
// // 'points in \code{cmbdf}
// // '
// // '@name minDist_internal2
// //'
// //'@name maxDist_internal1
// //'
// [[Rcpp::export]]
double maxDist_internal2(Rcpp::DataFrame cmbdf) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];

  double maxdot = 1;
  for ( int i = 0; i < n-1; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      // Find d(xi, xj)
      double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];

      if ( dot < maxdot ) maxdot = dot;
    }
  }

  return acos(std::min<double>(std::max<double>( maxdot, -1),1));
}






