#include <Rcpp.h>
using namespace Rcpp;



//'@title
//'car2sph
//'
//'@param df a data.frame with columns labelled x, y and z
//'
//'@return a data.frame with columns theta and phi for colatitude and
//'longitude in ranges \eqn{[0,pi]} and \eqn{[0,2pi]} respectively
//'
//'@name car2sph
//'
//'@export
// [[Rcpp::export]]
DataFrame car2sph(DataFrame df) {

  int n = df.nrow();
  NumericVector x = df["x"];
  NumericVector y = df["y"];
  NumericVector z = df["z"];
  NumericVector theta(n);
  NumericVector phi(n);

  for ( int i = 0; i < n; i++ )
  {

    theta[i] = acos(z[i]);
    double phi_i = std::atan2(y[i],x[i]);
    if ( phi_i < 0 )
    {
      phi[i] = 2*M_PI + phi_i;
    }
    else
    {
      phi[i] = phi_i;
    }

  }

  DataFrame sph = DataFrame::create( Named("theta") = theta,
                                     Named("phi") = phi);

  return sph;

}

//// YU GUANG'S VERSION
// DataFrame car2sph(DataFrame df) {
//
//   int n = df.nrow();
//   NumericVector x = df["x"];
//   NumericVector y = df["y"];
//   NumericVector z = df["z"];
//   NumericVector theta(n);
//   NumericVector phi(n);
//
//   for ( int i = 0; i < n; i++ )
//   {
//     // colatitude
//     theta[i] = M_PI/2.0 - acos(z[i]);
//
//     // Longitude
//     if ( z[i] == 1 || z[i] == -1 )
//     {
//       phi[i] = 0;
//     }
//     else
//     {
//       double t = std::min<double>(
//         std::max<double>( x[i]/sqrt(1-z[i]*z[i]), -1 ), 1);
//
//       if ( y[i] >= 0 )
//       {
//         phi[i] = acos(t);
//       }
//       else
//       {
//         phi[i] = 2*M_PI - acos(t);
//       }
//
//     }
//
//   }
//
//   DataFrame sph = DataFrame::create( Named("theta") = theta,
//                                      Named("phi") = phi);
//
//   return sph;
//
// }






//'@title
//'sph2car
//'
//'@param df a data.frame with columns labelled \code{theta} and \code{phi}
//'for colatitude and longitude respectively
//'
//'@return a data.frame with columns x, y, z (cartesian coordinates)
//'
//'@name sph2car
//'
//'@export
// [[Rcpp::export]]
DataFrame sph2car(DataFrame df) {

  int n = df.nrow();
  NumericVector theta = df["theta"];
  NumericVector phi = df["phi"];
  NumericVector x(n);
  NumericVector y(n);
  NumericVector z(n);

  for ( int i = 0; i < n; i++ )
  {
    x[i] = sin(theta[i])*cos(phi[i]);
    y[i] = sin(theta[i])*sin(phi[i]);
    z[i] = cos(theta[i]);
  }

  DataFrame xyz = DataFrame::create( Named("x") = x,
                                     Named("y") = y,
                                     Named("z") = z);

  return xyz;

}
