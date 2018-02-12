#include <Rcpp.h>
using namespace Rcpp;


// Note Some values returned differ from Yu Guang's version (below) by 2pi
//'@title
//'car2sph
//'
//'@param df a data.frame with columns labelled x, y and z
//'
//'@return a data.frame with columns lat and long for latitude and
//'longitude
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
  NumericVector lat(n);
  NumericVector lon(n);

  for ( int i = 0; i < n; i++ )
  {

    lat[i] = M_PI/2.0 - acos(z[i]);
    lon[i] = std::atan2(y[i],x[i]);

  }

  DataFrame sph = DataFrame::create( Named("lat") = lat,
                                     Named("long") = lon);

  return sph;

}

//// YU GUANG'S VERSION
// DataFrame car2sph(DataFrame df) {
//
//   int n = df.nrow();
//   NumericVector x = df["x"];
//   NumericVector y = df["y"];
//   NumericVector z = df["z"];
//   NumericVector lat(n);
//   NumericVector lon(n);
//
//   for ( int i = 0; i < n; i++ )
//   {
//     // Latitude
//     lat[i] = M_PI/2.0 - acos(z[i]);
//
//     // Longitude
//     if ( z[i] == 1 || z[i] == -1 )
//     {
//       lon[i] = 0;
//     }
//     else
//     {
//       double t = std::min<double>(
//         std::max<double>( x[i]/sqrt(1-z[i]*z[i]), -1 ), 1);
//
//       if ( y[i] >= 0 )
//       {
//         lon[i] = acos(t);
//       }
//       else
//       {
//         lon[i] = 2*M_PI - acos(t);
//       }
//
//     }
//
//   }
//
//   DataFrame sph = DataFrame::create( Named("lat") = lat,
//                                      Named("long") = lon);
//
//   return sph;
//
// }






//'@title
//'sph2car
//'
//'@param df a data.frame with columns labelled \code{lat} and \code{long}
//'
//'@return a data.frame with columns x, y, z (cartesian coordinates)
//'
//'@name sph2car
//'
//'@export
// [[Rcpp::export]]
DataFrame sph2car(DataFrame df) {

  int n = df.nrow();
  NumericVector lat = df["lat"];
  NumericVector lon = df["long"];
  NumericVector x(n);
  NumericVector y(n);
  NumericVector z(n);

  for ( int i = 0; i < n; i++ )
  {
    double theta = M_PI/2 - lat[i];
    x[i] = sin(theta)*cos(lon[i]);
    y[i] = sin(theta)*sin(lon[i]);
    z[i] = cos(theta);
  }

  DataFrame xyz = DataFrame::create( Named("x") = x,
                                     Named("y") = y,
                                     Named("z") = z);

  return xyz;

}
