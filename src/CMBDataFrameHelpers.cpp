#include <Rcpp.h>
using namespace Rcpp;


//'@title
//'pointInConvexPolygon
//'
//'@param df a data.frame with columns x, y, z for cartesian coordinates.
//'The rows represent points on the surface of a unit sphere
//'@param win a data.frame with columns x, y, z for cartesian coordinates.
//'The rows represent clockwise oriented vertices of a convex spherical
//'polygon that lies entirely within one open hemisphere of the unit sphere.
//'
//'@return a logical vector indicated which rows of \code{df}
//'lie within the spherical convex polygon determined by \code{win}
//'
//'@name pointInConvexPolygon
//'
//'@export
// [[Rcpp::export]]
LogicalVector pointInConvexPolygon(DataFrame df, DataFrame win)
{
  NumericVector x = df["x"];
  NumericVector y = df["y"];
  NumericVector z = df["z"];
  NumericVector Vx = win["x"];
  NumericVector Vy = win["y"];
  NumericVector Vz = win["z"];

  int n = df.nrow();
  int k = win.nrow();

  LogicalVector keep(n);
  keep = !keep;
  for ( int i = 0; i < n; i++ )
  {
    for ( int one = 0; one < k; one++ )
    {
      // v2 = v1 + 1 (cyclic)
      int two = (one + 1) % k;

      double det = x[i]*(Vy[one]*Vz[two] - Vz[one]*Vy[two])
                 - y[i]*(Vx[one]*Vz[two] - Vz[one]*Vx[two])
                 + z[i]*(Vx[one]*Vy[two] - Vy[one]*Vx[two]);


      if ( det < - 1e-14 )
      {
        keep(i) = FALSE;
        break;
      }

    }
  }

  return(keep);
}






//'@title
//'pointInDisc
//'
//'@param df a data.frame with columns x, y, z for cartesian coordinates.
//'The rows represent points on the surface of a unit sphere
//'@param win a data.frame with columns x, y, z for the cartesian coordinates
//'of a point on the unit sphere, representing a disc center, and column r for
//'the radius or that disc.
//'
//'@return a logical vector indicated which rows of \code{df}
//'lie within the spherical disc determined by \code{win}
//'
//'@name pointInDisc
//'
//'@export
// [[Rcpp::export]]
LogicalVector pointInDisc(DataFrame df, DataFrame win)
{
  NumericVector x = df["x"];
  NumericVector y = df["y"];
  NumericVector z = df["z"];
  double Vx = win["x"];
  double Vy = win["y"];
  double Vz = win["z"];
  double Vr = win["r"];

  int n = df.nrow();

  LogicalVector keep(n);
  for ( int i = 0; i < n; i++ )
  {
    if ( acos( Vx*x[i] + Vy*y[i] + Vz*z[i] ) <= Vr )
    {
      keep[i] = true;
    }
  }

  return(keep);
}










//'@title
//'maxDist_internal
//'
//'@param cmbdf a \code{data.frame} or \code{\link{CMBDataFrame}}
//'
//'@return the maximum distance between any of the
//'points in \code{cmbdf}
//'
//'@name maxDist_internal
//'
//'@export
// [[Rcpp::export]]
double maxDist_internal(Rcpp::DataFrame cmbdf) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  NumericVector I = cmbdf["I"];

  // Find the maximum distance between any pair of points
  double maxdist = 0;
  for ( int i = 0; i < n-1; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      // Find d(xi, xj)
      double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
      double dist = acos( std::max<double>( dot, -1) );

      if ( dist > maxdist ) maxdist = dist;
    }
  }

  return maxdist;
}







//'@title
//'minDist
//'
//'@param cmbdf a \code{data.frame} or \code{\link{CMBDataFrame}}
//'@param point a point on the unit sphere in cartesian coordinates
//'
//'@return the shortest distance from \code{point} to \code{cmbdf}
//'
//'@name minDist
//'
//'@export
// [[Rcpp::export]]
double minDist(Rcpp::DataFrame cmbdf, NumericVector point) {

  int n = cmbdf.nrow();
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  NumericVector I = cmbdf["I"];
  double px = point[1];
  double py = point[2];
  double pz = point[3];

  // Find the maximum distance between any pair of points
  double mindist = M_PI;
  for ( int i = 0; i < n; i++ )
  {
    // Find d(xi, xj)
    double dot = x[i]*px + y[i]*py + z[i]*pz;
    double dist = acos( std::max<double>( dot, -1) );

    if ( dist < mindist ) mindist = dist;
  }

  return mindist;
}



