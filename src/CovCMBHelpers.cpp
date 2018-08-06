#include <Rcpp.h>
using namespace Rcpp;





// // THIS SHOULD BECOME UNUSED ONCE covCMB_internal2 IS COMPLETE
// // IT REMAINS HERE TO DO SOME BENCHMARKING
// //
// //'@title covCMB_internal1
// //'
// //'
// //'see \code{\link{CMBDataFrame}}
// //'
// //'@name covCMB_internal1
// //'
// [[Rcpp::export]]
NumericVector covCMB_internal1(Rcpp::DataFrame cmbdf, NumericVector breaks) {

  int n = cmbdf.nrow();
  int nbreaks = breaks.length();
  int nbin = nbreaks;
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  NumericVector I = cmbdf["I"];
  // We use nbin + 1 to create space for the variance atom (zero bin)
  NumericVector C = NumericVector( nbin + 1 );
  NumericVector B = NumericVector( nbin + 1 );
  NumericMatrix out = NumericMatrix( nbin + 1, 2 ); //(nrow, ncolumn)

  // Deal with the variance atom case C(0) first (i = j)
  // These values also go into the first bin, C(1).
  for ( int i = 0; i < n; i++ )
  {
    double product = I[i]*I[i];
    C[0] = C[0] + product;
  }
  B[0] = n;

  // Now we deal with the case i != j
  // We only need i, j to cover the upper triangle distance matrix i < j
  for ( int i = 0; i < n-1; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      // Find d(xi, xj)
      double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
      double dist = acos( std::max<double>( dot, -1) );

      // Determine which bin d(xi, xj) falls into
      int bin = 1;
      for ( int b = 0; b < nbreaks; b++ )
      {
        if ( dist > breaks[b] ) bin++;
      }

      // Add Xi*Xj to the correct bin
      C[bin] = C[bin] + I[i]*I[j];

      // Increment the number of products in the bin
      // Note that the last bin is not the same size as it
      // contains all distances greater than max.dist
      B[bin] = B[bin] + 1;
    }

  }

  // Elementwise division, normalisation
  C = C/B;

  out( _ , 0) = C;
  out( _ , 1) = B;

  return out;
}



// //Unlike covCMB_internal1, this script acknowledges that there is no need
// //to transform with acos. The breaks are cos(r_i) where r_i is the radius.
// //The bins will have equal area provided that cos(r_i) - cos(r_{i+1}) is fixed.
// //Alternatively, the bins will have equal annular width if r_{i+1} - r_i is fixed,
// //but cos(r_i) must be passed in, regardless of which
// //metric is used to fix distance beforehand.
// //We must note that: For r in (0,pi), cos is a strictly decreasing function,
// //e.g. cos(0) > cos(max.dist)
// //'@title
// //'covCMB_internal2
// //'
// //'
// //'This function acknowledges that there is no need to transform with acos.
// //'The breaks are cos(r_i) where r_i is the radius.
// //'The bins will have equal area provided that cos(r_i) - cos(r_{i+1}) is fixed.
// //'Alternatively, the bins will have equal annular width if r_{i+1} - r_i is fixed,
// //'but cos(r_i) must be passed in, regardless of which
// //'metric is used to fix distance beforehand.
// //'We must note that: For r in (0,pi), cos is a strictly decreasing function,
// //'e.g. cos(0) > cos(max.dist)
// //'
// //'
// //'
// //'@name covCMB_internal2
// //'
// [[Rcpp::export]]
NumericVector covCMB_internal2(Rcpp::DataFrame cmbdf, NumericVector cos_breaks) {

  int n = cmbdf.nrow();
  int nbreaks = cos_breaks.length();
  // We have nbin bins + the zero bin (index 0) + the throw-away bin (index 1, dist > max.dist)
  int nbin = nbreaks;
  NumericVector x = cmbdf["x"];
  NumericVector y = cmbdf["y"];
  NumericVector z = cmbdf["z"];
  NumericVector I = cmbdf["I"];
  // We use nbin + 2 to create space for the variance atom (zero bin) and the throw away bin
  NumericVector C = NumericVector( nbin + 2 );
  NumericVector B = NumericVector( nbin + 2 );
  NumericMatrix out = NumericMatrix( nbin + 2, 2 ); //(nrow, ncolumn)

  // Deal with the variance atom case C(0) first (i = j)
  // These values also go into the first bin, C(1).
  for ( int i = 0; i < n; i++ )
  {
    C[0] = C[0] + I[i]*I[i];
  }
  B[0] = n;

  // Now we deal with the case i != j
  // We only need i, j to cover the upper triangle distance matrix i < j
  for ( int i = 0; i < n-1; i++ )
  {
    for ( int j = i+1; j < n; j++ )
    {
      // Find cosine of d(xi, xj)
      double cos_dist = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];

      // Determine which bin d(xi, xj) falls into.
      // If breaks consists of r_i rather than cos(r_i), then the
      // output vector (except the 0 bin) will need to be reversed (in R)
      // since cosine is a decreasing function.
      int bin = 1;
      for ( int b = 0; b < nbreaks; b++ )
      {
        // cos(r) > cos_breaks[b] iff r < breaks[b]
        if ( cos_dist > cos_breaks[b] )
        {
          bin++;
        }
        else
        {
          break;
        }
      }


      // Add Xi*Xj to the correct bin
      C[bin] = C[bin] + I[i]*I[j];

      // Increment the number of products in the bin
      // Note that the first bin is not the same size as it
      // contains all distances greater than max.dist
      B[bin] = B[bin] + 1;
    }

  }

  // Elementwise division, normalisation
  C = C/B;

  out( _ , 0) = C;
  out( _ , 1) = B;

  return out;
}
