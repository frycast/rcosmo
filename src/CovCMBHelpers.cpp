#include <Rcpp.h>
using namespace Rcpp;





// THIS SHOULD BECOME UNUSED ONCE covCMB_internal2 IS COMPLETE
// IT REMAINS HERE TO DO SOME BENCHMARKING
//
//'@title covCMB_internal1
//'
//'
//'see \code{\link{CMBDataFrame}}
//'
//'@name covCMB_internal1
//'
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



//Unlike covCMB_internal1, this script acknowledges that there is no need
//to transform with acos. The breaks are cos(r_i) where r_i is the radius.
//The bins will have equal area provided that cos(r_i) - cos(r_{i+1}) is fixed.
//Alternatively, the bins will have equal annular width if r_{i+1} - r_i is fixed,
//but cos(r_i) must be passed in, regardless of which
//metric is used to fix distance beforehand.
//We must note that: For r in (0,pi), cos is a strictly decreasing function,
//e.g. cos(0) > cos(max.dist)
//'@title
//'covCMB_internal2
//'
//'
//'This function acknowledges that there is no need to transform with acos.
//'The breaks are cos(r_i) where r_i is the radius.
//'The bins will have equal area provided that cos(r_i) - cos(r_{i+1}) is fixed.
//'Alternatively, the bins will have equal annular width if r_{i+1} - r_i is fixed,
//'but cos(r_i) must be passed in, regardless of which
//'metric is used to fix distance beforehand.
//'We must note that: For r in (0,pi), cos is a strictly decreasing function,
//'e.g. cos(0) > cos(max.dist)
//'
//'
//'
//'@name covCMB_internal2
//'
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






// //'@title
// //'geoDistList
// //'@description
// //'Create a list of all geodesic distances between points on the unit
// //'sphere corresponding to the rows of the data.frame cmbdf.
// //'
// //'@param cmbdf a data.frame whose first 3 columns represent the cartesian
// //'coordinates x,y,z of nrow(cmbdf) points on a unit sphere.
// //'
// //'@return
// //'Let \eqn{x_i, x_j} be the points represented by row i and row j of
// //'cmbdf, where i < j, and let L denote the output list
// //'\code{L <- geoDistList(cmbdf)}. Then the distance \eqn{d(x_i,x_j)}
// //'is stored in \code{L[[i]][j-i]}.
// //'
// //'@name geoDistList
// //'
// //'@export
// // [[Rcpp::export]]
// Rcpp::List geoDistList(Rcpp::DataFrame cmbdf) {
//
//   int n = cmbdf.nrow();
//   NumericVector x = cmbdf["x"];
//   NumericVector y = cmbdf["y"];
//   NumericVector z = cmbdf["z"];
//   Rcpp::List L(n-1);
//
//   for ( int i = 0; i < n-1; i++ )
//   {
//     std::vector<double> * dists = new std::vector<double>( n-i-1 );
//
//     for ( int j = i+1; j < n; j++ )
//     {
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       (*dists)[j-i-1] = acos( std::max<double>( dot, -1) );
//     }
//
//     L[i] = *dists;
//
//     delete dists;
//   }
//
//   return L;
//
// }
//
//
//
//
// //'@title
// //'distBinList
// //'@description
// //'Find all geodesic distances between points on the unit
// //'sphere corresponding to the rows of the data.frame cmbdf,
// //'then categorise these distances according to the intervals
// //'given by the breaks argument.
// //'
// //'@param cmbdf a data.frame whose first 3 columns represent the cartesian
// //'coordinates x,y,z of nrow(cmbdf) points on a unit sphere.
// //'@param breaks a vector, sorted from lowest to highest,
// //'specifying the break points for the intervals that are
// //'used to categorise the geodesic distances. The intervals are
// //'open at left and closed at right.
// //'
// //'@return
// //'Let \eqn{x_i, x_j} be the points represented by row i and row j of
// //'cmbdf, where\eqn{i < j}, and let \eqn{L} denote the output list
// //'\code{L <- distBinList(cmbdf, breaks)}. Suppose the distance
// //'\eqn{d(x_i,x_j)} falls into the \eqn{k^{th}} interval determined
// //'by \code{breaks}, then \code{L[[i]][j-i]} will hold the value \eqn{k}.
// //'
// //'@name distBinList
// //'
// //'@export
// // [[Rcpp::export]]
// Rcpp::List distBinList(Rcpp::DataFrame cmbdf, NumericVector breaks) {
//
//   int n = cmbdf.nrow();
//   NumericVector x = cmbdf["x"];
//   NumericVector y = cmbdf["y"];
//   NumericVector z = cmbdf["z"];
//   Rcpp::List L(n-1);
//
//   for ( int i = 0; i < n-1; i++ )
//   {
//     std::vector<int> * binned = new std::vector<int>( n-i-1 );
//
//     for ( int j = i+1; j < n; j++ )
//     {
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       //(*dists)[j-i-1] = acos( std::max<double>( dot, -1) );
//       double dist = acos( std::max<double>( dot, -1) );
//
//       int bin = 0;
//       for ( int b = 0; b < breaks.length(); b++ )
//       {
//         if ( dist > breaks[b] ) bin++;
//       }
//
//       (*binned)[j-i-1] = bin;
//     }
//
//     L[i] = *binned;
//
//     delete binned;
//   }
//
//   return L;
//
// }





//
//
// //'@title
// //'covCMB_internal2
// //'
// //'@name covCMB_internal2
// //'
// //'@export
// // [[Rcpp::export]]
// NumericVector covCMB_internal2(Rcpp::DataFrame cmbdf, int nbin) {
//
//   int nbreaks = nbin - 1;
//   NumericVector breaks( nbreaks );
//   int n = cmbdf.nrow();
//   NumericVector x = cmbdf["x"];
//   NumericVector y = cmbdf["y"];
//   NumericVector z = cmbdf["z"];
//   NumericVector I = cmbdf["I"];
//   // We use nbin + 1 to create space for the variance atom (zero bin)
//   NumericVector C = NumericVector( nbin + 1 );
//   NumericVector B = NumericVector( nbin + 1 );
//
//   // Find the maximum distance between any pair of points
//   double maxdist = 0;
//   for ( int i = 0; i < n-1; i++ )
//   {
//     for ( int j = i+1; j < n; j++ )
//     {
//       // Find d(xi, xj)
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       double dist = acos( std::max<double>( dot, -1) );
//
//       if ( dist > maxdist ) maxdist = dist;
//     }
//   }
//
//   // Create breaks sequence using maxdist and number of bins
//   for (int i = 0;  i < nbreaks; i++ )
//   {
//     breaks[i] = (i+1)*maxdist/nbin;
//   }
//
//   // Deal with the variance case C(0) first
//   for ( int i = 0; i < n; i++ )
//   {
//     C[0] = C[0] + I[i]*I[i];
//   }
//   B[0] = n;
//
//   // We only need i, j to cover the upper triangle distance matrix i < j
//   for ( int i = 0; i < n-1; i++ )
//   {
//     for ( int j = i+1; j < n; j++ )
//     {
//       // Find d(xi, xj)
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       double dist = acos( std::max<double>( dot, -1) );
//
//       // Determine which bin d(xi, xj) falls into
//       int bin = 1;
//       for ( int b = 0; b < nbreaks; b++ )
//       {
//         if ( dist > breaks[b] ) bin++;
//       }
//
//       // Add Xi*Xj to the correct bin
//       C[bin] = C[bin] + I[i]*I[j];
//
//       // Increment the number of products in the bin
//       B[bin] = B[bin] + 1;
//     }
//
//   }
//
//   // Elementwise division, normalisation
//   C = C/B;
//
//   return C;
// }
//
//
//












//
// NumericVector covCMB_internal1(Rcpp::DataFrame cmbdf, NumericVector breaks) {
//
//   int n = cmbdf.nrow();
//   int nbin = breaks.length() - 1;
//   NumericVector x = cmbdf[0];
//   NumericVector y = cmbdf[1];
//   NumericVector z = cmbdf[2];
//   NumericVector I = cmbdf[3];
//   NumericVector C = NumericVector( nbin + 1 );
//   NumericVector B = NumericVector( nbin + 1 );
//
//   for ( int i = 0; i < n-1; i++ )
//   {
//     for ( int j = i+1; j < n; j++ )
//     {
//       // Find d(xi, xj)
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       double dist = acos( std::max<double>( dot, -1) );
//
//       // Determine which bin d(xi, xj) falls into
//       int bin = 0;
//       for ( int b = 0; b <= nbin; b++ )
//       {
//         if ( dist > breaks[b] ) bin++;
//       }
//
//       // Add Xi*Xj to the correct bin
//       C[bin] = C[bin] + I[i]*I[j];
//
//       // Increment the number of products in the bin
//       B[bin] = B[bin] + 1;
//     }
//
//   }
//
//   // Elementwise division, normalisation
//   C = C/B;
//
//   return C;
// }
//
//
// NumericVector covCMB_internal2(Rcpp::DataFrame cmbdf, int nbin) {
//
//   NumericVector breaks( nbin + 1 );
//   int n = cmbdf.nrow();
//   NumericVector x = cmbdf[0];
//   NumericVector y = cmbdf[1];
//   NumericVector z = cmbdf[2];
//   NumericVector I = cmbdf[3];
//   NumericVector C = NumericVector( nbin + 1 );
//   NumericVector B = NumericVector( nbin + 1 );
//
//   // Find the maximum distance between any pair of points
//   double maxdist = 0;
//   for ( int i = 0; i < n-1; i++ )
//   {
//     for ( int j = i+1; j < n; j++ )
//     {
//       // Find d(xi, xj)
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       double dist = acos( std::max<double>( dot, -1) );
//
//       if ( dist > maxdist ) maxdist = dist;
//     }
//   }
//
//   // Create breaks sequence using maxdist and number of bins
//   for (int i = 0; i <= nbin + 1; i++ )
//   {
//     breaks[i] = i*maxdist/nbin;
//   }
//
//   for ( int i = 0; i < n-1; i++ )
//   {
//     for ( int j = i+1; j < n; j++ )
//     {
//       // Find d(xi, xj)
//       double dot = x[i]*x[j] + y[i]*y[j] + z[i]*z[j];
//       double dist = acos( std::max<double>( dot, -1) );
//
//       // Determine which bin d(xi, xj) falls into
//       int bin = 0;
//       for ( int b = 0; b <= nbin; b++ )
//       {
//         if ( dist > breaks[b] ) bin++;
//       }
//
//       // Add Xi*Xj to the correct bin
//       C[bin] = C[bin] + I[i]*I[j];
//
//       // Increment the number of products in the bin
//       B[bin] = B[bin] + 1;
//     }
//
//   }
//
//   // Elementwise division, normalisation
//   C = C/B;
//
//   return C;
// }
