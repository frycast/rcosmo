//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;







// HELPER FUNCTION FOR nest2ring
// mkpix2xyC calculates the vector of x and y in the face from pixel number
// for nested order.
//
// OUTPUTS:
//  pix2$x  - pix index for x
//
//  pix2$y  - pix index for y
// This is a helper function
// [[Rcpp::export]]
IntegerMatrix mkpix2xyC(int nside = 1024) {

  IntegerMatrix pix2xy(nside, 2);

  //-----pix2x <- matrix(rep(nside,0),ncol=nside)
  //-----pix2y <- matrix(rep(nside,0),ncol=nside)
  int jpix;
  int ix;
  int iy;
  int ip;
  int id = 0;

  for ( int kpix = 0;  kpix < nside; kpix++ ) {
    jpix = kpix;
    ix = 0;
    iy = 0;
    ip = 1;
    // bit position in x and y
    while ( !(jpix==0) ) {

      // bit value in k//pix, for ix
      id = jpix % 2;
      jpix = jpix/2;
      ix = id*ip+ix;
      // bit value in kpix, for iy
      id = jpix % 2;
      jpix = jpix/2;
      iy = id*ip+iy;
      // next bit in x and y
      ip = 2*ip;

    }

    // kpix in 0:31
    pix2xy(kpix,0) = ix;
    // kpix in 0:31
    pix2xy(kpix,1) = iy;
  }

  return pix2xy;
}





















//' @title Convert nest to ring ordering
//'
//' @description
//' Convert from "nested" to "ring" ordering
//'
//' \code{nest2ring} computes the HEALPix pixel index
//' in the "ring" ordering scheme from the pixel index
//' in the "nested" ordering scheme.
//'
//' @param nside is the HEALPix nside parameter.
//'
//' @param pix is the set or subset of pixel indices at nside.
//' If pix is left blank then all pixels are converted.
//'
//' @return the output is the corresponding set of pixel in
//' the ring ordering scheme.
//'
//' @examples
//' # compute HEALPix indices in the ring ordering scheme
//' nside <- 8
//' pix <-c(1,2,23)
//' nest2ring(nside,pix)
//'
//' @name nest2ring
//' @export
// [[Rcpp::export]]
IntegerVector nest2ring(int nside, IntegerVector pix ) {

  // number of pix at nside
  int nPix = 12*nside*nside;

  int N = pix.length();

  IntegerVector jrll = IntegerVector::create(
    2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 );
  IntegerVector jpll = IntegerVector::create(
    1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 );

  IntegerMatrix pix2xy = mkpix2xyC(); // ---------- IN THE ORIGINAL nside = 1024 always, why?

  // number of pixels in a face
  int npface = nside*nside;
  int nl4 = 4*nside;

  // find the face number
  // face number in 0:11
  IntegerVector face_num(N);
  for (int k = 0; k < N; k++)
  {
    face_num[k] = (int) ( (pix[k]-1)/npface );
  }
  IntegerVector ipf(N);
  for (int k = 0; k < N; k++)
  {
    ipf[k] = (int) ( ((int)pix[k]-1) % npface );
  }

  // finds the x,y on the face (starting from the lowest corner)
  // from the pixel number
  IntegerVector ix(N);
  IntegerVector iy(N);
  int scalemlv = 1;
  int ismax = 4;

  int n1 = 1024;
  for (int i = 0; i <= ismax; i++) {
    IntegerVector ip_low(N);
    for (int k = 0; k < N; k++)
    {
      ip_low[k] = (int)ipf[k] % n1;
      ix[k] = (int) (ix[k] + scalemlv*pix2xy(ip_low[k],0) );
      iy[k] = (int) (iy[k] + scalemlv*pix2xy(ip_low[k],1) );
      ipf[k] = (int) (ipf[k]/n1);
    }
    scalemlv = scalemlv*32;
  }
  for (int k = 0; k < N; k++)
  {
    ix[k] = (int) (ix[k] + scalemlv*pix2xy(ipf[k],0) );
    iy[k] = (int) (iy[k] + scalemlv*pix2xy(ipf[k],1) );
  }


  // transform to (horizontal, vertical) coordinates
  // 'vertical' in 0:2*(nside-1)
  IntegerVector jrt = ix + iy;
  // 'horizontal' in -nside+1:nside-1
  IntegerVector jpt = ix - iy;

  // Find z coordinate on S^2
  // ring number in 1:4*nside-1
  IntegerVector jr(N);
  for (int k = 0; k < N; k++ )
  {
    jr[k] = (int) (jrll[face_num[k]]*nside - jrt[k] - 1);
  }

  // initialisation
  IntegerVector nr(N);
  IntegerVector kshift(N);
  IntegerVector n_before(N);

  for ( int k = 0; k < N; k++ )
  {
    // north pole area
    if ( jr[k] < nside )
    {
      n_before[k] = (int) ( 2*jr[k]*(jr[k] - 1) );
      kshift[k] = 0;
      nr[k] = jr[k];
    }
    // equatorial area
    else if ( (jr[k] >= nside) && (jr[k] <= 3*nside) )
    {
      n_before[k] = (int) (2*nside*(2*jr[k] - nside - 1));
      kshift[k] = ((int)jr[k] - nside) % 2;
      nr[k] = nside;
    }
    // south pole area
    else if ( jr[k] > 3*nside )
    {
      int nrS = (int)(nl4 - jr[k]);
      n_before[k] = (int)(nPix - 2*nrS*(nrS + 1));
      kshift[k] = 0;
      nr[k] = nrS;
    }
  }


  // computes the phi coordinate on S^2, in [0,2*pi)
  // 'phi' number in the ring in 1:4*nr
  IntegerVector jp(N);
  for (int k = 0; k < N; k++) {
    jp[k] = (int)((jpll[face_num[k]]*nr[k] + jpt[k] + 1 + kshift[k])/2);
    if (jp[k] > nl4) {
      jp[k] = (int) (jp[k] - nl4);
    } else if (jp[k] < 1) {
      jp[k] = (int) (jp[k] + nl4);
    }
  }

  // index in 0:nPix-1
  IntegerVector ipring = (n_before + jp);

  return ipring;
}












// //'@title
// //'pix2coords_internal
// //'@description
// //'Converts HEALPix pixel scheme to spherical or
// //'Cartesian coordinates.
// //'
// //'@param nside The number of cuts to a HEALPix base resolution pixel.
// //'@param nested Set to TRUE for NESTED ordering scheme and FALSE for RING.
// //'@param spix Optional integer or vector of sample pixel indices.
// //'@param cartesian Set to FALSE to output spherical coordinates
// //'or else TRUE for cartesian.
// //'
// //'@details
// //'This is a place holder
// //'
// //'@return A matrix with columns theta and phi (in that order), or
// //' x, y, z (if cartesian = TRUE). Theta (in [0,pi]) is the colatitude
// //' in radians measured from the North Pole and phi (in [0, 2*pi])
// //' is the longitude in radians measured Eastward. The remaining 3 columns
// //' returned are i, j, and p which represent the HEALPix ring index,
// //' pixel-in-ring index, and pixel index respectively.
// //'
// //'@name pix2coords_internal
// //'
// [[Rcpp::export]]
NumericMatrix pix2coords_internal(int nside = 0,
                   bool nested = true,
                   Rcpp::Nullable<Rcpp::IntegerVector> spix = R_NilValue,
                   bool cartesian = false){
  int Npix = 12*nside*nside;
  double fact1 = 1.5*nside;
  double fact2 = 3.0*nside*nside;
  int nl2 = 2*nside;
  int nl4 = 4*nside;
  double z = 0;
  double phi = 0;
  int i = 0;
  int j = 0;
  IntegerVector pvec;
  int p;
  int N;

  if (spix.isNull()) {
    N = Npix;
    pvec = IntegerVector(Npix);
    for (int p = 0; p < Npix; ++p){
      pvec[p] = p;
    }
  } else {
    IntegerVector sp(spix.get());
    // We used to sort spix before HPDataFrame but I
    // am not sure if this is necessary anymore:
    //std::sort(sp.begin(),sp.end());
    N = sp.length();
    pvec = sp - 1;

    for (int spi = 0; spi < N; ++spi) {
      if (sp[spi] > Npix || sp[spi] <= 0) {
        throw std::invalid_argument("sample pixel is out of range");
      }
    }
// We banned duplicate pixel indices before HPDataFrame was created,
// but I'm not sure if this is necessary anymore:
//    for (int spi = 0; spi < N-1; ++spi) {
//      if (sp[spi] == sp[spi+1]) {
//        throw std::invalid_argument("duplicate pixel indices not allowed");
//      }
//    }
  }

  NumericMatrix ang(N,2);
  NumericMatrix car(N,3);
  NumericMatrix * out;

  if (cartesian == false)
  {
    out = &ang;
  }
  else
  {
    out = &car;
  }

  // Regional boundary pixel indices for Ring ordering scheme
  int bpiRingNP = 2*(nside-1)*nside - 1;
  int bpiRingNE = (nside + 1)*4*nside + bpiRingNP;
  int bpiRingSE = nside*nside*4 + bpiRingNE;

  if (nested == true){
    pvec = nest2ring( nside, pvec + 1 ) - 1;
  }

  for (int k = 0; k < N; ++k){
    p = pvec[k];

    if (p <= bpiRingNP){ // North Polar pixel

      double ph = (p+1)/2.0;
      i = floor(sqrt(ph-sqrt(floor(ph)))) + 1;
      j = p + 1 - 2*i*(i-1);
      z = 1 - i*i/fact2;
      phi = (j - 0.5)*PI/(2.0*i);

    } else if (p <= bpiRingSE){

      int pp = p - 2*nside*(nside - 1);
      i = pp/nl4 + nside;
      j =  pp % nl4 + 1;

      double s = 0.5*(1 + ((i + nside) % 2));
      z =  (nl2 - i)/fact1;
      phi = (j - s)*PI/(2.0*nside);

    } else { // South Polar pixel

      int ps = Npix - p;
      double ph = ps/2.0;
      i = floor(sqrt(ph-sqrt(floor(ph)))) + 1;
      j = 4*i + 1 - (ps - 2*i*(i-1));
      z = i*i/fact2 - 1;
      phi = (j - 0.5)*PI/(2*i);

    }

    if ( cartesian == false ) {

      ang(k,0) = acos(z);
      ang(k,1) = phi;
      //ang(k,2) = i;
      //ang(k,3) = j;
      //ang(k,4) = p+1;

    } else {

      double sth = sqrt(1-z)*sqrt(1+z);
      car(k,0) = sth*cos(phi);
      car(k,1) = sth*sin(phi);
      car(k,2) = z;
      //car(k,3) = i;
      //car(k,4) = j;
      //car(k,5) = p+1;

    }
  }

  return *out;
}






// //'@title
// //'car2sph
// //'
// //'@param df a data.frame with columns labelled x, y and z
// //'
// //'@return a data.frame with columns theta and phi for colatitude and
// //'longitude in ranges \eqn{[0,pi]} and \eqn{[0,2pi]} respectively
// //'
// //'@name car2sph
// //'
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

    if ( phi_i >= 2*M_PI-1e-13 )
    {
      phi_i = phi_i - 2*M_PI;
    }

    if ( phi_i < -1e-13 )
    {
      phi[i] = 2*M_PI + phi_i;
    }
    else if ( phi_i < 0 )
    {
      phi[i] = 0;
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






// //'@title
// //'sph2car
// //'
// //'@param df a data.frame with columns labelled \code{theta} and \code{phi}
// //'for colatitude and longitude respectively
// //'
// //'@return a data.frame with columns x, y, z (cartesian coordinates)
// //'
// //'@name sph2car
// //'
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
    if ( theta[i] == 0 )
    {
      x[i] = 0;
      y[i] = 0;
      z[i] = 1;
    }
    else if ( theta[i] == M_PI )
    {
      x[i] = 0;
      y[i] = 0;
      z[i] = -1;
    }
    else
    {
      x[i] = sin(theta[i])*cos(phi[i]);
      y[i] = sin(theta[i])*sin(phi[i]);
      z[i] = cos(theta[i]);
    }
  }

  DataFrame xyz = DataFrame::create( Named("x") = x,
                                     Named("y") = y,
                                     Named("z") = z);

  return xyz;

}






// //'@title
// //'pointInConvexPolygonHP
// //'
// //'@param nside the nside parameter at which to find pixels
// //'@param nested Set to TRUE for NESTED ordering scheme and FALSE for RING
// //'@param win a data.frame with columns x, y, z for cartesian coordinates
// //'The rows represent clockwise oriented vertices of a convex spherical
// //'polygon that lies entirely within one open hemisphere of the unit sphere
// //'@param spix Optional integer or vector of sample pixel indices. If \code{spix}
// //'is unspecified then all pixels at \code{nside} are used
// //'
// //'@return a logical vector indicated which pixels in \code{spix}
// //'lie within the spherical convex polygon determined by \code{win}
// //'
// //'@name pointInConvexPolygonHP
// //'
// [[Rcpp::export]]
LogicalVector pointInConvexPolygonHP(int nside, bool nested, DataFrame win,
                                     Rcpp::Nullable<Rcpp::IntegerVector> spix = R_NilValue)
{
  NumericMatrix crds = pix2coords_internal(nside, nested, spix, true);

  NumericVector x = crds(_,0);
  NumericVector y = crds(_,1);
  NumericVector z = crds(_,2);
  NumericVector Vx = win["x"];
  NumericVector Vy = win["y"];
  NumericVector Vz = win["z"];

  int n = crds.nrow();
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




// //'@title
// //'pointInDiscHP
// //'
// //'@param nside the nside parameter at which to find pixels
// //'@param nested Set to TRUE for NESTED ordering scheme and FALSE for RING
// //'@param win a data.frame with columns x, y, z for the cartesian coordinates
// //'of a point on the unit sphere, representing a disc center, and column r for
// //'the radius or that disc
// //'@param spix Optional integer or vector of sample pixel indices. If \code{spix}
// //'is unspecified then all pixels at \code{nside} are used
// //'
// //'@return a logical vector indicated which pixels in \code{spix}
// //'lie within the spherical disc determined by \code{win}
// //'
// //'@name pointInDiscHP
// //'
// [[Rcpp::export]]
LogicalVector pointInDiscHP(int nside, bool nested, DataFrame win,
                            Rcpp::Nullable<Rcpp::IntegerVector> spix = R_NilValue)
{
  NumericMatrix crds = pix2coords_internal(nside, nested, spix, true);

  NumericVector x = crds(_,0);
  NumericVector y = crds(_,1);
  NumericVector z = crds(_,2);
  double Vx = win["x"];
  double Vy = win["y"];
  double Vz = win["z"];
  double Vr = win["r"];

  int n = crds.nrow();

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









// //'@title
// //'pointInConvexPolygon
// //'
// //'@param df a data.frame with columns x, y, z for cartesian coordinates.
// //'The rows represent points on the surface of a unit sphere
// //'@param win a data.frame with columns x, y, z for cartesian coordinates.
// //'The rows represent clockwise oriented vertices of a convex spherical
// //'polygon that lies entirely within one open hemisphere of the unit sphere.
// //'
// //'@return a logical vector indicated which rows of \code{df}
// //'lie within the spherical convex polygon determined by \code{win}
// //'
// //'@name pointInConvexPolygon
// //'
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






// //'@title
// //'pointInDisc
// //'
// //'@param df a data.frame with columns x, y, z for cartesian coordinates.
// //'The rows represent points on the surface of a unit sphere
// //'@param win a data.frame with columns x, y, z for the cartesian coordinates
// //'of a point on the unit sphere, representing a disc center, and column r for
// //'the radius or that disc.
// //'
// //'@return a logical vector indicated which rows of \code{df}
// //'lie within the spherical disc determined by \code{win}
// //'
// //'@name pointInDisc
// //'
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









