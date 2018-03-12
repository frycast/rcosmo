//THIS SHOULD BE IMPROVED BY NOT NEEDING TO GENERATE ALL COORDINATES WHEN ONLY A HANDFUL ARE REQUIRED
// I.E. WHEN SPIX IS SPECIFIED.


//Includes/namespaces
#include <Rcpp.h>
#include <bitset>
using namespace Rcpp;

//'@title
//'pix2coords
//'@description
//'Converts HEALPix pixel scheme to spherical or
//'Cartesian coordinates.
//'
//'@param Nside The number of cuts to a HEALPix base resolution pixel.
//'@param Nest Set to TRUE for NESTED ordering scheme and FALSE for RING.
//'@param spix Optional integer or vector of sample pixel indices.
//'@param cartesian Set to FALSE to output spherical coordinates
//'or else TRUE for cartesian.
//'
//'@details
//'This is a place holder
//'
//'@return A matrix with columns theta and phi (in that order), or
//' x, y, z (if cartesian = TRUE). Theta (in [0,pi]) is the colatitude
//' in radians measured from the North Pole and phi (in [0, 2*pi])
//' is the longitude in radians measured Eastward. The remaining 3 columns
//' returned are i, j, and p which represent the HEALPix ring index,
//' pixel-in-ring index, and pixel index respectively.
//'
//'@name pix2coords


//Disused cumulative sum helper function.
//// [[Rcpp::export]]
//NumericVector cumsumC(NumericVector x) {
//  int n = x.size();
//  NumericVector out(n);

//  out[0] = x[0];
//  for(int i = 1; i < n; ++i) {
//    out[i] = out[i - 1] + x[i];
//  }
//  return out;
//}


std::string DecToBin(int number)
{
  if ( number == 0 ) return "0";
  if ( number == 1 ) return "1";

  if ( number % 2 == 0 ) {
    return DecToBin(number / 2) + "0";
  } else {
    return DecToBin(number / 2) + "1";
  }
}


int BinToDec(std::string number)
{
  int result = 0, pow = 1;
  for ( int i = number.length() - 1; i >= 0; --i, pow <<= 1 )
    result += (number[i] - '0') * pow;

  return result;
}



//' @export
// [[Rcpp::export]]
NumericMatrix pix2coords(int Nside = 0,
                       bool Nest = true,
                       Rcpp::Nullable<Rcpp::IntegerVector> spix = R_NilValue,
                       bool cartesian = false){
  int Npix = 12*Nside*Nside;
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
    std::sort(sp.begin(),sp.end());
    N = sp.length();
    pvec = sp - 1;

    // Check that spix is valid
    if (N > Npix) {
      throw std::invalid_argument("too many sample pixels, Nside is too small");
    }
    for (int spi = 0; spi < N; ++spi) {
      if (sp[spi] > Npix || sp[spi] <= 0) {
        throw std::invalid_argument("sample pixel is out of range");
      }
    }
    for (int spi = 0; spi < N-1; ++spi) {
      if (sp[spi] == sp[spi+1]) {
        throw std::invalid_argument("duplicate pixel indices not allowed");
      }
    }
  }

  NumericMatrix ang(N,5);
  NumericMatrix car(N,6);
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
  int bpiRingNP = 2*(Nside-1)*Nside - 1;
  int bpiRingNE = (Nside + 1)*4*Nside + bpiRingNP;
  int bpiRingSE = Nside*Nside*4 + bpiRingNE;

  if (Nest == false){
    for (int k = 0; k < N; ++k){
      p = pvec[k];

      if (p <= bpiRingNP){ // North Polar pixel

        double ph = (p+1)/2.0;
        i = floor(sqrt(ph-sqrt(floor(ph)))) + 1;
        j = p + 1 - 2*i*(i-1);
        z = 1 - i*i/(3.0*Nside*Nside);
        phi = PI/(2*i)*(j - 0.5);

      }else if(p <= bpiRingSE){ // Equatorial pixel

        int pp = p - 2*Nside*(Nside - 1);
        i = floor(pp/(4.0*Nside)) + Nside;
        j =  pp % (4*Nside) + 1;
        int s = (i - Nside + 1) % 2;
        z =  4.0/3 - 2.0*i/(3*Nside);
        phi = PI/(2*Nside)*(j - s/2.0);

      }else{ // South Polar pixel

        int ps = Npix - p - 1;
        double ph = (ps+1)/2.0;
        i = floor(sqrt(ph-sqrt(floor(ph)))) + 1;
        j = 4*i - (ps - 2*i*(i-1));
        z = i*i/(3.0*Nside*Nside) - 1;
        phi = PI/(2*i)*(j - 0.5);

      }

      ang(k,0) = acos(z);
      ang(k,1) = phi;
      ang(k,2) = i;
      ang(k,3) = j;
      ang(k,4) = p+1;
    }
  } else {
  // Then NEST = TRUE
    int p = 0;
    int k = 0;

    // Iterate through base resolution pixels
    for (int f = 0; f <= 11; ++f){

        // Iterate through pixels nested in f
        for (int pp = 0; pp < Nside*Nside; ++pp) {
          double F1 = floor(f/4.0) + 2;
          double F2 = 2*(f % 4) - ((int)floor(f/4.0) % 2) + 1;

          // Convert pp to binary and thus get its coordinates x, y.
          std::string bin = DecToBin(pp);
          std::string xbin = "";
          std::string ybin = "";
          for (int l = bin.length() - 1; l >= 1; l -= 2) {
            xbin.insert(0,1,bin[l]);
            ybin.insert(0,1,bin[l-1]);
          }
          if (bin.length() % 2 == 1){
            xbin.insert(0,1,bin[0]);
          }
          int x = BinToDec(xbin);
          int y = BinToDec(ybin);

          int v = x + y;
          int h = x - y;
          i = F1*Nside - v - 1;
          int s = (i - Nside + 1) % 2;

          if (f < 4 && v > Nside - 1) { //North Polar pixel

              j = (F2*i + h + s)/2;
              z = 1 - i*i/(3.0*Nside*Nside);
              phi = PI/(2*i)*(j - s/2.0);

          } else if (f > 7 && v < Nside - 1) { //South Polar pixel

              int is = 4*Nside - i;
              int ss = (is - Nside + 1) % 2;
              j = (F2*is + h + ss)/2;
              z = is*is/(3.0*Nside*Nside) - 1;
              phi = PI/(2*is)*(j - ss/2.0);

          } else { //Equatorial pixel

              j = (F2*Nside + h + s)/2;
              z =  4.0/3 - 2.0*i/(3*Nside);
              phi = PI/(2*Nside)*(j - s/2.0);

          }

          // In case spix is not null we only add the sample pixels
          if (p == pvec[k]) {
            ang(k,0) = acos(z);
            ang(k,1) = phi;
            ang(k,2) = i;
            ang(k,3) = j;
            ang(k,4) = p+1;
            k += 1;
          }

          p += 1;
        }
    }

  }

  // Convert everything to cartesian if required
  if (cartesian == true)
  {
    for (int k = 0; k < N; k++)
    {
      double tempTheta = ang(k,0);
      double tempPhi = ang(k,1);
      int tempi = ang(k,2);
      int tempj = ang(k,3);
      int tempp = ang(k,4);
      car(k,0) = cos(tempPhi)*sin(tempTheta);
      car(k,1) = sin(tempPhi)*sin(tempTheta);
      car(k,2) = cos(tempTheta);
      car(k,3) = tempi;
      car(k,4) = tempj;
      car(k,5) = tempp;
    }
  }

  return *out;
}
