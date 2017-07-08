//Includes/namespaces
#include <Rcpp.h>
#include <bitset>
using namespace Rcpp;

//'@title
//'pix2angC
//'@description
//'Converts HEALPix pixel scheme to spherical coordinates.
//'
//'@param Nside The number of cuts to a HEALPix base resolution pixel.
//'
//'@param Nest A boolean, TRUE if ordering scheme is NESTED and FALSE otherwise.
//'
//'@details
//'This is a place holder
//'@name pix2angC


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
NumericMatrix pix2angC(int Nside, bool Nest) {
  int Npix = 12*Nside*Nside;
  NumericMatrix ang(Npix, 4); // Only set to 4 for debugging to give i, j

  // Regional boundary pixel indices for Ring ordering scheme
  int bpiRingNP = 2*(Nside-1)*Nside - 1;
  int bpiRingNE = (Nside + 1)*4*Nside + bpiRingNP;
  int bpiRingSE = Nside*Nside*4 + bpiRingNE;

  double z = 0;
  double phi = 0;
  int i = 0;
  int j = 0;

  if (Nest == FALSE){
    for (int p = 0; p < Npix; ++p){

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

      ang(p,0) = acos(z);
      ang(p,1) = phi;
      ang(p,2) = i;
      ang(p,3) = j;
    }
  } else {
  // Then NEST = TRUE
    int p = 0;
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
          for (int k = bin.length() - 1; k >= 1; k -= 2) {
            xbin.insert(0,1,bin[k]);
            ybin.insert(0,1,bin[k-1]);
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

        ang(p,0) = acos(z);
        ang(p,1) = phi;
        ang(p,2) = i;
        ang(p,3) = j;
        p += 1;
        }
    }

  }

  return ang;
}
