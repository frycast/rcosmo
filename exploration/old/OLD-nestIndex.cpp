//Includes/namespaces
#include <Rcpp.h>
#include <bitset>
using namespace Rcpp;

//'@title
//'nestIndexC
//'@description
//'Provides the index needed to convert between NEST and RING
//'ordering schemes for a given Nside.
//'
//'@param Nside The number of cuts to a HEALPix base resolution pixel.
//'
//'@details
//'This is a place holder
//'
//'@return 
//'A matrix where column 1 contains the isolattitude ring numbers (i),
//'column 2 contains the pixel-in-ring numbers (j).
//'
//'@name nestIndexC


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

//'//' @export
// [[Rcpp::export]]
NumericVector nestIndexC(int Nside = 0){
 int Npix = 12*Nside*Nside;
 NumericMatrix ind(Npix,3);
 
 int i = 0;
 int j = 0;
 
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
         
       } else if (f > 7 && v < Nside - 1) { //South Polar pixel
         
         int is = 4*Nside - i;
         int ss = (is - Nside + 1) % 2;
         j = (F2*is + h + ss)/2;
         
       } else { //Equatorial pixel
         
         j = (F2*Nside + h + s)/2;
         
       }

      ind(k,0) = i;
      ind(k,1) = j;
      ind(k,2) = k+1;
      k += 1;
    }
 }
 
 return(ind);
}