#include <RcppEigen.h>
#include "systematic.h"
using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
IntegerVector systematic_resampling(const NumericVector & weights, int ndraws, double u){
  RNGScope scope;
  // int nparticles = weights.size();
  IntegerVector ancestors(ndraws);
  u = u / ndraws;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < ndraws; k++){
    while(csw < u){
      j++;
      csw += weights(j);
    }
    u = u + 1. / ndraws;
    ancestors(k) = j + 1;
  }
  return ancestors;
}
