#include <RcppEigen.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector coxprocess_loglikelihood(const NumericMatrix & x, const IntegerVector & counts, double area){
  int N = x.nrow();
  int d = x.ncol();
  NumericVector output(N);
  for (int n = 0; n < N; n++){
    double cumsum = 0;
    for (int i = 0; i < d; i++){
      cumsum += x(n,i) * counts(i) - area * exp(x(n,i));
    }
    output[n] = cumsum;
  }
  return output;
}

