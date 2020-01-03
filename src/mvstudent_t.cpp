#include <RcppEigen.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector dmvstudent_t_cholesky_inverse(const NumericMatrix & x, 
                                            const double & degree, 
                                            const NumericVector & mean, 
                                            const Eigen::MatrixXd & cholesky_inverse, 
                                            const double & constant, 
                                            const double & factor){
  int N = x.nrow();
  int d = x.ncol();
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::MatrixXd xcentered(x_);
  for(int j = 0; j < d; j++){
    for(int i = 0; i < N; i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = (1 / degree) * (cholesky_inverse.transpose() * xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = - factor * log(1 + results(i)) + constant;
  }
  return wrap(results);
}

//' @export
// [[Rcpp::export]]
NumericMatrix grad_dmvstudent_t(const NumericMatrix & x, 
                                const double & degree,
                                const NumericVector & mean, 
                                const NumericMatrix & precision,
                                const double & factor) {
  int N = x.nrow(), d = x.ncol();
  NumericMatrix output(N,d);
  for(int n = 0; n < N; ++n){
    double quadsum = 0;
    for(int i = 0; i < d; ++i){
      double total = 0;
      for(int j = 0; j <d; ++j){
        total += ( x(n,j) - mean(j)) * precision(i,j);
      }
      output(n,i) = total;
      quadsum += total * ( x(n,i) - mean(i) );
    }
    for(int k = 0; k < d; ++k){
      output(n,k) = factor * output(n,k) / (1 + quadsum / degree);  
    }
  }
  return(output);
}
