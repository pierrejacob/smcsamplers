#include <RcppEigen.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector dmvnorm_(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance){
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::LLT<Eigen::MatrixXd> lltofcov(covariance_);
  Eigen::MatrixXd lower = lltofcov.matrixL();
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = lower.diagonal().array().log().sum();;
  double cst = - (halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * lower.triangularView<Eigen::Lower>().solve(xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}

// [[Rcpp::export]]
NumericVector dmvnorm_cholesky_inverse(const NumericMatrix & x, const NumericVector & mean, const Eigen::MatrixXd & cholesky_inverse){
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = cholesky_inverse.diagonal().array().log().sum();;
  double cst = - (-halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * (cholesky_inverse.transpose() * xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}


// [[Rcpp::export]]
NumericMatrix grad_dmvnorm_precision(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & precision) {
  int N = x.nrow(), d = x.ncol();
  NumericMatrix output(N,d);
  for (int n = 0; n < N; ++n){
    for(int i = 0; i < d; ++i){
      double total = 0;
      for(int j = 0; j <d; ++j){
        total += (mean(j) - x(n,j)) * precision(i,j);
      }
      output(n,i) = total;
    }
  }
  return(output);
}

// [[Rcpp::export]]
NumericMatrix rmvnorm_(int nsamples, const NumericVector & mean, const NumericMatrix & covariance){
  RNGScope scope;
  int ncols = covariance.cols();
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  Eigen::MatrixXd cholesky_covariance(covariance_.llt().matrixU());
  Eigen::MatrixXd Y(nsamples, ncols);
  for(int i = 0; i < ncols; i++){
    Y.col(i) = as<Eigen::ArrayXd>(rnorm(nsamples));
  }
  Y = Y * cholesky_covariance;
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nsamples; i++){
      Y(i,j) = Y(i,j) + mean(j);
    }
  }
  return wrap(Y);
}

// [[Rcpp::export]]
NumericMatrix rmvnorm_cholesky_(int nsamples, const NumericVector & mean, const Eigen::MatrixXd & cholesky){
  RNGScope scope;
  int ncols = cholesky.cols();
  Eigen::MatrixXd Y(nsamples, ncols);
  for(int i = 0; i < ncols; i++){
    Y.col(i) = as<Eigen::ArrayXd>(rnorm(nsamples));
  }
  Y = Y * cholesky;
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nsamples; i++){
      Y(i,j) = Y(i,j) + mean(j);
    }
  }
  return wrap(Y);
}


