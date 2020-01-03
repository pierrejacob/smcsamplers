#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Variance components model on baseball dataset
// Problem settings
const int dimension = 11;
const double nobservations = 10.0;
const double prior_gamma = 0.01;
const double prior_gamma_0 = 4.0;
const double prior_delta = 1.0;
const double inv_prior_delta = 1.0 / prior_delta;
const double prior_alpha = 1.802;

// parameters (beta, theta) are stored as a vector x
// using the ordering x[0] = beta, x[1:10] = theta

// Gamma distribution
double gamma_logdensity(double x, double alpha, double beta){
  double output = alpha * std::log(beta) - std::lgamma(alpha) + (alpha - 1.0) * std::log(x) - beta * x;
  return(output);
}

// Prior
//' @rdname nuclear_logartificialprior
//' @title Evaluate failure counts model artificial prior density on nuclear pumps dataset
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec nuclear_logartificialprior(arma::mat x){
  int N = x.n_rows;
  double logdensity, beta, theta_i;
  arma::vec output(N);
  
  for (int n = 0; n < N; n++){
    beta = x(n, 0);
    logdensity = 0;
    logdensity = gamma_logdensity(beta, prior_gamma_0, prior_delta);
    for (int i = 1; i < dimension; i++){
      theta_i = x(n, i);
      logdensity += gamma_logdensity(theta_i, prior_alpha, beta);
    }
    output(n) = logdensity;
    
  }
  return(output);
}

// Sample from artificial prior distribution
//' @rdname nuclear_sample_artificialprior
//' @title Sample from failure counts model artificial prior distribution
//' @param n number of samples
//' @return samples
//' @export
// [[Rcpp::export]]
arma::mat nuclear_sample_artificialprior(int N){
  arma::mat output(N, dimension);
  double beta, inv_beta;
  
  for (int n = 0; n < N; n++){
    beta = R::rgamma(prior_gamma_0, inv_prior_delta);
    output(n, 0) = beta;
    for (int i = 1; i < dimension; i++){
      inv_beta = 1.0 / beta;
      output(n, i) = R::rgamma(prior_alpha, inv_beta);
    }
  }
  
  return(output);
}

// Likelihood
// load nuclear pump dataset 
arma::vec nuclear_failures(){
  arma::vec output(nobservations);
  output(0) = 5.0;
  output(1) = 1.0;
  output(2) = 5.0;
  output(3) = 14.0;
  output(4) = 3.0;
  output(5) = 19.0;
  output(6) = 1.0;
  output(7) = 1.0;
  output(8) = 4.0;
  output(9) = 22.0;
  return(output);
}
arma::vec dataset_failures = nuclear_failures();

arma::vec nuclear_times(){
  arma::vec output(nobservations);
  output(0) = 94.3;
  output(1) = 15.7;
  output(2) = 62.9;
  output(3) = 126;
  output(4) = 5.24;
  output(5) = 31.4;
  output(6) = 1.05;
  output(7) = 1.05;
  output(8) = 2.1;
  output(9) = 10.5;
  return(output);
}
arma::vec dataset_times = nuclear_times();

//' @rdname nuclear_logartificiallikelihood
//' @title Evaluate failure counts model artificial loglikelihood function on nuclear pumps dataset
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec nuclear_logartificiallikelihood(arma::mat x){
  int N = x.n_rows;
  double loglikelihood, theta_i, beta;
  arma::vec output(N);
  
  for (int n = 0; n < N; n++){
    loglikelihood = 0;
    for (int i = 1; i < dimension; i++){
      theta_i = x(n, i);
      loglikelihood += dataset_failures(i-1) * std::log(theta_i) - theta_i * dataset_times(i-1);
    }
    beta = x(n, 0);
    loglikelihood += (prior_gamma - prior_gamma_0) * std::log(beta);
    output(n) = loglikelihood;
  }
  
  return(output);
}

// Sample from conditonal distribution of beta given theta
//' @rdname nuclear_sample_beta
//' @title Sample from conditonal distribution of beta given theta
//' @param sum_theta vector of length N containg sum of theta values
//' @param lambda an inverse temperature 
//' @return samples
//' @export
// [[Rcpp::export]]
arma::vec nuclear_sample_beta(arma::vec sum_theta, double lambda){
  int N = sum_theta.n_elem;
  arma::vec output(N);
  double shape, rate, scale;
  
  shape = prior_gamma_0 + prior_alpha * nobservations + lambda * (prior_gamma - prior_gamma_0);
  for (int n = 0; n < N; n++){
    rate = prior_delta + sum_theta(n);
    scale = 1.0 / rate;
    output(n) = R::rgamma(shape, scale);
  }
  
  return(output);
}

// Sample from conditonal distribution of theta given beta
//' @rdname nuclear_sample_theta
//' @title Sample from conditonal distribution of theta given beta
//' @param N number of samples
//' @param beta vector of length N 
//' @param lambda an inverse temperature 
//' @return samples
//' @export
// [[Rcpp::export]]
arma::mat nuclear_sample_theta(arma::vec beta, double lambda){
  int N = beta.n_elem;
  arma::mat output(N, nobservations);
  arma::vec shape(nobservations);
  double rate, scale;
  
  for (int k = 0; k < nobservations; k++){
    shape(k) = prior_alpha + lambda * dataset_failures(k);
  }
  
  for (int n = 0; n < N; n++){
    for (int k = 0; k < nobservations; k++){
      rate = beta(n) + lambda * dataset_times(k);
      scale = 1.0 / rate;
      output(n, k) = R::rgamma(shape(k), scale);
    }
  }
  
  return(output);
}
