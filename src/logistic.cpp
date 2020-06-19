#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List pgg_m_sigma_(const Eigen::Map<Eigen::MatrixXd>  & omega,
                  const Eigen::Map<Eigen::MatrixXd>  & X,
                  const Eigen::Map<Eigen::MatrixXd>  & invB,
                  const Eigen::Map<Eigen::VectorXd>  & KTkappaplusinvBtimesb){
  int n = X.rows();
  int p = X.cols();
  // The matrix A stores XT Omega X + B^{-1}, that is, Sigma^{-1}
  Eigen::MatrixXd A(p,p);
  for (int j1 = 0; j1 < p; j1 ++){
    for (int j2 = j1; j2 < p; j2 ++){
      A(j1,j2) = invB(j1, j2);
      for (int i = 0; i < n; i++){
        A(j1,j2) = A(j1,j2) + X(i,j1) * X(i,j2) * omega(i);
      }
      A(j2,j1) = A(j1,j2);
    }
  }
  Eigen::LLT<Eigen::MatrixXd> lltofA(A);
  Eigen::MatrixXd lower = lltofA.matrixL();
  Eigen::VectorXd x = lltofA.solve(KTkappaplusinvBtimesb);
  return List::create(Named("m")=x,
                      Named("Sigma_inverse") = A,
                      Named("Cholesky_inverse") = lower,
                      Named("Cholesky") = lower.inverse());
}


// [[Rcpp::export]]
Rcpp::List logistic_loglikelihood_gradient(const Eigen::MatrixXd & betas, const Eigen::ArrayXd & Y, const Eigen::MatrixXd & X){
  int p = X.cols();
  int n = X.rows();
  int M = betas.cols();
  Eigen::MatrixXd xbeta = X * betas;
  Eigen::ArrayXXd expxbeta = xbeta.array().exp();
  Eigen::ArrayXXd gradients(p,M);
  gradients.setZero(p,M);
  Eigen::ArrayXd evals(M);
  // Eigen::MatrixXd tX = X.transpose();
  evals.setZero();
  for (int i = 0; i < n; i ++){
    evals += Y(i) * xbeta.row(i).array() - (1 + expxbeta.row(i)).log();
    for (int j = 0; j < M; j ++){
      gradients.col(j) += X.row(i).array() * (Y(i) - expxbeta(i,j) / (1. + expxbeta(i,j)));
    }
  }
  return List::create(Named("gradients") = wrap(gradients), Named("logls") = wrap(evals));
}

// [[Rcpp::export]]
List logistic_loglikelihood_gradient_temperlast(const Eigen::MatrixXd & betas, const Eigen::ArrayXd & Y, const Eigen::MatrixXd & X, int delta, double gamma){
  int p = X.cols();
  int n = X.rows();
  int M = betas.cols();
  double eval = 0;
  Eigen::MatrixXd xbeta = X * betas;
  Eigen::ArrayXXd expxbeta = xbeta.array().exp();
  Eigen::ArrayXXd gradients(p,M);
  gradients.setZero(p,M);
  Eigen::ArrayXd evals(M);
  // Eigen::MatrixXd tX = X.transpose();
  evals.setZero();
  for (int i = 0; i < (n - delta); i ++){
    evals += Y(i) * xbeta.row(i).array() - (1 + expxbeta.row(i)).log();
    for (int j = 0; j < M; j ++){
      gradients.col(j) += X.row(i).array() * (Y(i) - expxbeta(i,j) / (1. + expxbeta(i,j)));
    }
  }
  for (int i = n - delta; i < n; i ++){
    evals += gamma * (Y(i) * xbeta.row(i).array() - (1 + expxbeta.row(i)).log());
    for (int j = 0; j < M; j ++){
      gradients.col(j) += gamma * (X.row(i).array() * (Y(i) - expxbeta(i,j) / (1. + expxbeta(i,j))));
    }
  }
  return List::create(Named("gradients") = wrap(gradients), Named("logls") = wrap(evals));
}



