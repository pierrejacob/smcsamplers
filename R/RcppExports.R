# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @export
coxprocess_loglikelihood <- function(x, counts, area) {
    .Call('_smcsamplers_coxprocess_loglikelihood', PACKAGE = 'smcsamplers', x, counts, area)
}

pgg_m_sigma_ <- function(omega, X, invB, KTkappaplusinvBtimesb) {
    .Call('_smcsamplers_pgg_m_sigma_', PACKAGE = 'smcsamplers', omega, X, invB, KTkappaplusinvBtimesb)
}

logistic_loglikelihood_gradient <- function(betas, Y, X) {
    .Call('_smcsamplers_logistic_loglikelihood_gradient', PACKAGE = 'smcsamplers', betas, Y, X)
}

logistic_loglikelihood_gradient_temperlast <- function(betas, Y, X, delta, gamma) {
    .Call('_smcsamplers_logistic_loglikelihood_gradient_temperlast', PACKAGE = 'smcsamplers', betas, Y, X, delta, gamma)
}

#' @export
eigenMapMatMult <- function(A, B) {
    .Call('_smcsamplers_eigenMapMatMult', PACKAGE = 'smcsamplers', A, B)
}

dmvnorm_ <- function(x, mean, covariance) {
    .Call('_smcsamplers_dmvnorm_', PACKAGE = 'smcsamplers', x, mean, covariance)
}

dmvnorm_cholesky_inverse <- function(x, mean, cholesky_inverse) {
    .Call('_smcsamplers_dmvnorm_cholesky_inverse', PACKAGE = 'smcsamplers', x, mean, cholesky_inverse)
}

grad_dmvnorm_precision <- function(x, mean, precision) {
    .Call('_smcsamplers_grad_dmvnorm_precision', PACKAGE = 'smcsamplers', x, mean, precision)
}

eval_and_grad_dmvnorm_precision <- function(x, mean, precision, cholesky_inverse) {
    .Call('_smcsamplers_eval_and_grad_dmvnorm_precision', PACKAGE = 'smcsamplers', x, mean, precision, cholesky_inverse)
}

rmvnorm_ <- function(nsamples, mean, covariance) {
    .Call('_smcsamplers_rmvnorm_', PACKAGE = 'smcsamplers', nsamples, mean, covariance)
}

rmvnorm_cholesky_ <- function(nsamples, mean, cholesky) {
    .Call('_smcsamplers_rmvnorm_cholesky_', PACKAGE = 'smcsamplers', nsamples, mean, cholesky)
}

#' @export
dmvstudent_t_cholesky_inverse <- function(x, degree, mean, cholesky_inverse, constant, factor) {
    .Call('_smcsamplers_dmvstudent_t_cholesky_inverse', PACKAGE = 'smcsamplers', x, degree, mean, cholesky_inverse, constant, factor)
}

#' @export
grad_dmvstudent_t <- function(x, degree, mean, precision, factor) {
    .Call('_smcsamplers_grad_dmvstudent_t', PACKAGE = 'smcsamplers', x, degree, mean, precision, factor)
}

#' @export
systematic_resampling <- function(weights, ndraws, u) {
    .Call('_smcsamplers_systematic_resampling', PACKAGE = 'smcsamplers', weights, ndraws, u)
}

