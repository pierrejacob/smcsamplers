#' @rdname compute_ess
#' @title Computing the effective sample size
#' @param gamma_current the current inverse temperature
#' @param gamma_proposal the proposed inverse temperature
#' @param logweights_IS the log ratio of the target density over the proposal density
#' @export
compute_ess <- function(gamma_current, gamma_proposal, logweights_IS){
  ### calculate the log importance sampling weights
  logweights <- (gamma_proposal - gamma_current) * logweights_IS
  ### subtract the maximum to overcome numerical issues
  maxlogweights <- max(logweights)
  weights <- exp(logweights - maxlogweights)
  ### normalize the importance sampling weights
  normweights <- weights / sum(weights)
  ### calculate the effective sample size
  ess <- (1 / sum(normweights^2)) / nparticles
  criterion_condition <- FALSE
  ### handle the boundary case separately
  if(gamma_proposal == 1){
    ### check if the ESS is larger than the given threshold
    if(ess > ess_criterion){
      ### calculate the increment of the log normalizing constant estimate
      log_ratio_normconst <- log(mean(weights)) + maxlogweights
      criterion_condition <- TRUE
      return(list(ess = ess, log_ratio_normconst = log_ratio_normconst, normweights = normweights, criterion_condition = criterion_condition))
    }else{
      ### otherwise we shall decrease the proposed inverse temperature
      return(list(criterion_condition = criterion_condition))
    }
    ### handle the intermediate steps
  }else{
    ### binary search the next inverse temperature until the ESS is equal to the given threshold.
    if(abs(ess - ess_criterion) < 10^(-2)){
      ### calculate the increment of the log normalizing constant estimate
      log_ratio_normconst <- log(mean(weights)) + maxlogweights
      criterion_condition <- TRUE
      return(list(ess = ess, log_ratio_normconst = log_ratio_normconst, normweights = normweights, criterion_condition = criterion_condition))
    }else{
      return(list(criterion_condition = criterion_condition, ess = ess))
    }
  }
}
