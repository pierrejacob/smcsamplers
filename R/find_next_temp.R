#' @rdname find_next_temp
#' @title Binary search the next inverse temperature
#' @param gamma_current the current inverse temperature
#' @param ess_criterion the designated effective sample size threshold
#' @param logweights_IS the log ratio of the target density over the proposal density
#' @export
find_next_temp <- function(gamma_current, ess_criterion, logweights_IS){
  ### initialize the proposed inverse temperature to be 1
  gamma_proposal <- 1
  ### compute the effective sample size
  ess_result <- compute_ess(gamma_current, gamma_proposal, logweights_IS)
  if(ess_result$criterion_condition) {
    return(list(ess = ess_result$ess, gamma_next = gamma_proposal, log_ratio_normconst = ess_result$log_ratio_normconst,
                normweights = ess_result$normweights))
  }else{
    ### update the proposed inverse temperature to be the middle point
    gamma_lower <- gamma_current
    gamma_upper <- 1
    gamma_proposal <- (gamma_lower + gamma_upper) / 2
    ### calculate the effective sample size again
    ess_result <- compute_ess(gamma_current, gamma_proposal, logweights_IS)
  }

  while(!ess_result$criterion_condition){
    ess_current <- ess_result$ess
    ### if the effective sample size is larger than the given threshold
    if(ess_current > ess_criterion){
      ### increase the proposed inverse temperature by moving up the lower bound
      gamma_lower <- gamma_proposal
    }else{
      ### otherwise decrease the proposed inverse temperature by moving down the upper bound
      gamma_upper <- gamma_proposal
    }
    ### update the proposed inverse temperature to be the middle point
    gamma_proposal <- (gamma_lower + gamma_upper) / 2
    ### calculate the effective sample size again
    ess_result <- compute_ess(gamma_current, gamma_proposal, logweights_IS)
  }

  return(list(ess = ess_result$ess, gamma_next = gamma_proposal, log_ratio_normconst = ess_result$log_ratio_normconst,
              normweights = ess_result$normweights))
}



