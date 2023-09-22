#' @title Functions related to Poisson regression model with log link
#' 
#' @description This file defines functions related to Poisson regression model
#' with logistic link, 
#' including (log of) pseudo-variance and pseudo-response 
#' that are needed to be calculated during iterative 
#' fitting of the (approximated) weighted linear regression model,
#' and the exact and approximated log-likelihood.
#' 
#' The functions below output vectors of the same length of y
#' 
#' #' \code{compute_logw2_poisson} computes the log-pseudo-variance
#' @keywords internal
#' 
compute_logw2_poisson <- function(eta) {
    return(-eta)
}

compute_psdresponse_poisson <- function(eta, y, eta_tol = -50) {
  if (length(eta) != length(y)) 
    stop("Dimensions of input eta and y do not match")
  
  clip_eta <- ifelse(eta < eta_tol, eta_tol, eta)
  
  return(eta + y * exp(-clip_eta) - 1)
  # return(eta + y * exp(-eta) - 1)
}


compute_loglik_exact_poisson <- function(eta, y, eta_tol = 50) {

  if (length(eta) != length(y))
      stop("Dimensions of input eta and y do not match")
  
  clip_eta <- ifelse(eta > eta_tol, eta_tol, eta)
  ## to avoid extremely large exp(eta)
  
  # res <- y * eta - exp(eta) + lfactorial(y)
  res <- y * eta - exp(clip_eta) + lfactorial(y)
  
  return(res)
}
## Always be unexpectedly large!


compute_loglik_apprx_poisson <- function(eta, y) {
    logw2 <- compute_logw2_poisson(eta)
    zz <- compute_psdresponse_poisson(eta, y)  # pseudo-response
    res <- - 1 / 2 * (eta - zz)^2 * exp(-logw2)
    return(res)
}
