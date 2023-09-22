#

# ELBO = Expected loglik + KL(g(b)||q(b)). 
# This function is to compute the second part of ELBO of the lth WSER model 
#
#' @param sigma02 a scalar, prior variance of coefficient 
#' @param pie     a p-dim vector of prior inclusion probabilities
#' @param alpha   a p-dim vector of posterior inclusion probabilities
#' @param mu      a p-dim vector of posterior means
#' @param sigma12 a p-dim vector of posterior variances
#' 
#' @keywords internal
#'
KL_prior_v_post <- function(sigma02, pie, alpha, mu, sigma12) {
    # check dimension
    # if (length(pie) != length(alpha))
    #     stop("Lengths of input pie and alpha do not match")
    # if (length(mu) != length(sigma12) | length(mu) != length(alpha))
    #     stop("Lengths of input mu, sigma12 and alpha do not match")
  
  alpha_   <- clipped(alpha)
  pie_     <- clipped(pie)
  sigma02_ <- clipped(sigma02)
  sigma12_ <- clipped(sigma12)

  res <- alpha_/2 * (1 + log(sigma12_) - log(sigma02_) 
                     - (mu^2 + sigma12_) / sigma02_) + 
         alpha_ * (log(pie_) - log(alpha_))

  # res <- alpha/2 * (1 + log(clipped(sigma12)) - log(clipped(sigma02)) 
  #                   - (mu^2 + sigma12) / sigma02) +
  #         alpha * (log(clipped(pie)) - log(clipped(alpha)))
  return(rowSums(res))
}


#' @description compute the (estimated) ELBO for the overall model.
#' @keywords internal
#' 
get_objective <- function(X, y, gs, model, abn.rm = FALSE) {
    eta <- compute_Xb(X, colSums(gs$alpha * gs$mu))
    
    # check if there is abnormal points.
    if (abn.rm) {
      eta_sub <- remove_abnormal(gs$abn_subjects, eta)
      y_sub <- remove_abnormal(gs$abn_subjects, y)
    } else {
      eta_sub <- eta
      y_sub   <- y
    }
    
    loglik_exact <- sum(model$.compute_loglik_exact(eta_sub, y_sub))
    sigma12 <- gs$mu2 - gs$mu^2
    KL <- sum(KL_prior_v_post(gs$sigma02, gs$pie, gs$alpha, gs$mu, sigma12))
    return(loglik_exact - KL)
}



#' @description clip the inputs to avoid feeding 0 into log()
#' @keywords internal
#' 
clipped <- function(values, tol = 1e-16) {
    ifelse(values > tol, values, tol)
}
