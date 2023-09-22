#' @title weighted_single_effect_regression_glm
#' 
#' @description
#' Function \code{weighted_single_effect_regression_glm()} is the weighted single
#' effect regression model dedicated to GLM, in which the the intercept (offset) 
#' term in each simple GLM model is always included. 
#' 
#' When fitting each simple GLM model, 
#' we treat the intercept/offset as an auxiliary variable; 
#' that is, we include the intercept/offset term (with all 1) in the model, 
#' but we do not care about its estimated coefficient. 
#' This intercept/offset term is used to attenuate the "significance" of 
#' the non-effective variables, 
#' as in the pilot tests where we run \code{glm(y ~ -1 + X1, family = ...)},
#' the non-effective variables tend to be biased away from the null,
#' and their estimated regression coefficients are always larger in abs values.
#' 
#' @param coef_prior_variance the prior variance of regression coefficients. 
#' It is a scalar, and is currently designed to be consistent across all 
#' regression coefficients. 
#' REQUIRE MORE FLEXIBILITY!!

# set.seed(Sys.Date())
# 
# nn <- 100
# pp <- 2
# # X <- matrix(rnorm(nn * pp), ncol = pp)
# X <- matrix(rnorm(nn * pp), ncol = pp)
# 
# bb <- 2
# effective_idx <- 1
# 
# Eta <- X[, effective_idx, drop = F] %*% as.matrix(bb)
# y <- rpois(nn, exp(Eta))
# 
# ## Run Poisson-IRLS...
# tol <- 1e-3
# maxIters <- 100
# 
# # loglik_hist <- rep(NA, times = maxIters + 1)
# # loglik_hist[1] <- -Inf
# loglik_raw <- loglik_wlm <- rep(NA, times = 7)
# loglik_raw[1] <- loglik_wlm[1] <- -Inf
# 
# col_idx <- c(1, 2)
# X_input_ <- X[, col_idx, drop=F]
# 
# beta1 <- beta2 <- as.matrix(rep(0, times = length(col_idx)))
# eta1  <- eta2  <- X_input_ %*% beta1

# weighted_single_effect_regression_glm(y, X)

weighted_single_effect_regression_glm <- 
  function(y, X, weights = 1, coef_prior_variance = 1, 
           prior_inclusion_prob = NULL) {
    
  p <- ncol(X)
  
  # Check weights
  if (length(weights) == 1) {
    if (weights != 1 / nrow(X))
      warnings("Assign equal weight to each subject")
    weights <- rep(1 / nrow(X), times = nrow(X))
  } else {
    if (length(weights) != nrow(X))
      stop("Dimensions of input weights and X (1st dim) do not match")
    if (length(weights) != length(y))
      stop("Dimensions of input weights and y do not match")
  }
  weights <- as.numeric(weights)
  
  # Check coefficient prior variance 
  ## REQUIRE MORE FLEXIBILITY!!
  if (length(coef_prior_variance) != 1) {
    stop("Please specify one consistent coefficient prior variance")  
  } 
  
  # Check prior_inclusion_probability
  if (is.null(prior_inclusion_prob)) {
    prior_inclusion_prob <- rep(1 / p, times = p)
  } else {
    if (length(prior_inclusion_prob) != p)
      stop("The length of prior inclusion probability does not match")
  }
  
  
  # Ignore fitted results of the intercept/offset
  shat2 <- betahat <- rep(NA, times = p)
  post_variance <- post_mean <- rep(NA, times = p)
  for (j in 1 : p)  {  ## Need accelaration!!
    X_j <- cbind(1, X[,j, drop = F])  # (n by 2)
    
    XtWX_j <- t(X_j) %*% diag(weights) %*% X_j
    XtWy_j <- t(X_j) %*% diag(weights) %*% y
    
    shat2_j    <- solve(XtWX_j)
    shat2[j]   <- shat2_j[2, 2]
    betahat[j] <- (shat2_j %*% XtWy_j)[2]
    
    postV_j          <- solve(XtWX_j + diag(1/coef_prior_variance, nrow = 2))
    post_variance[j] <- postV_j[2, 2]
    post_mean[j]     <- (postV_j %*% XtWy_j)[2]
  }
  
  # compute log-ABF for each variable
  zscore2 <- betahat^2 / shat2
  logABF <- 1 / 2 * (
    (log(shat2) - log(shat2 + coef_prior_variance)) +
      zscore2 * coef_prior_variance / (coef_prior_variance + shat2)
  )
  # deal with special case of infinite shat2 (e.g. happens if X does not vary)
  logABF[is.infinite(shat2)] <- 0
  maxlogABF <- max(logABF)
  
  # w is proportional to ABF, but subtract max for numerical stability
  w <- exp(logABF - maxlogABF)
  # Update the posterior estimates
  # Posterior prob for each variable
  w_weighted <- w * prior_inclusion_prob
  alpha <- w_weighted / sum(w_weighted)
  # post_variance <- 1 / (1/shat2 + 1/coef_prior_variance)  # posterior variance
  # post_mean <- (1/shat2) * post_variance * betahat         # posterior mean (point estimate)
  post_mean2 <- post_variance + post_mean^2                # second moment
  
  # ABF for single effect model
  logABF_model <- maxlogABF + log(sum(w_weighted))  # = log(sum(ABF x prior_weights))
  
  return(list(alpha = alpha,
              mu = post_mean,
              mu2 = post_mean2,
              betahat = betahat,
              shat2 = shat2,
              logABF = logABF,
              logABF_model = logABF_model
              ))
}


