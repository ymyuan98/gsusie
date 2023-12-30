#' @title weighted single effect regression
#'
#' @description
#' The WSER function is to compute the posterior distribution of the regression
#' coefficients of a WSER model.
#'
#'
#' @param weights a (n x 1) array, weights of each subject.
#' For GLM using iterative weighted linear regression,
#' \eqn{weights = exp(-logw2)}
#'
#' @param optimize_V The optimization method to use for fitting the prior
#' variance.
#'
#' @param check_null_threshold Scalar specifying the threshold on the log-scale
#' to compare likelihood between current estimate and zero (null).
#'
#' @return A list with the following elements:
#'
#' \item{alpha}{Vector of posterior inclusion probabilities;
#'  \code{alpha[i]} is posterior probability that the ith coefficient is non-zero.}
#'
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#'
#' \item{mu2}{Vector of posterior second moments (conditional on inclusion).}
#'
#' \item{sigma12}{Vector of posterior variance (conditional on inclusion).}
#'
#' \item{logABF}{Vector of log of asymptotic Bayes factor for each variable.}
#'
#' \item{logABF_model}{Log of asymptotic Bayes factor for the
#'          weighted single effect regression.}
#'
#' \item{V}{Prior variance (after optimization if \code{optimize_V != "none"}).}

#' @keywords internal
#'
weighted_single_effect_regresion <-
    function(y, X,
             weights = 1,
             prior_inclusion_prob = NULL,
             coef_prior_variance = 1,
             optimize_V = c("none", "optim", "EM", "simple"),
             check_null_threshold = 0) {

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

    # Check prior_inclusion_probability
    if (is.null(prior_inclusion_prob)) {
        prior_inclusion_prob <- rep(1 / p, times = p)
    } else {
        if (length(prior_inclusion_prob) != p)
            stop("The length of prior inclusion probability does not match")
    }

    # Update the MLE
    XtWX <- colSums(sweep(X * X, 1, weights, "*"))
    XtWY <- colSums(sweep(X, 1, weights * y, "*"))
    shat2 <- 1 / XtWX
    betahat <- XtWY / XtWX  # = shat2 * XtWY

    optimize_V <- match.arg(optimize_V)  # not applicable yet confirm this part!
    if (optimize_V == "none") {
      V <- coef_prior_variance
    } else if (optimize_V != "EM"){
      V <- optimize_prior_variance(optimize_V, betahat, shat2,
                                   prior_inclusion_prob,
                                   alpha = NULL, post_mean2 = NULL,
                                   V_init = coef_prior_variance,
                                   check_null_threshold = check_null_threshold)
    }

    # compute log-ABF for each variable
    zscore2 <- betahat^2 / shat2
    logABF <- 1 / 2 * ((log(shat2) - log(shat2 + V)) + zscore2 * V / (V + shat2))
    # deal with special case of infinite shat2 (e.g. happens if X does not vary)
    logABF[is.infinite(shat2)] <- 0
    maxlogABF <- max(logABF)

    # w is proportional to ABF, but subtract max for numerical stability
    w <- exp(logABF - maxlogABF)  # w = ABF / ABFmax

    # Update the posterior estimates
    # Posterior prob for each variable
    w_weighted <- w * prior_inclusion_prob
    alpha <- w_weighted / sum(w_weighted)
    post_variance <- 1 / (1/shat2 + 1/V)  # posterior variance
    post_mean <- XtWY / (1/shat2 + 1/V)   # posterior mean
    post_mean2 <- post_variance + post_mean^2               # second moment

    # ABF for WSER model
    logABF_model <- maxlogABF + log(sum(w_weighted))  # = log(sum(ABF x prior_weights))

    if (optimize_V == "EM") {
      V <- optimize_prior_variance(optimize_V, betahat, shat2,
                                   prior_inclusion_prob,
                                   alpha, post_mean2,
                                   check_null_threshold = check_null_threshold)
    }

    return(list(alpha = alpha,
                mu = post_mean,
                mu2 = post_mean2,
                betahat = betahat,
                logABF = logABF,
                logABF_model = logABF_model,
                V = V
                ))
}


#' Estimate prior variance

#' In this function, betahat represents the MLE,
#' and shat2 represents the corresponding variance.
optimize_prior_variance <- function(optimize_V, betahat, shat2,
                                    prior_inclusion_prob,
                                    alpha = NULL, post_mean2 = NULL,
                                    V_init = NULL,
                                    check_null_threshold = 0) {
  V = V_init

  if (optimize_V != "simple") {
    if (optimize_V == "optim") {
      logV <- optim(par = log(max(c(betahat^2 - shat2, 1), na.rm = T)),
                    fn = neg.optimfunc.logscale,
                    betahat = betahat, shat2 = shat2,
                    prior_inclusion_prob = prior_inclusion_prob,
                    method = "Brent", lower = -30, upper = 15)$par
      # If the estimated one is worse than the current one, don't change it
      if (neg.optimfunc.logscale(logV, betahat = betahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob) >
          neg.optimfunc.logscale(log(V), betahat = betahat, shat2 = shat2,
                                 prior_inclusion_prob = prior_inclusion_prob)) {
        logV <- log(V)
      }
      V <- exp(logV)
    }
    else if (optimize_V == "EM") {
      V <- sum(alpha * post_mean2)  # second-order of beta_js (WHY?!!)
    } else
      stop("Invalid option for optimize_V method")
  }
  ## if (optimize_V == "simple"), V is compared with `check_null_threshold`
  ## without any other updates;
  ## if (optimize_V == "none"), V is always the pre-assigned coefficient prior
  ## variance and is not compared with any other values.

  # Following instructions at
  # https://github.com/stephenslab/susieR/blob/master/R/single_effect_regression.R
  # set V exactly 0 if that beats the numerical value by check_null_threshold
  # in loglik. It means that for parsimony reasons we set estimate of V to zero
  # if its numerical estimate is only "negligibly" different from zero.
  if (optimfunc.logscale(0, betahat, shat2, prior_inclusion_prob) +
      check_null_threshold >=
      optimfunc.logscale(V, betahat, shat2, prior_inclusion_prob))
    V <- 0

  return(V)
}


#' The following is the log-scale of the optimization goal
#' as a function of prior variance V.
optimfunc.logscale <- function(V, betahat, shat2, prior_inclusion_prob) {

 # compute log-ABF for each variable
  zscore2 <- betahat^2 / shat2
  logABF <- 1 / 2 * ((log(shat2) - log(shat2 + V)) +
                       zscore2 * V / (V + shat2))
  # deal with special case of infinite shat2 (e.g. happens if X does not vary)
  logABF[is.infinite(shat2)] <- 0
  maxlogABF <- max(logABF)

  # w is proportional to ABF, but subtract max for numerical stability
  w <- exp(logABF - maxlogABF)

  # Update the posterior estimates
  # Posterior prob for each variable
  w_weighted <- w * prior_inclusion_prob
  sum_w_weighted <- sum(w_weighted)

  obj <- log(sum_w_weighted) + maxlogABF  # logABF for a WSER

  return(obj)
}


#' The following is the negative of the objective function
neg.optimfunc.logscale <- function(lV, betahat, shat2, prior_inclusion_prob) {
  return(-optimfunc.logscale(exp(lV), betahat, shat2, prior_inclusion_prob))
}




