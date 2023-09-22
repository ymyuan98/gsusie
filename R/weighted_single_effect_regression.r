#' @title weighted single effect regression
#' 
#' @param weights a (n x 1) array, weights of each subject. 
#' For GLM using iterative weighted linear regression, 
#' \eqn{weights = exp(-logw2)}
#' 
#' @return A list with the following elements: 
#' 
#' \item{alpha}{Vector of posterior inclusion probabilities; 
#'       \code{alpha[i]} is posterior probability that the ith coefficient is non-zero.}
#' 
#' \item{mu}{Vector of posterior means (conditional on inclusion).}
#' 
#' \item{mu2}{Vector of posterior second moments (conditional on inclusion).} 
#' 
#' \item{sigma12}{Vector of posterior variance (conditional on inclusion).} 
#'  
#' \item{logABF}{Vector of log of asymptotic Bayes factor for each variable.}
#' 
#' \item{logABF_model}{log of asymptotic Bayes factor for the 
#'          weighted single effect regression}
#' 
#' Any additional output? 

#' @keywords internal
#'
weighted_single_effect_regresion <-
    function(y, X, 
             weights = 1, 
             coef_prior_variance = 1,
             prior_inclusion_prob = NULL) {

    p <- ncol(X)

    # optimize_V <- match.arg(optimize_V)  # not applicable yet

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
    # XtWX <- colSums(X * X * weights)     # sweep
    # XtWY <- colSums(X * c(weights * y))  # sweep
    XtWX <- colSums(sweep(X * X, 1, weights, "*"))
    XtWY <- colSums(sweep(X, 1, weights * y, "*"))
    shat2 <- 1 / XtWX
    betahat <- XtWY / XtWX  # = shat2 * XtWY

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
    post_variance <- 1 / (1/shat2 + 1/coef_prior_variance)  # posterior variance
    # post_mean <- (1/shat2) * post_variance * betahat      # posterior mean
    post_mean <- XtWY / (1/shat2 + 1/coef_prior_variance)   # posterior mean
    post_mean2 <- post_variance + post_mean^2               # second moment

    # ABF for single effect model
    logABF_model <- maxlogABF + log(sum(w_weighted))  # = log(sum(ABF x prior_weights))

    return(list(alpha = alpha,
                mu = post_mean,
                mu2 = post_mean2,
                betahat = betahat,
                logABF = logABF,
                logABF_model = logABF_model
                ))
}
