#' @title Update each effect once
#' @param X an n by p matrix of regression variables.
#' @param y an n vector of response variable.
#' @param gs a gsusie fit.
#' @param model a list containing functions required for
#' the specific GLM.
#' @param abnormal_proportion a value between 0 and 1.
#' If the number of detected abnormal subjects exceeds
#' \eqn{abnormal_proportion * nrow(X)}, stop fitting the model.
#' 

#' @keywords internal

update_each_effect <- function(X, y, gs, model, 
                               prior_var,
                               abnormal_proportion = 0.5
                              ) {

    if (abnormal_proportion <= 0 | abnormal_proportion >= 1)
        stop("Input abnormal proportion should be a value between 0 and 1.")

    # Update the current linear predictor
    eeta <- compute_Xb(X, colSums(gs$alpha * gs$mu))

    # Update the overall log-pseudo-variance
    llogw2 <- model$.logw2(eeta)
    weights <- exp(-llogw2)    ## May result in Inf or NAN!!
    gs$abn_subjects <- check_abnormal_subjects(weights)

<<<<<<< Updated upstream
    # Update the pseudo-response
    zz <- model$.zz(eeta, y)
    # check any abnormal points based on log-pseudo-responses
=======
  # Update the overall log-pseudo-variance
  llogw2 <- model$.logw2(gs$Xr)
  weights <- exp(-llogw2)    ## May result in Inf or NAN!!
  gs$abn_subjects <- check_abnormal_subjects(weights)
  message(paste(quantiles(weights, c(0, 0.025, 0.975, 1)), collapse = " "))

  gs$abn_subjects <- append(gs$abn_subjects, which(weights > 1000))
  sub_idx <- !((1 : nrow(X)) %in% gs$abn_subjects)
>>>>>>> Stashed changes

    # Update overall residuals
    rr <- zz - eeta

    # Repeat for each effect to update
    L <- nrow(gs$alpha)

    if (L > 0) {

        # Remove all abnormal points
        if (length(gs$abn_subjects) > abnormal_proportion * nrow(X))
            stop("Too much abnormal subjects detected!")
        X_sub      <- remove_abnormal(gs$abn_subjects, X)
        rr_sub     <- remove_abnormal(gs$abn_subjects, rr)
        weights_sub <- remove_abnormal(gs$abn_subjects, weights)

        for (l in 1 : L) {

            # Residuals belonging to Model l
            rrl <- rr_sub + compute_Xb(X_sub, gs$alpha[l, ] * gs$mu[l, ])

            res <- weighted_single_effect_regresion(rrl, X_sub,
                                    weights = weights_sub,
                                    coef_prior_variance = prior_var,
                                    prior_inclusion_prob = gs$pie)

            # Update the variational estimate of the posterior mean
            gs$mu[l, ]           <- res$mu
            gs$mu2[l, ]          <- res$mu2
            gs$alpha[l, ]        <- res$alpha
            gs$betahat[l, ]      <- res$betahat
            gs$lbf_variable[l, ] <- res$logABF
            gs$lbf[l]            <- res$logABF_model

            # Update overall residuals
            rr_sub <- rrl - compute_Xb(X_sub, gs$mu[l, ] * gs$alpha[l, ])

        }
    }

    return(gs)  # return a gsusie model after running one outer loop
}
