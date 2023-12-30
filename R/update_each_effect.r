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
                               estimate_prior_variance = FALSE,
                               estimate_prior_method = "optim",
                               check_null_threshold = 0,
                               abnormal_proportion = 0.5,
                               robust_estimation = FALSE,
                               robust_method = "simple",
                               simple_outlier_fraction = 0.01,
                               simple_outlier_thres = NULL,
                               robust_tuning_method = NULL
                              ) {

  if (!estimate_prior_variance) estimate_prior_method = "none"
  if (!robust_estimation)       robust_method = "none"

  if (abnormal_proportion <= 0 | abnormal_proportion >= 1)
      stop("Input abnormal proportion should be a value between 0 and 1.")

  # Update the current linear predictor
  gs$Xr <- compute_Xb(X, colSums(gs$alpha * gs$mu))

  # Update the pseudo-response
  zz <- model$.zz(gs$Xr, y)

  # Update the overall log-pseudo-variance
  llogw2 <- model$.logw2(gs$Xr)
  weights <- exp(-llogw2)    ## May result in Inf or NAN!!
  # check any abnormal points based on log-pseudo-variance
  gs$abn_subjects <- check_abnormal_subjects(weights)

  # Update overall residuals
  rr <- zz - gs$Xr


  # Robust estimation
  if (robust_estimation) {
    if (robust_method == "simple") {
      # "simple" method bases on weights
      gs$imp_weights <- robust_importance_weights(weights,
                              robust_method = "simple",
                              simple_outlier_fraction = simple_outlier_fraction,
                              simple_outlier_thres = simple_outlier_thres)
    } else if (robust_method == "huber") {
      # "huber" method bases on residual rr
      gs$imp_weights <- robust_importance_weights(rr,
                              robust_method = "huber",
                              robust_tuning_method = robust_tuning_method,
                              previous_imp_weights = gs$imp_weights)
    } else if (robust_method == "bisquare") {
      # "bisquare" method bases on residual rr
      gs$imp_weights <- robust_importance_weights(rr,
                              robust_method = "bisquare",
                              robust_tuning_method = robust_tuning_method,
                              previous_imp_weights = gs$imp_weights)
    } else {
      stop("Invalid option for robust_importance_weights method!")
    }

    weights <- sweep(weights, 1, gs$imp_weights, "*")
  }


  # Repeat for each effect to update
  L <- nrow(gs$alpha)

  if (L > 0) {
      # Remove all abnormal points
      # where should I put this part?!
      if (length(gs$abn_subjects) > abnormal_proportion * nrow(X))
          stop("Too much abnormal subjects detected!")
      X_sub       <- remove_abnormal(gs$abn_subjects, X)
      rr_sub      <- remove_abnormal(gs$abn_subjects, rr)
      weights_sub <- remove_abnormal(gs$abn_subjects, weights)

      for (l in 1 : L) {

          # Residuals belonging to Model l
          rrl <- rr_sub + compute_Xb(X_sub, gs$alpha[l, ] * gs$mu[l, ])

          res <- weighted_single_effect_regresion(rrl, X_sub,
                                  weights = weights_sub,
                                  prior_inclusion_prob = gs$pie,
                                  coef_prior_variance = gs$V[l],
                                  optimize_V = estimate_prior_method,
                                  check_null_threshold)

          # Update the variational estimate of the posterior mean
          gs$mu[l, ]           <- res$mu
          gs$mu2[l, ]          <- res$mu2
          gs$alpha[l, ]        <- res$alpha
          gs$betahat[l, ]      <- res$betahat
          gs$lbf_variable[l, ] <- res$logABF
          gs$lbf[l]            <- res$logABF_model
          gs$V[l]              <- res$V

          # Update overall residuals
          rr_sub <- rrl - compute_Xb(X_sub, gs$mu[l, ] * gs$alpha[l, ])
      }
  }

    return(gs)  # return a gsusie model after running one outer loop
}
