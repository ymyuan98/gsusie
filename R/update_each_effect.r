#' @title Update each single effect at a time
#'
#' @param X an n by p matrix of regression variables.
#'
#' @param y an n vector of response variable.
#'
#' @param gs a GSuSiE fit.
#'
#' @param model a list containing functions required for
#' the specific generalized linear model.
#'
#' @param abnormal_proportion a value between 0 and 1.
#' If the number of detected abnormal subjects exceeds
#' \eqn{abnormal_proportion * nrow(X)}, stop fitting the model.
#'
#' @param estimate_prior_variance Boolean. If \code{TRUE}, the coefficient prior
#' variance is estimated and \code{coef_prior_variance} is then used as an
#' initial value for the optimization (if provided). If \code{FALSE}, the prior
#' variance for each of the \code{maxL} effects is fixed and determined by the
#' value supplied to \code{coef_prior_variance}.
#'
#' @param estimate_prior_method The method used for estimating prior variance.
#' When \code{estimate_prior_method="simple"} is used, the likelihood at the
#' specified prior variance is compared to the likelihood at a variance of zero,
#' and the setting with the larger likelihood is retained.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'  compare the estimated with the null, and set the prior variance to zero
#'  unless the log-likelihood using the estimate is larger by this threshold
#'  amount.
#'
#' @param robust_estimation If \code{robust_estimation=TRUE},
#' robust estimation method is applied on the coefficient estimation.
#'
#' @param robust_method If \code{robust_method= "simple"},
#' then the outliers defined by \code{simple_outlier_fraction} or
#' \code{simple_outlier_thres} are simply removed in each iteration.
#' If \code{robust_method="huber"}, Huber weights are additionally
#' multiplied to the pseudo-weights in each iteration.
#' If \code{robust_method = "bisquare"}, (Tukey's) bisquare weights are
#' additionally multipled to the pseudo-weights in each iteration.
#'
#' @param simple_outlier_fraction a value between 0 and 1, indicating the
#' fraction of outliers: we define the
#' \code{simple_outlier_fraction}\eqn{\times 100%} of subjects with the highest
#' absolute values (inverse of pseudo-variance) as outliers.
#' By default, \code{simple_outlier_fraction=NULL} does not set any outlier
#' fraction.
#'
#' @param simple_outlier_thres a real value, indicating the outliers whose
#' inverse of pseudo-variance exceed this threshold to be removed from the
#' current iteration.
#'
#' @param robust_tuning_method If \code{robust_tuning_method = "M"},
#' M-estimation is performed. If \code{robust_tuning_method = "S"},
#' S-estimation is performed.
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
          # gs$betahat[l, ]      <- res$betahat
          gs$lbf_variable[l, ] <- res$logABF
          gs$lbf[l]            <- res$logABF_model
          gs$V[l]              <- res$V

          # Update overall residuals
          rr_sub <- rrl - compute_Xb(X_sub, gs$mu[l, ] * gs$alpha[l, ])
      }
  }

    return(gs)  # return a gsusie model after running one outer loop
}
