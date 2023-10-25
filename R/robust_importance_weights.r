#' @title Robust estimations
#'
#' @description
#' In GLM fitting, outliers, i.e. unexpected values, may occur in both
#' observations and intermediate procedures, such as weights
#' (i.e. inverse of pseudo-variance) and pseudo-responses. Hence,
#' robust estimations are be considered to "down-weigh" the influential points.
#'
#' @param values a vector that may contains outliers.
#' May be weights () or residuals(?!)
#'
#' @param method If \code{method = "none"}, then all subjects are assigned
#' (adjuct) weights 1, i.e. all subjects are included in coefficient estimation.
#' If \code{method = "simple"}, then \code{simple_outlier_fraction}\eqn{\times 100%}
#' subjects with the highest (absolute) weights are assigned adjuct weights 0,
#' i.e. those subjects are excluded from coefficient estimation in the current
#' iteration; other subjects are assigned adjuct weights 1.
#' If \code{method = "huber"}, huber weightings are assigned to each subject
#' based on the "pseudo-residuals" \code{rr}.
#' The Huber weight function is defined in
#' [http://users.stat.umn.edu/~sandy/courses/8053/handouts/robust.pdf]
#'
#' @param simple_outlier_fraction a value between 0 and 1, indicating the
#' fraction of outliers: we define the
#' \code{simple_outlier_fraction}\eqn{\times 100%} of subjects with the highest
#' absolute values as outliers
#'
#' @param huber_tuning_k tuning parameter in the Huber weight function.
#' By default, \eqn{\text{huber_tuning_k} = 1.345 \times \hat{\sigma}},
#' where \eqn{\hat{\sigma} = \operatorname{MAR}/0.6745} and
#' \eqn{\operatorname{MAR}} is the median absolute residual.
#'
#' @param bisquare_tuning_k tuning parameter in the Bisquare weight function.
#' By default, \eqn{\text{bisquare_tuning_k} = 4.685 \times \hat{\sigma}}
#' and \eqn{\hat{\sigma}} is estimated using the same equation above.
#'
#' @param previous_imp_weights importance weights in the previous iteration,
#' used to calculate the importance weights in the current iteration for the
#' S-estimation / MM-estimation. If \param{previous_imp_weights} is NULL,
#' it indicates the first iteration, and thus for S-estimation and MM-estimation
#' need to be initialized.
#'
#'
#' @return {imp_weights} The imp weights.
#' If a subject is to be removed from the current iteration, its adjuct_weights
#' is 0.
#'
#' @keywords internal
#'

robust_importance_weights <- function(
    values,
    robust_method = c("none", "huber", "simple", "bisquare"),
    simple_outlier_fraction = NULL,
    simple_outlier_thres = NULL,
    # huber_tuning_k = NULL,
    huber_tuning_k = c("M", "S", "vSD", "vMd"),
    bisquare_tuning_k = NULL,
    previous_imp_weights = NULL  # if is NULL, S-estimation needs to be initialized.
){

  robust_method <- match.arg(robust_method)
  huber_tuning_k <- match.arg(huber_tuning_k)

  if (robust_method == "none") {
    imp_weights <- rep(1, times = length(values))

  }
  else if (robust_method == "simple") {
    if (!is.null(simple_outlier_fraction)) {
      if (simple_outlier_fraction <= 0 | simple_outlier_fraction >= 1)
      {stop("simple_outlier_fraction should be a value between 0 and 1")}

      thres <- quantile(values, 1 - simple_outlier_fraction)
      imp_weights <- 1 * (values <= thres)
    }
    else if(!is.null(simple_outlier_thres)) { # is.null(simple_outlier_fraction)
      thres <- simple_outlier_thres
      imp_weights <- 1 * (values <= thres)
    }
    else { # is.null(simple_outlier_fraction) && is.null(simple_outlier_thres)
      stop("Please specify a threshold to simple-robust estimation!")
    }
  }
  else if (robust_method == "huber"){
    # value: residual

    # if (is.null(huber_tuning_k)) {
    if (huber_tuning_k == "M") {  # M-Estimation
      hat_sigma_r <- median(abs(values - median(values))) / 0.6745
      huber_tuning_k <- 1.345 * hat_sigma_r
      # super conservative
      imp_weights <- huber_weight_(values, huber_tuning_k)
    }
    else if (huber_tuning_k == "S") { # S-Estimation
      if (is.null(previous_imp_weights)) {
        # iteration = 1, initialization
        hat_sigma_r <- median(abs(values - median(values))) / 0.6745
        huber_tuning_k <- 1.345 * hat_sigma_r
        imp_weights <- huber_weight_(values, huber_tuning_k)

      } else { # for iteration > 1
        hat_sigma_r <- sqrt(
          sum(sweep(previous_imp_weights, 1, values^2, "*")) /
            (length(values) * 0.199)
          )
        sdr <- values / hat_sigma_r  # standardized residuals u
        imp_weights <- huber_loss_(sdr, 1.345) / sdr^2
      }
    }
    else if(huber_tuning_k == "vSD") {
      huber_tuning_k <- 1.345 * sd(values)
      imp_weights <- huber_weight_(values, huber_tuning_k)

    }
    else if (huber_tuning_k == "vMd") {
      huber_tuning_k <- 1.345 * median(abs(values))
      imp_weights <- huber_weight_(values, huber_tuning_k)
    }


  } else if (robust_method == "bisquare") {
    if (is.null(bisquare_tuning_k)) {
      bisquare_tuning_k <- 4.685 * median(abs(values)) / 0.6745
    }
    imp_weights <- bisquare_weight_(values, bisquare_tuning_k)

  }

  return(imp_weights)

}

huber_loss_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, r^2/2, k*abs(r) - k^2/2)
  return(out)
}

huber_weight_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, 1, k / abs(r))
  return(out)
}

bisquare_weight_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, (1 - (r/k)^2)^2, 0)
  return(out)
}






