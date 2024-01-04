#' @title Robust estimations
#'
#' @description
#' When fitting a generalized linear model, outliers, i.e. unexpected values,
#' may occur in both observations and intermediate procedures, such as weights
#' (i.e. inverse of pseudo-variance) and pseudo-responses. Hence,
#' robust estimations are be considered to "down-weigh" the influential points.
#'
#' @param values a vector that may contains outliers.
#' May be weights () or residuals(?!)
#'
#' @param robust_method If \code{method = "none"}, then all subjects are
#' included in the coefficient estimation. If \code{method = "simple"},
#' then the outliers defined by \code{simple_outlier_fraction} or
#' \code{simple_outlier_thres} are exlucded from the current iteration of
#' coefficient estimation. If \code{method = "huber"}, huber weights are
#' assigned to each subject based on the \dQuote{pseudo-residuals} \code{rr}.
#'
#' @param simple_outlier_fraction a value between 0 and 1, indicating the
#' fraction of outliers: we define the
#' \code{simple_outlier_fraction} percent of subjects with the highest
#' absolute values (inverse of pseudo-variance) as outliers.
#' By default, \code{simple_outlier_fraction=NULL} does not set any outlier
#' fraction.
#'
#' @param simple_outlier_thres a real value, indicating the outliers whose
#' inverse of pseudo-variance exceed this threshold to be removed from the
#' current iteration.
#'
#' @param robust_tuning_method If \code{robust_tuning_method="M"},
#' then M-estimation is performed. If \code{robust_tuning_method="S"},
#' then S-estimation is performed.
#'
#' @param previous_imp_weights importance weights in the previous iteration,
#' used to calculate the importance weights in the current iteration for the
#' S-estimation. If \param{previous_imp_weights} is NULL, it indicates the
#' first iteration, and thus S-estimation need to be initialized.
#'
#' @returns {imp_weights} An vector of length n, the importance weights.
#' If a subject is to be removed from the current iteration,
#' its importance weights is 0.
#'
#' @keywords internal
#'

robust_importance_weights <- function(
    values,
    robust_method = c("none", "huber", "simple", "bisquare"),
    simple_outlier_fraction = NULL,
    simple_outlier_thres = NULL,
    robust_tuning_method = c("M", "S"),
    previous_imp_weights = NULL
){

  robust_method <- match.arg(robust_method)
  robust_tuning_method <- match.arg(robust_tuning_method)

  if (robust_method == "none") {
    imp_weights <- rep(1, times = length(values))

  }

  else if (robust_method == "simple") {
    # value: weight
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
    else {
      stop("Please specify either 'simple_outlier_fraction'",
           "or 'simple_outlier_thres' to perform simple robust estimation!")
    }
  }
  else if (robust_method == "huber"){
    # value: residual
    huber_tuning_k <- 1.345

    if (robust_tuning_method == "M") {  # M-Estimation
      hat_sigma_r <- median(abs(values - median(values))) / 0.6745
      std_resid <- values /  hat_sigma_r  # standardized residuals u
      imp_weights <- huber_weight_(std_resid, huber_tuning_k)
    }
    else if (robust_tuning_method == "S") { # S-Estimation
      if (is.null(previous_imp_weights)) {
        # iteration = 1, initialization
        hat_sigma_r <- median(abs(values - median(values))) / 0.6745
        std_resid <- values /  hat_sigma_r  # standardized residuals u
        imp_weights <- huber_weight_(std_resid, huber_tuning_k)
      }
      else {
        # iteration > 1
        hat_sigma_r <- sqrt(
          sum(sweep(as.matrix(previous_imp_weights), 1, values^2, "*")) /
            (length(values) * 0.199)
          )
        std_resid <- values / hat_sigma_r
        imp_weights <-
          huber_loss_(std_resid, huber_tuning_k) / std_resid^2
      }
    }

  }
  else if (robust_method == "bisquare") {
    bisquare_tuning_k <- 4.685

    if (robust_tuning_method == "M") {
      hat_sigma_r <- median(abs(values - median(values))) / 0.6745
      std_resid <- values / hat_sigma_r
      imp_weights <- bisquare_weight_(std_resid, bisquare_tuning_k)

    } else if (robust_tuning_method == "S") {
      if (is.null(previous_imp_weights)) {
        # iteration = 1, initialization

        # bisquare_tuning_k <- 1.547  ## too strict?!
        hat_sigma_r <- median(abs(values - median(values))) / 0.6745
        std_resid <- values / hat_sigma_r
        imp_weights <- bisquare_weight_(std_resid, bisquare_tuning_k)
      }
      else {  # iteration > 1
        hat_sigma_r <- sqrt(
          sum(sweep(as.matrix(previous_imp_weights), 1, values^2, "*")) /
            (length(values) * 0.199)
        )
        std_resid <- values / hat_sigma_r
        imp_weights <-
          bisquare_loss_(std_resid, bisquare_tuning_k) / std_resid^2
      }
    }

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

bisquare_loss_ <- function(r, k) {
  out <- ifelse(abs(r) <= k,
                k^2/6 * (1 - (1-(r/k)^2)^3),
                k^2 / 6)
  return(out)
}

bisquare_weight_ <- function(r, k) {
  out <- ifelse(abs(r) <= k, (1 - (r/k)^2)^2, 0)
  return(out)
}



