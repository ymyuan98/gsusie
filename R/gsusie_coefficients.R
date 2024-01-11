#' @title Summarize posterior estimations of regression coefficients
#'
#' @description
#' This function returns, by default, \code{top_n=10} coefficients with the
#' highest PIPs, sorted in decreasing order.
#'
#'
#' @param object a gsusie fit
#'
#' @param vars an array of numerical indices or variable names
#'    presented in the output. If \code{var_names = NULL},
#'    by default the output contains coefficients of all variables
#'    if \code{ncol(X)<= 10} or \code{top_n} variables with highest PIPs
#'    if \code{ncol(X) > 10}.
#'
#' @param decreasing logical, whether to list the variables
#'    in decreasing order of PIP values.
#'
#' @param top_n numeric, the number of coefficients to be displayed. By default
#' \code{top_n=min(10, length(object$pip))}.
#'
#' @param cred_int logical, whether to return equal-tailed credible intervals
#'
#' @param coverage numeric, between 0 and 1, the coverage probability
#'    of credible intervals.
#'
#' @param digits integer indicating the number of decimal places to be used.
#'    if \code{digits == NULL}, the output will not be rounded.
#'
#' @importFrom stats qnorm
#'
#' @export
#'
gsusie_coefficients <- function(object, vars = NULL,
                                decreasing = TRUE,
                                top_n = 10,
                                cred_int = TRUE,
                                coverage = 0.95,
                                digits = 4) {

  if (is.null(object$pip) || is.null(object$mu) || is.null(object$mu2))
    stop("Cannot print GSuSiE coefficients because ",
         "either PIP, mu, or mu2 is not available.")

  # extract variable names
  if (is.null(names(object$pip))) {
    var_names <- paste0("X", 1:length(object$pip))
  } else {
    var_names <- names(object$pip)
  }

  out <- data.frame(
    PIP       = object$pip,
    post_mean = gsusie_get_posterior_mean(object),
    post_sd   = gsusie_get_posterior_sd(object)
  )
  if (cred_int) {
    if (!is.numeric(coverage) | coverage < 0 | coverage > 1) {
      stop("Input probability should between 0 and 1")
    }

    out$ci_lower <- out$post_mean + qnorm((1 - coverage)/2) * out$post_sd
    out$ci_upper <- out$post_mean - qnorm((1 - coverage)/2) * out$post_sd
  }

  rownames(out) <- var_names

  if (is.null(vars)) {
    if (length(object$pip) <= top_n) {
      out <- out[order(out$PIP, decreasing = decreasing), ]
    } else { # output the top n variables with the highest PIPs
      out <- out[order(out$PIP, decreasing = decreasing)[1:top_n], ]
    }
  } else {
    out <- out[vars, ]
    out[order(out$PIP, decreasing = decreasing), ]
  }

  if (is.null(digits)) {
    return(out)
  } else {
    return(round(out, digits))
  }
}
