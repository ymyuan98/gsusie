#' @title Summarize G-SuSiE Fit.
#'
#' @description \code{summary} method for the \dQuote{gsusie} class
#' Almost the same as in
#' [https://github.com/stephenslab/susieR/blob/master/R/summary.susie.R],
#' except that we change the output class into "summary.gsusie"...(?)
#'
#' @param object A gsusie fit
#'
#' @param \dots Additional arguments passed to the generic \code{summary}
#'  or \code{print.summary} method.
#'
#' @return a list containing a data frame of variables and a data frame
#'  of credible sets.
#'
#' @method summary gsusie
#'
#' @export summary.gsusie
#'
#' @export
#'
summary.gsusie <- function (object, ...) {
  if (is.null(object$sets))
    stop("Cannot summarize GSuSiE object because credible set information ",
         "is not available")
  variables <- data.frame(cbind(1 : length(object$pip), object$pip, -1))
  colnames(variables) <- c("variable", "variable_prob", "cs")
  rownames(variables) <- NULL
  if (object$null_index > 0)
    variables <- variables[-object$null_index,]
  if (!is.null(object$sets$cs)) {
    cs <- data.frame(matrix(NA, length(object$sets$cs), 5))
    colnames(cs) <- c("cs", "cs_loge_abf", "cs_avg_r2", "cs_min_r2", "variable")
      for (i in 1:length(object$sets$cs)) {
        variables$cs[variables$variable %in% object$sets$cs[[i]]] <-
          object$sets$cs_index[[i]]
        cs$cs[i] <- object$sets$cs_index[[i]]
        cs$cs_loge_abf[i] <- object$lbf[cs$cs[i]]
        cs$cs_avg_r2[i] <- object$sets$purity$mean.abs.corr[i]^2
        cs$cs_min_r2[i] <- object$sets$purity$min.abs.corr[i]^2
        cs$variable[i] <- paste(object$sets$cs[[i]],collapse=",")
      }
      variables <- variables[order(variables$variable_prob,decreasing = TRUE),]
  } else
    cs <- NULL
  out <- list(vars = variables,cs = cs)
  class(out) <- c("summary.gsusie", "list")   ##! summary.susie -> summary.gsusie
  return(out)
}

#' @rdname summary.gsusie
#'
#' @param x A (g)susie summary.
#'
#' @method print summary.gsusie
#'
#' @export print.summary.gsusie
#'
#' @export
#'
print.summary.gsusie <- function (x, ...) {
  cat("\nVariables in credible sets:\n\n")
  print.data.frame(x$vars[which(x$vars$cs > 0),],row.names = FALSE)
  cat("\nCredible sets summary:\n\n")
  print.data.frame(x$cs,row.names = FALSE)
}


#' @rdname summary.gsusie
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
#' @param cred_int logical, whether to return equal-tailed credible intervels
#'
#' @param coverage numeric, between 0 and 1, the coverage probability
#'    of credible intervals.
#'
#' @param digits integer indicating the number of decimal places to be used.
#'    if \code{digits == NULL}, the output will not be rounded.
#'
#' @method summarize gsusie coefficients
#'
#' @export coefficients.gsusie
#'
#' @export
#'
gsusie_coefficients <- function(object, vars = NULL,
                                decreasing = TRUE,
                                top_n = 10,
                                cred_int = TRUE,
                                coverage = 0.95,
                                digits = 4, ...) {

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
