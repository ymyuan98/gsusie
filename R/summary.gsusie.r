#' @title Summarize G-SuSiE Fit.
#'
#' @description \code{summary} method for the \dQuote{gsusie} class
#' Almost the same as in
#' <https://github.com/stephenslab/susieR/blob/master/R/summary.susie.R>,
#' except that we change the output class into "summary.gsusie".
#'
#' @param object A GSuSiE fit
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
  class(out) <- c("summary.gsusie", "list")   ## summary.susie -> summary.gsusie
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



