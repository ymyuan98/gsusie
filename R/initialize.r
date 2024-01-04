#' @title Initialize a GSuSiE object using regression coefficients
#'
# Set default susie initialization
#' @keywords internal
#'
init_setup <- function(n, p, L, family,
                       prior_inclusion_prob,
                       coef_prior_variance,
                       null_weight) {

  if (!is.numeric(coef_prior_variance) || any(coef_prior_variance < 0))
    stop("Prior variance should be (a) positive number(s)")
  if (length(coef_prior_variance) == 1)
    coef_prior_variance <- rep(coef_prior_variance, times = L)  # of length L

  if (is.null(prior_inclusion_prob)) {
      prior_inclusion_prob <- rep(1 / p, p)
  } else {
      if (all(prior_inclusion_prob == 0))
          stop("Prior weight should greater than 0 for at least one variable.")
      prior_inclusion_prob <- prior_inclusion_prob / sum(prior_inclusion_prob)
  }
  if (length(prior_inclusion_prob) != p)
      stop("Prior weights must have length p")
  if (p < L)
      L <- p
  gs = list(alpha   = matrix(1/p, nrow = L, ncol = p),   ## posterior inclusion probability
            mu      = matrix(0, nrow = L, ncol = p),     ## posterior mean estimate
            mu2     = matrix(0, nrow = L, ncol = p),     ## posterior mean^2 estimate
            betahat = matrix(0, nrow = L, ncol = p),     ## MLE estimates
            lbf     = rep(as.numeric(NA), L),
            lbf_variable = matrix(as.numeric(NA), L, p),
            pie     = prior_inclusion_prob,              ## prior inclusion probability
            sigma02 = coef_prior_variance,               ## prior variance of coefficient b.
            family  = family,
            Xr      = rep(0, n),                         ## linear predictors
            V       = coef_prior_variance,               ## estimated prior variance
            abn_subjects = NULL                          ## indices of abnormal subjects
            )
  if (is.null(null_weight)) {
      gs$null_index <- 0
  } else {
      gs$null_index <- p
  }
  class(gs) <- "gsusie"
  return(gs)
}

