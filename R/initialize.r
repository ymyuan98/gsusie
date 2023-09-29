#' @title Initialize a gsusie object using regression coefficients (?)
#'
#'
# Set default susie initialization
#' @keywords internal
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
            V       = coef_prior_variance,               ## estimated prior variance(?!)
            abn_subjects = NULL                          ## indices of abnormal subjects
            #   logw2   = rep(0, n),   ## monitor (pseudo-)variance of each subject
            )
  if (is.null(null_weight)) {
      gs$null_index <- 0
  } else {
      gs$null_index <- p
  }
  class(gs) <- "gsusie"
  return(gs)
}


init_finalize <- function(gs, X = NULL, Xr = NULL) {
    # if (length(gs$V) == 1)
    #     gs$V = rep(gs$V, nrow(gs$alpha))   ###(?)  How to initialize V?!

    # if (nrow(gs$alpha) != length(gs$V))
    #     stop("Input prior variance must have length of nrow of alpha in ",
    #         "input object")

    # Check prior variance
    # if (!is.numeric(gs$V))
    #     stop("Input prior variance must be numeric")
    # if (!all(gs$V >= 0))
    #     stop("Prior variance must be non-negative")
    if (!all(dim(gs$mu) == dim(gs$mu2)))
        stop("Dimensions of mu and mu2 in input object do not match")
    if (!all(dim(gs$mu) == dim(gs$alpha)))
        stop("Dimensions of mu and alpha in input object do not match")

    # Update Xr: linear predictor.
    if (!missing(Xr))
        gs$Xr <- Xr
    if (!missing(X))
        gs$Xr <- compute_Xb(X, colSums(gs$mu * gs$alpha))

    # Reset KL and lbf
    gs$lbf <- rep(as.numeric(NA), nrow(gs$alpha))
    class(gs) <- "gsusie"
    return(gs)
}
