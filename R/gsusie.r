#' @title GSuSiE(gsusie)

#' @description This .r file defines a gsusie function.
#' The codes as well as the arguments and descriptions are referenced from
#' [https://github.com/stephenslab/susieR/blob/master/R/susie.R].
#' The main difference is that we include \emph{family} that describes the relationship
#' between the error distribution and link function.
#'

#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param maxL Maximum number of non-zero effects in the susie
#'   regression model. If L is larger than the number of covariates, p,
#'   L is set to p.
#'
#' @param family A description of error distribution and link function
#' to be used in the model (similar to the same argument in \code{glm}).
#' So far, only binomial distribution with logit link and Poisson
#' distribution with log link are developed.
#' !!! Refererence to `glm` codes!!!
#'
#' @param scaled_prior_variance The prior variance, divided by
#'  \code{var(y)} (or by \code{(1/(n-1))yty} for
#'  \code{susie_suff_stat}); that is, the prior variance of each
#'  non-zero element of b is \code{var(y) * scaled_prior_variance}. The
#'  value provided should be either a scalar or a vector of length
#'  \code{L}. If \code{estimate_prior_variance = TRUE}, this provides
#'  initial estimates of the prior variances.
#'
#' @param prior_inclusion_prob A vector of length p, in which each entry
#'  gives the prior probability that corresponding column of X has a
#'  nonzero effect on the outcome, y.
#'
#' @param estimate_prior_variance If \code{estimate_prior_variance = TRUE},
#' the coefficient prior variance is estimated (this is a hyperparameter for
#' each of the L effects). If provided, (the first element of)
#' \code{coef_prior_variance} is then used as an initial value for the
#' optimization (if it is a vector). When \code{estimate_prior_variance = FALSE},
#' the prior variance for each of the L effects is determined by the value
#' supplied to \code{coef_prior_variance}.
#'
#' @param estimate_prior_method
#'
#' @param coef_prior_variance NULL, or a scalar specifying the prior variance
#' of coefficient, or a vector of length \code{maxL}. If
#' \code{coef_prior_variance} is a scalar, assume this value is consistent in
#' all \code{maxL} sub-models.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'  compare the estimated with the null, and set the prior variance to zero
#'  unless the log-likelihood using the estimate is larger by this threshold
#'  amount. For example, if you set
#'  \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'  towards zero when the difference in log-likelihoods is small. A
#'  note of caution that setting this to a value greater than zero may
#'  lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'  estimated value to \code{prior_tol} at the end of the computation,
#'  and exclude a single effect from PIP computation if the estimated prior
#'  vairnace is smaller than this tolerance value.
#'
#'  #' Check the logic of parameters related to coef_prior_variance
#'
#' @param null_weight Prior probability of no effect (a number between
#'  0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting (or equivalently
#'   standardize XtX and Xty to have the same effect). Any column of \code{X}
#'   that has zero variance is not standardized.
#'
#' @param coverage A number between 0 and 1 specifying the
#'   \dQuote{coverage} of the estimated confidence sets.
#'
#' @param min_abs_corr Minimum absolute correlation allowed in a
#'   credible set. The default, 0.5, corresponds to a squared
#'   correlation of 0.25, which is a commonly used threshold for
#'   genotype data in genetic studies.
#'
#' @param max_iters Maximum number of iterations.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param abnormal_proportion a value between 0 and 1.
#' If the number of detected abnormal subjects exceeds
#' \eqn{abnormal_proportion * nrow(X)}, stop fitting the model.
#'
#' @param robust_estimation If \code{robust_estimation=TRUE},
#' robust estimation method is applied on the coefficient estimation.
#'
#' @param robust_method If \code{robust_est_method = "simple"},
#' then the 1% observations with the highest weights (i.e. inverse of
#' pseudo-variance) during coefficient estimation in each iteration.
#' If \code{robust_method = "huber"}, huber weightings are additionally
#' multiplied to the original weights in each iteration. (TBC?!)


gsusie <- function(X, y,
                  maxL = min(10, ncol(X)),
                  family = c("binomial", "poisson"),
                  prior_inclusion_prob = NULL,
                  estimate_prior_variance = TRUE,
                  estimate_prior_method = c("optim", "EM", "simple"),
                  coef_prior_variance = 1,
                  check_null_threshold = 0,
                  prior_tol = 1e-9,
                  robust_estimation = FALSE,
                  robust_method = c("simple", "huber"),
                  simple_outlier_fraction = 0.01,
                  simple_outlier_thres = NULL,
                  huber_tuning_k = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  max_iters = 500,
                  na.rm = FALSE,
                  tol = 1e-2,
                  n_purity = 100,
                  abnormal_proportion = 0.5,
                  track_fit = FALSE,
                  verbose = FALSE) {

  estimate_prior_method = match.arg(estimate_prior_method)
  robust_method = match.arg(robust_method)

  # Check input X
  if (!(is.double(X) & is.matrix(X)) & !inherits(X, "CsparseMatrix") &
      is.null(attr(X,"matrix.type"))) {
          stop("Input X must be a double-precision matrix, or a sparse matrix, or ",
              "a trend filtering matrix")
      }

  if (is.numeric(null_weight) && null_weight == 0) null_weight <- NULL

  if (!is.null(null_weight) && is.null(attr(X, "matrix.type"))) {
      if (!is.numeric(null_weight))
          stop("Null weight must be numeric")
      if (null_weight < 0 || null_weight >= 1)
          stop("Null weight must be between 0 and 1")
      if (missing(prior_inclusion_prob)) {
          prior_inclusion_prob <- c(rep(1/ncol(X) * (1 - null_weight),ncol(X)),
                                    null_weight)
      } else {
          prior_inclusion_prob <- c(prior_inclusion_prob * (1-null_weight),
                                    null_weight)
      }
      X <- cbind(X, 0)
  }
  if (anyNA(X)) stop("Input X must not contain missing values")
  if (anyNA(y)) {
      if (na.rm) {
          samples_kept <- which(!is.na(y))
          y <- y[samples_kept]
          X <- X[samples_kept, ,drop = F]
      } else
          stop("Input y must not contain missing values")
  }
  p <- ncol(X)
  n <- nrow(X)

  # Set two attributes for matrix X: attr(X,'scaled:center') is a
  # p-vector of column means of X if center=TRUE, a p vector of zeros
  # otherwise; 'attr(X,'scaled:scale') is a p-vector of column
  # standard deviations of X if scale=TRUE, a p vector of ones
  # otherwise.
  colstats <- compute_colstats(X, center = standardize, scale = standardize)
  attr(X, "scaled:center") <- colstats$cm
  attr(X, "scaled:scale")  <- colstats$csd
  if (standardize) {  # standardize the input predictors
    if (any(attr(X, "scaled:scale") == 0)) {
      const_idx <- which(attr(X, "scaled:scale") == 0)
      X[, -const_idx] <- scale(X[, -const_idx])
    } else {
      const_idx <- NULL
      X <- scale(X)
    }
  }

  # Check GLM family
  model <- list()
  model$type <- match.arg(family, c("binomial", "poisson"))  # weighted LeastSquare?! ADD THIS OPTION!
  switch(model$type,
        "binomial" =
        {
            model$.logw2 <- compute_logw2_logistic
            model$.zz    <- compute_psdresponse_logistic
            model$.compute_loglik_exact  <- compute_loglik_exact_logistic
            model$.compute_loglik_apprx  <- compute_loglik_apprx_logistic
        },
        "poisson" =
        {
            model$.logw2 <- compute_logw2_poisson
            model$.zz    <- compute_psdresponse_poisson
            model$.compute_loglik_exact  <- compute_loglik_exact_poisson
            model$.compute_loglik_apprx  <- compute_loglik_apprx_poisson
        })

  # Check coef_prior_variance
  if (!is.null(coef_prior_variance)) {
    if (length(coef_prior_variance) == 1 ||
        length(coef_prior_variance) == maxL) {
      prior_var <- coef_prior_variance
    } else {
      stop("Length of input coef_prior_variance should be either 1 or maxL!")
    }
  }

  gs <- init_setup(n, p, maxL, family, prior_inclusion_prob,
                   prior_var, null_weight)
  # Initialize elbo to NA
  elbo <- rep(as.numeric(NA), max_iters + 1)
  elbo[1] <- -Inf
  #Initialize loglik_exact and loglik_apprx to NA
  loglik_exact <- loglik_apprx <- rep(as.numeric(NA), max_iters + 1)
  loglik_exact[1] <- loglik_apprx[1] <- -Inf

  for (tt in 1 : max_iters){
    if (track_fit) {
      tracking <- list()
      tracking[[tt]] <- gsusie_slim(gs_res)
    }

    gs <- update_each_effect(X, y, gs, model,
                             estimate_prior_variance = estimate_prior_variance,
                             estimate_prior_method   = estimate_prior_method,
                             check_null_threshold    = check_null_threshold,
                             abnormal_proportion     = abnormal_proportion,
                             robust_estimation       = robust_estimation,
                             robust_method           = robust_method,
                             simple_outlier_fraction = simple_outlier_fraction,
                             simple_outlier_thres = simple_outlier_thres,
                             huber_tuning_k = huber_tuning_k
                             )

    eta_cur <- compute_Xb(X, colSums(gs$mu * gs$alpha))
    loglik_exact[tt+1] <- sum(model$.compute_loglik_exact(eta_cur, y))
    loglik_apprx[tt+1] <- sum(model$.compute_loglik_apprx(eta_cur, y))
    elbo[tt+1] <- get_objective(X, y, gs, model)

    if (verbose) {
      print(tt)

      if (!is.null(gs$abn_subjects)) {
        cat("Abnormal subjects in this round: \n")
        print(gs$abn_subjects)
        print(paste0("Number of abnormal points: ", length(gs$abn_subjects)))
      }

      cat("ELBO:", elbo[tt+1], "\n")
      cat("Loglik:", loglik_exact[tt+1], "\n")
    }

    if (abs(elbo[tt + 1] - elbo[tt]) < tol) {
        gs$converged <- TRUE
        break
    }
  }

  elbo <- elbo[2:(tt + 1)]
  loglik_exact <- loglik_exact[2:(tt + 1)]
  loglik_apprx <- loglik_apprx[2:(tt + 1)]

  gs$elbo <- elbo
  gs$loglik_exact <- loglik_exact
  gs$loglik_apprx <- loglik_apprx
  gs$niter <- tt

  if (is.null(gs$converged)) {
      warning(paste("IBSS algorithm did not converge in",
                    max_iters, "iterations"))
      gs$converged <- FALSE
  }
  if (track_fit) gs$trace <- tracking
  gs$abn_subjects <- unique(sort(gs$abn_subjects))

  # GSuSiE CS and PIP
  if (!is.null(coverage) && !is.null(min_abs_corr)) {
      gs$sets <- gsusie_get_cs(gs, coverage = coverage, X = X,
                             min_abs_corr = min_abs_corr,
                             n_purity = n_purity)
      gs$pip <- gsusie_get_pip(gs, prune_by_cs = FALSE, prior_tol = prior_tol)
  }

  if (!is.null(colnames(X))) {
    variable_names <- colnames(X)
    if (!is.null(null_weight)) {
      variable_names[length(variable_names)] <- "null"
      names(gs$pip) <- variable_names[-p]
    } else {
      names(gs$pip) <- variable_names
    }
  } else {
    variable_names <- paste0("X", 1:p)
    if (!is.null(null_weight)) {  ## Why it is here?
      variable_names[length(variable_names)] <- "null"
      names(gs$pip) <- variable_names[-p]
    } else {
      names(gs$pip) <- variable_names
    }
  }
  colnames(gs$alpha) <- variable_names
  colnames(gs$mu)    <- variable_names
  colnames(gs$mu2)   <- variable_names
  colnames(gs$betahat)      <- variable_names
  colnames(gs$lbf_variable) <- variable_names

  # For prediction
  gs$X_column_scale_factors  <- attr(X, "scaled:scale")
  gs$X_column_center_factors <- attr(X, "scaled:center")

  return(gs)
}
