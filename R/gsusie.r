#' @title Generalized Sum of Single-effects (GSuSiE)
#'
#' @description
#' Performs a sparse Bayesian variable selection on generalized
#' linear models using the genaralized sum of single-effects (GSuSiE) model.
#' This function is currently designed for \emph{Poisson} and \emph{logistic}
#' regression models, where the response y can be a count type or binary
#' variable. Codes as well as the arguments and descriptions are referred from
#' [https://github.com/stephenslab/susieR/blob/master/R/susie.R].
#' In brief, this function transforms the GLMs into iterative reweighted
#' least square problems, and then (at each iteration) tries to find the
#' regression coefficients $\beta$ such that the log-likelihood function
#' \eqn{\ell(X,y|\beta,\nu^2)=-\frac{1}{2}\sum_{i=1}^n\left(\frac{z_i-x_{i.}
#' \beta}{\nu}\right)^2} is maximized. Following the \dQuote{susie assumptions},
#' \eqn{\beta = \sum_{l=1}^L \beta_l}, where each \eqn{\beta_l} is a vector of
#' length p with one non-zero element. The value of $L$ is fixed and should be
#' a value of the reasonable upper bound on the number of non-zero effects
#' to be detected.
#'
#' @details
#' To fit GSuSiE on a Poisson regression model, it is recommended to perform
#' robust estimation. The robust estimation helps address the
#' potential negative effect of outliers especially in count data.
#' M/S-estimation with Huber weights is suggested (by setting
#' \code{robust_estimation=TRUE}, \code{robust_method="huber"}, and
#' \code{robust_tuning_method="M"}.
#' To fit GSuSiE on a logistic regression model, it is recommended not to
#' perform robust estimation (by default \code{robust_estimation=FALSE}).
#'
#' @param X An n by p matrix of covariates.
#'
#' @param y The observed responses, a vector of length n.
#'
#' @param family A description of error distribution and link function used in
#' the model. So far, only \emph{binomial distributions with logit link}
#' (\code{family="binomial"}) and \emph{Poisson distributions with log link}
#' (\code{family="poisson"}) are developed.
#'
#' @param maxL Maximum number of non-zero effects in the gsusie model. By
#' default, \code{maxL=10}. If the number of covariates, \emph{p}, is greater
#' than 10, then set \code{maxL=p}.
#'
#' @param prior_inclusion_prob A vector of length p, in which each entry
#'  gives the prior probability that corresponding column of X has a
#'  nonzero effect on the outcome, y.
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
#' @param coef_prior_variance NULL, or a scalar specifying the prior variance
#' of coefficient, or a vector of length \code{maxL}. If it is a scalar, then
#' this value is consistent in all \eqn{\beta_l}.
#'
#' @param check_null_threshold When the prior variance is estimated,
#'  compare the estimated with the null, and set the prior variance to zero
#'  unless the log-likelihood using the estimate is larger by this threshold
#'  amount. For example, if you set \code{check_null_threshold = 0.1},
#'  this will "nudge" the estimate towards zero when the difference in
#'  log-likelihoods is small. A note of caution that setting this to a value
#'  greater than zero may lead the IBSS fitting procedure to occasionally
#'  decrease the ELBO.
#'
#' @param prior_tol When the prior variance is estimated, compare the
#'  estimated value to \code{prior_tol} at the end of the computation,
#'  and exclude a single effect from PIP computation if the estimated prior
#'  vairnace is smaller than this tolerance value.
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
#' @param null_weight Prior probability of no effect (a number between
#'  0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting (or equivalently
#'   standardize XtX and Xty to have the same effect).
#'   Any column of \code{X} that has zero variance is not standardized.
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
#' @param na.rm Drop any missing values in y from both X and y.
#'
#' @param tol A small, non-negative number specifying the convergence
#'   tolerance for the IBSS fitting procedure. The fitting procedure
#'   will halt when the difference in the variational lower bound, or
#'   \dQuote{ELBO} (the objective function to be maximized), is
#'   less than \code{tol}.
#'
#' @param n_purity Passed as argument \code{n_purity} to [gsusie_get_cs]
#'
#' @param abnormal_proportion A value between 0 and 1. Abnormal data point
#' may arise when transforming the GLM into (iterative re)weighted least square
#' problems. Currently, the abnormal points are data points whose
#' pseudo-variances or pseudo-responses are Inf or NAN. If an abnormal point
#' is detected, it is removed from the current iteration.
#' If the number of detected abnormal subjects exceeds
#' \code{abnormal_proportion * nrow(X)}, stop fitting the model.
#'
#' @param track_fit If \code{track_fit=TRUE}, \code{trace} is also returned
#' containing detailed information about the estimates at each iteration of the
#' IBSS fitting procedure.
#'
#' @param verbose If \code{verbose=TRUE}, the algorithm's progress, and a
#' summary of the optimization settings, are printed to the console.
#'
#' @returns A \code{"gsusie"} object with some or all of the following elements:
#'
#' \item{alpha} A \code{maxL} by p matrix of posterior inclusion probabilities.
#'
#' \item{mu} A \code{maxL} by p matrix of posterior means, conditional on
#' inclusion.
#'
#' \item{mu2} A \code{maxL} by p matrix of posterior second moments, conditional
#' on inclusion.
#'
#' \item{lbf} log-Bayes Factor for each single effect.
#'
#' \item{lbf_variable} log-Bayes Factor for each variable and single effect.
#'
#' \item{family} A generalized linear model family, either \dQuote{binomial} or
#' \dQuote{poisson}.
#'
#' \item{Xr} A vector of length n, equal to \eqn{X %\*% colSums(alpha \* mu)}
#'
#' \item{V} Prior variance of the non-zero elements of \eqn{\beta}
#'
#' \item{converged} \code{TRUE} or \code{FALSE} indicating whether the IBSS
#' converged to a solution within the chosen tolerance level.
#'
#' \item{elbo} The value of variational lower bound, or \dQuote{ELBO}, achieved
#' at each iteration of the IBSS fitting procedure.
#'
#' \item{loglik_exact} The value of exact log-likelihood function of the GLM.
#'
#' \item{loglik_apprx} The value of approximated log-likelihood function of the
#' weighted least square model transformed from the corresponding GLM.
#'
#' \item{niter} Number of IBSS iterations that were performed.
#'
#' \item{sets} Credible sets estimated from model fit.
#'
#' \item{pip} A vector of length p giving marginal posterior inclusion
#' probabilities for all p covariates.
#'
#' \item{X_column_scale_factor} A vector of length p giving the scale factor of
#' each column.
#'
#' \item{X_column_center_factor} A vector of length p giving the center factor
#' of each column
#'
#' @examples
#' ## A Poisson regression case ------------------------------------------------
#' set.seed(20231130)
#'
#' # Generative data model
#' nn <- 1000
#' pp <- 10
#' X <- matrix(rnorm(nn * pp), ncol = pp)
#' X[,1:2] <- MASS::mvrnorm(nn, mu = c(0,0),
#'                          Sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2))
#' X[,5:7] <- MASS::mvrnorm(nn, mu = rep(0, 3),
#'                          Sigma = matrix(c(1, 0.6, 0.9,
#'                                           0.6, 1, 0.75,
#'                                           0.9, 0.75, 1),
#'                                         nrow = 3, byrow = T))
#' effect_idx <- c(1, 6)
#' Eta <- scale(X[,effect_idx, drop=F] %*% as.matrix(c(-2, 0.2)))
#' y <- rpois(nn, exp(Eta))
#' plot(y)
#' plot(exp(scale(log1p(y))))
#'
#' ## Vanilla G-SuSiE
#' res_gs <- gsusie(cbind(X, 1), exp(scale(log1p(y))), family = "poisson")
#' summary(res_gs)
#'
#' ## Robust G-SuSiE for Poisson regression
#' res_gs <- gsusie(cbind(X, 1), exp(scale(log1p(y))), family = "poisson",
#' robust_estimation = T, robust_method = "huber", robust_tuning_method = "M")
#' summary(res_gs)
#'
#' ## A logistic regression case -----------------------------------------------
#' set.seed(20240103)
#' # Generative data model
#' nn <- 1000
#' pp <- 10
#' X <- matrix(rnorm(nn * pp), ncol = pp)
#'
#' effect_idx <- c(1, 6)
#' bb <- rnorm(2)
#' Eta <- scale(X[,effect_idx, drop=F] %*% as.matrix(bb))
#' expit <- function(eta) {
#'   ifelse(eta > 0, 1 / (1 + exp(-eta)), exp(eta) / (1 + exp(eta)))
#' }
#' y <- rbinom(nn, 1, expit(Eta))
#' ## Vanilla GSuSiE for logistic regression
#' res_gs <- gsusie(X, y, family = "binomial")
#' summary(res_gs)
#' gsusie_coefficients(res_gs)
#'
#' @export
#'

gsusie <- function(X, y,
                   family = c("binomial", "poisson"),
                   maxL = min(10, ncol(X)),
                   prior_inclusion_prob = NULL,
                   estimate_prior_variance = TRUE,
                   estimate_prior_method = c("optim", "EM", "simple"),
                   coef_prior_variance = 1,
                   check_null_threshold = 0,
                   prior_tol = 1e-9,
                   robust_estimation = FALSE,
                   robust_method = c("huber", "simple", "bisquare"),
                   simple_outlier_fraction = NULL,
                   simple_outlier_thres = NULL,
                   robust_tuning_method = c("M", "S"),
                   null_weight = 0,
                   standardize = TRUE,
                   coverage = 0.95,
                   min_abs_corr = 0.5,
                   max_iters = 100,
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
  model$type <- match.arg(family, c("binomial", "poisson"))
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
                             simple_outlier_thres    = simple_outlier_thres,
                             robust_tuning_method    = robust_tuning_method
                             )

    eta_cur <- compute_Xb(X, colSums(gs$mu * gs$alpha))
    loglik_exact[tt+1] <- sum(model$.compute_loglik_exact(eta_cur, y))
    loglik_apprx[tt+1] <- sum(model$.compute_loglik_apprx(eta_cur, y))
    elbo[tt+1] <- get_objective(X, y, gs, model)

    if (verbose) {
      if (!is.null(gs$abn_subjects)) {
        cat("Abnormal subjects in this round: \n")
        print(gs$abn_subjects)
      }

      cat("ELBO:", elbo[tt+1], "\n")
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

  # Names
  if (!is.null(colnames(X))) {
    variable_names <- colnames(X)
  } else {
    variable_names <- paste0("X", 1:p)
  }
  if (!is.null(null_weight)) {
    variable_names[length(variable_names)] <- "null"
    names(gs$pip) <- variable_names[-p]
  } else {
    names(gs$pip) <- variable_names
  }
  colnames(gs$alpha) <- variable_names
  colnames(gs$mu)    <- variable_names
  colnames(gs$mu2)   <- variable_names
  colnames(gs$lbf_variable) <- variable_names

  # Drop redundant elements
  gs$pie     <- NULL
  gs$sigma02 <- NULL
  # gs$betahat <- NULL

  # For prediction
  gs$X_column_scale_factors  <- attr(X, "scaled:scale")
  gs$X_column_center_factors <- attr(X, "scaled:center")

  return(gs)
}
