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
#' @param coef_prior_variance NULL, or a scalar specifying the prior variance
#' of coefficient. Assume \code{coef_prior_variance} is
#' consistent in all \eqn{maxL} submodels if it is not NULL.
#' Besides, if \code{coef_prior_variance} is specified,
#' \code{grid_opt_prior_variance} should be "FALSE".
#'
#' @param grid_opt_prior_variance default FALSE.
#' If \code{grid_opt_prior_variance = TRUE},
#' use grid search (?) to auto-tune (?) the prior variance of coefficient.
#'
#' @param grid_prior_variance_value an array specifying the value of
#' \code{coef_prior_variance} when \code{grid_opt_prior_variance = TRUE}.
#'
#' Default: a fixed value of \code{coef_prior_variance=0.1} with
#' \code{grid_opt_prior_variance = FALSE}
#'
#' @param null_weight Prior probability of no effect (a number between
#'  0 and 1, and cannot be exactly 1).
#'
#' @param standardize If \code{standardize = TRUE}, standardize the
#'   columns of X to unit variance prior to fitting (or equivalently
#'   standardize XtX and Xty to have the same effect). Note that
#'   \code{scaled_prior_variance} specifies the prior on the
#'   coefficients of X \emph{after} standardization (if it is
#'   performed). If you do not standardize, you may need to think more
#'   carefully about specifying \code{scaled_prior_variance}. Whatever
#'   your choice, the coefficients returned by \code{coef} are given for
#'   \code{X} on the original input scale. Any column of \code{X} that
#'   has zero variance is not standardized.
#'
#' @param intercept If \code{intercept = TRUE}, the intercept is
#'   fitted; it \code{intercept = FALSE}, the intercept is set to
#'   zero. Setting \code{intercept = FALSE} is generally not
#'   recommended.
#'
#' But it seems containing/removing intercept does not help finding the
#' effect variables in GLM.
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
#' @param check_null_threshold When the prior variance is estimated,
#'   compare the estimated with the null, and set the prior variance to zero
#'   unless the log-likelihood using the estimate is larger by this threshold
#'   amount. For example, if you set
#'   \code{check_null_threshold = 0.1}, this will "nudge" the estimate
#'   towards zero when the difference in log-likelihoods is small. A
#'   note of caution that setting this to a value greater than zero may
#'   lead the IBSS fitting procedure to occasionally decrease the ELBO.
#'
#' @param aggregate_method Method of aggregating multiple gsusie fits when
#' \code{grid_opt_prior_variance=TRUE}. If \code{aggregate_method="best"},
#' the output is the `gsusie` fit with the maximum ELBO;
#' if \code{aggregate_method="weighted_sum"}, the output is the sum of
#' `gsusie` results fitted under each of \code{grid_prior_variance_value}
#' weighted by their ELBOs.


gsusie <- function(X, y,
                  maxL = min(10, ncol(X)),
                  family = c("binomial", "poisson"),
                  prior_inclusion_prob = NULL,
                  coef_prior_variance = 1,
                  grid_opt_prior_variance = FALSE,
                  grid_prior_variance_value = NULL,
                  null_weight = 0,
                  standardize = TRUE,
                  intercept = FALSE,
                  coverage = 0.95,
                  min_abs_corr = 0.5,
                  max_iters = 500,
                  na.rm = FALSE,
                  check_null_threshold = 0,
                  tol = 1e-2,
                  n_purity = 100,
                  # sensitive_tol = 50,
                  abnormal_proportion = 0.5,
                  # aggregate_method = c("best", "weighted_sum"),
                  track_fit = FALSE,
                  verbose = FALSE) {

    # Check hyperparameter (coef_prior_variance & grid_opt_prior_variance)
    if (grid_opt_prior_variance) {
        if (!is.null(coef_prior_variance)) {
            stop("Please set grid_opt_prior_variance as FALSE ",
                 "or set coef_prior_variance as NULL")
        }
        if (is.null(grid_prior_variance_value)) {
            cat("By default `grid_prior_variance_value` = c(0.1, 1).\n")
            grid_prior_variance_value <- c(0.01, 0.1, 1, 10)
        }
        prior_var_array <- grid_prior_variance_value   ## prior_var is an array
    } else { # grid_opt_prior_variance == FALSE
        if (is.null(coef_prior_variance)) {
            stop("Please specify a coef_prior_variance or change ",
                "grid_opt_prior_variance into TRUE.")
        }
        prior_var_array <- coef_prior_variance   ## prior_var is a scalar
    }

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
            prior_inclusion_prob <- c(rep(1/ncol(X) * (1 - null_weight),ncol(X)), null_weight)
        } else {
            prior_inclusion_prob <- c(prior_inclusion_prob * (1-null_weight), null_weight)
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

    # Check input y.
    n <- nrow(X)

    ## Center and scale input (is it necessary in GLM?)
    # if (intercept)
    #     y <- y - mean_y

    # Set three attributes for matrix X: attr(X,'scaled:center') is a
    # p-vector of column means of X if center=TRUE, a p vector of zeros
    # otherwise; 'attr(X,'scaled:scale') is a p-vector of column
    # standard deviations of X if scale=TRUE, a p vector of ones
    # otherwise; 'attr(X,'d') is a p-vector of column sums of
    # X.standardized^2,' where X.standardized is the matrix X centered
    # by attr(X,'scaled:center') and scaled by attr(X,'scaled:scale').

    colstats <- compute_colstats(X, center = intercept, scale = standardize)
    attr(X, "scaled:center") <- colstats$cm
    attr(X, "scaled:scale")  <- colstats$csd
    if (standardize) {  # standardize the input predictors ## Put it here?
      if (any(attr(X, "scaled:scale") == 0)) {
        intercept_idx <- which(attr(X, "scaled:scale") == 0)
        X[, -intercept_idx] <- scale(X[, -intercept_idx])
      } else {
        intercept_idx <- NULL
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

    gs <- list()  ## aggregation of gs results
    for (prior_var in prior_var_array) {

        gs_res <- init_setup(n, p, maxL, family,
                         prior_inclusion_prob,
                         prior_var,
                         null_weight)

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

          gs_res <- update_each_effect(X, y, gs_res, model, prior_var,
                                       abnormal_proportion) ##!

          eta_cur <- compute_Xb(X, colSums(gs_res$mu * gs_res$alpha)) ##!
          loglik_exact[tt+1] <- sum(model$.compute_loglik_exact(eta_cur, y))
          loglik_apprx[tt+1] <- sum(model$.compute_loglik_apprx(eta_cur, y))
          ## check! exp(-logw2)!!! may result in Inf!
          elbo[tt+1] <- get_objective(X, y, gs_res, model)


          if (verbose) {
            print(tt)
            if (!is.null(gs_res$abn_subjects)) {
              cat("Abnormal subjects in this round: \n")
              print(gs_res$abn_subjects)
              print(paste0("Number of abnormal points: ",
                           length(gs_res$abn_subjects)))
            }
            # print(paste0("objective:", get_objective(X, y, gs_res, model)))
            cat("ELBO:", elbo[tt+1], "\n")
            cat("Loglik:", loglik_exact[tt+1], "\n")
          }


          if (abs(elbo[tt + 1] - elbo[tt]) < tol) {
              gs_res$converged <- TRUE
              break
          }

        }

        elbo <- elbo[2:(tt + 1)]
        loglik_exact <- loglik_exact[2:(tt + 1)]
        loglik_apprx <- loglik_apprx[2:(tt + 1)]

        gs_res$elbo <- elbo
        gs_res$loglik_exact <- loglik_exact
        gs_res$loglik_apprx <- loglik_apprx
        gs_res$niter <- tt


        if (is.null(gs_res$converged)) {
            warning(paste("IBSS algorithm did not converge in", max_iters,
                            "iterations when coef_prior_variance =", prior_var))
            gs_res$converged <- FALSE
        }

        if (track_fit) gs_res$trace <- tracking

        gs_res$abn_subjects <- sort(gs_res$abn_subjects)
        gs <- append(gs, list(gs_res))  # append fitted results as a list

        rm(gs_res)
    }

    # summarize the results
    # gs <- aggregate_gsusie_results(gs, method = aggregate_method)
    gs <- aggregate_gsusie_results(gs)
    gc()


    # GSuSiE CS and PIP
    if (!is.null(coverage) && !is.null(min_abs_corr)) {
        gs$sets <- gsusie_get_cs(gs, coverage = coverage, X = X,
                               min_abs_corr = min_abs_corr,
                               n_purity = n_purity)
        gs$pip <- gsusie_get_pip(gs, prune_by_cs = FALSE)
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
      # if (!is.null(intercept_idx)) {
      #   variable_names <- append(paste0("X", 1:(p-1)), "(Intercept)",
      #                            after = (intercept_idx - 1))
      # } else {
      #   variable_names <- paste0("X", 1:p)
      # }
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
