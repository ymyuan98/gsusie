#' @rdname gsusie_get_methods

#' @title Inferences From Fitted G-SuSiE / SuSiE Model

#' @description These functions access basic properties or
#' draw inference from a fitted susie/gsusie model.
#' The codes and descriptions below are copied from
#' [https://github.com/stephenslab/susieR/blob/master/R/susie_utils.R]


#' @return \code{susie_get_objective} returns the evidence lower bound
#' (ELBO) achieved by the fitted susie model and, optionally, at each
#' iteration of the IBSS fitting procedure.cs
#'
#' \code{susie_get_prior_variance} returns the (estimated or fixed)
#' prior variance parameters.
#'
#' \code{susie_get_posterior_mean} returns the posterior mean for the
#' regression coefficients of the fitted susie model.
#'
#' \code{susie_get_posterior_sd} returns the posterior standard
#' deviation for coefficients of the fitted susie model.
#'
#' \code{susie_get_niter} returns the number of model fitting
#' iterations performed.
#'
#' \code{susie_get_pip} returns a vector containing the posterior
#' inclusion probabilities (PIPs) for all variables.
#'
#' \code{susie_get_lfsr} returns a vector containing the average lfsr
#' across variables for each single-effect, weighted by the posterior
#' inclusion probability (alpha).
#'
#' \code{susie_get_posterior_samples} returns a list containing the
#' effect sizes samples and causal status with two components: \code{b},
#' an \code{num_variables} x \code{num_samples} matrix of effect
#' sizes; \code{gamma}, an \code{num_variables} x \code{num_samples}
#' matrix of causal status random draws.
#'
#' \code{susie_get_cs} returns credible sets (CSs) from a susie fit,
#' as well as summaries of correlation among the variables included in
#' each CS. If desired, one can filter out CSs that do not meet a
#' specified \dQuote{purity} threshold; to do this, either \code{X} or
#' \code{Xcorr} must be supplied. It returns a list with the following
#' elements:
#'
#' \item{cs}{A list in which each list element is a vector containing
#'   the indices of the variables in the CS.}
#'
#' \item{coverage}{The nominal coverage specified for each CS.}
#'
#' \item{purity}{If \code{X} or \code{Xcorr} iis provided), the
#'   purity of each CS.}
#'
#' \item{cs_index}{If \code{X} or \code{Xcorr} is provided) the index
#'   (number between 1 and L) of each reported CS in the supplied susie
#'   fit.}
#'
#'
#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_objective <- function(res, last_only = TRUE, warning_tol = 1e-6) {
    if (is.null(res$elbo)) stop("Please specify a case.")

    if (!all(diff(res$elbo) <= (-1 * warning_tol)))
        warning("Objective is decreasing")
    if (last_only)
        return(res$elbo[length(res$elbo)])
    else
        return(res$elbo)
}


#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_log_pseudo_variance <- function(res, X = NULL, scale = TRUE) {
  # Any warning_tol for outliers?

  if (is.null(res$alpha) || is.null(res$mu)) stop("Please specify a case.")
  if (is.null(X)) stop("Please enter data matrix X.")

  model$type <- match.arg(res$family, c("binomial", "poisson"))
  switch(model$type,
         "binomial" = {model$.logw2 <- compute_logw2_logistic},
         "poisson"  = {model$.logw2 <- compute_logw2_poisson}
         )

  if (scale) {
    # TO deal with the all-one column
    X_column_scale_factors <- res$X_column_scale_factors
    X_column_scale_factors[X_column_scale_factors == 0] <- 1

    X_scaled <- sweep(X,        2, res$X_column_center_factors, "-")
    X_scaled <- sweep(X_scaled, 2, X_column_scale_factors,      "/")
  }

  eeta <- compute_Xb(X, colSums(res$alpha * res$mu))
  llogw2 <- model$.logw2(eeta)
  return(llogw2)
}


#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_posterior_mean <- function(res, prior_tol = 1e-9) {

    # Drop the single-effects with estimated prior of zero(?)
    if (is.numeric(res$V)) {
        include_idx <- which(res$V > prior_tol)
    } else {
        include_idx <- 1 : nrow(res$alpha)
    }

    # Extract relevant rows from alpha matrix
    if (length(include_idx) > 0){
      posterior_mean <- colSums((res$alpha * res$mu)[include_idx,,drop = FALSE])
      const_idx <- which(res$X_column_scale_factors == 0)
      posterior_mean[-const_idx] <-
        posterior_mean[-const_idx] / res$X_column_scale_factors[-const_idx]
      return(posterior_mean)
    }
    else
        return(numeric(ncol(res$mu)))
}

#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_posterior_sd <- function(res, prior_tol = 1e-9) {

    # Drop the single-effects with estimated prior of zero(?)
    if (is.numeric(res$V)) {
        include_idx <- which(res$V > prior_tol)
    } else {
        include_idx <- 1 : nrow(res$alpha)
    }

    #Now extract relevent rows from alpha matrix.
    if (length(include_idx) > 0) {
      posterior_var <- colSums((res$alpha * res$mu2 -
                          (res$alpha * res$mu)^2)[include_idx,, drop = FALSE])

      const_idx <- which(res$X_column_scale_factors == 0)
      posterior_var[-const_idx] <-
        posterior_var[-const_idx] / res$X_column_scale_factors[-const_idx]

      return(sqrt(posterior_var))

    } else {
       return(numeric(ncol(res$mu)))
    }
}


#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_niter <- function(res) {
    return(res$niter)
}

#' @rdname gsusie_get_methods
#' @export
#'
gsusie_get_prior_variance <- function(res) {
    res$V
}


#' @rdname gsusie_get_methods
#' @importFrom stats pnorm
#' @export
#'
gsusie_get_lfsr <- function(res) {

    pos_prob <- pnorm(0, mean = t(res$mu),
                         sd = sqrt(res$mu2 - res$mu^2))
    neg_prob <- 1 - pos_prob
    return(1 - rowSums(res$alpha * t(pmax(pos_prob, neg_prob))))
}


#' @rdname gsusie_get_methods
#'
#' @param res A susie fit, an output from \code{\link{susie}}
#' or \code{\link{gsusie}}.
#'
#' @param num_samples The number of draws from the posterior
#' distribution.
#'
#' @importFrom stats rmultinom
#' @importFrom stats rnorm
#'
#' @export
#'
gsusie_get_posterior_samples <- function(res, num_samples) {

    # Remove effects having estimated prior variance equals zero
    if (is.numeric(res$V))
        include_idx <- which(res$V > 1e-9)
    else
        include_idx <- 1 : nrow(res$alpha)

    ## Remember to deal with the all-one column
    X_column_scale_factors <- res$X_column_scale_factors
    X_column_scale_factors[X_column_scale_factors == 0] <- 1
    posterior_mean <- sweep(res$mu, 2,
                            X_column_scale_factors, "/")
    posterior_sd <- sweep(sqrt(res$mu2 - (res$mu)^2), 2,
                            X_column_scale_factors, "/")


    pip <- res$alpha
    L <- nrow(pip)
    p <- ncol(pip)
    bb_samples <- matrix(as.numeric(NA), p, num_samples)
    gamma_samples <- matrix(as.numeric(NA), p, num_samples)
    for (sample_i in 1 : num_samples) {
        bb <- 0
        if (length(include_idx) > 0) {
            for (l in include_idx) {
                gamma_l <- rmultinom(1, 1, pip[l, ])
                effect_size <- rnorm(1,
                                     mean = posterior_mean[l, which(gamma_l != 0)],
                                     sd = posterior_sd[l, which(gamma_l != 0)])
                bb_l <- gamma_l * effect_size
                bb <- bb + bb_l
            }
        }
        bb_samples[, sample_i] <- bb
        gamma_samples[, sample_i] <- as.numeric(bb != 0)
    }
    return(list(bb = bb_samples, gamma = gamma_samples))
}

#' @rdname gsusie_get_methods
#'
#' @param X n by p matrix of values of p variables (covariates) in
#'      n samples. When provided, correlation between variables will be
#'      computed and used to remove CSs whose minimum correlation among
#'      variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables (covariates).
#'      When provided, it will be used to remove CSs whose minimum correlation
#'      among variables is smaller than \code{min_abs_corr}.
#'
#' @param coverage A number between 0 and 1 specifying desired coverage
#'      of each CS.
#'
#' @param min_abs_corr A "purity" threshold for the CS. Any CS that contains
#'      a pair of variables with correlation less than this threshold will be
#'      filtered out and not reported.
#'
#' @param dedup If \code{dedup = TRUE}, remove duplicate CSs.
#'
#' @param squared If \code{squared = TRUE}, report min, mean and median of
#'      squared correlation instead of the absolute correlation.
#'
#' @param check_symmetric If \code{check_symmetric = TRUE}, perform a check
#'      for symmetry of matrix \code{Xcorr} when \code{Xcorr} is provided
#'      (not \code{NULL}).
#'
#' @param n_purity The maximum number of credible set (CS) variables used in
#'      calculating the correlation (\dQuote{purity}) statistic.
#'      When the number of variables included in the CS is greater than
#'      this number, the CS variables are randomly subsampled.
#'
#' @param use_rfast Use the Rfast package for purity calculations.
#'      By default \code{use_rfast = TRUE} if the Rfast package is installed.
#'
#' @export
#'
gsusie_get_cs <- function(res, X = NULL, Xcorr = NULL, coverage = 0.95,
                         min_abs_corr = 0.5, dedup = TRUE, squared = FALSE,
                         check_symmetric = TRUE, n_purity = 100, use_rfast) {

    if (!is.null(X) && !is.null(Xcorr))
        stop("Only one of X or Xcorr should be specified")
    if (check_symmetric) {
        if (!is.null(Xcorr) && !is.symmetric_matrix(Xcorr)) {
            warning_message("X is not symmetric; forcing Xcorr to be symmetric",
                            "by replacing Xcorr with (Xcorr + t(Xcorr))/2")
            Xcorr <- Xcorr + t(Xcorr)
            Xcorr <- Xcorr / 2
        }
    }

    null_index <- 0
    include_idx <- rep(TRUE, nrow(res$alpha))
    if (!is.null(res$null_index)) null_index <- res$null_index
    if (is.numeric(res$V)) include_idx <- (res$V > 1e-9)

    # L x P binary matrix
    status <- in_CS(res$alpha, coverage)

   # L-list of CS positions
    cs <- lapply(1 : nrow(status), function(i) which(status[i, ] != 0))
    claimed_coverage <- sapply(1 : length(cs),
                              function(i) sum(res$alpha[i,][cs[[i]]]))
    include_idx <- include_idx * (lapply(cs, length) > 0)

    # https://github.com/stephenslab/susieR/issues/21
    if (dedup)
        include_idx <- include_idx * (!duplicated(cs))
    include_idx <- as.logical(include_idx)
    if (sum(include_idx) == 0){
        return(list(cs = NULL,
                    coverage = NULL,
                    requested_coverage = coverage))
    }
    cs <- cs[include_idx]
    claimed_coverage <- claimed_coverage[include_idx]

    # Compute and filter by "purity"
    if (missing(use_rfast))
        use_rfast <- requireNamespace("Rfast", quietly = TRUE)
    if (is.null(Xcorr) && is.null(X)) {
        names(cs) <- paste0("L", which(include_idx))
        return(list(cs = cs,
                    coverage = claimed_coverage,
                    requested_coverage = coverage))
    } else {
        purity <- NULL
        for (i in 1 : length(cs)) {
            if (null_index > 0 && null_index %in% cs[[i]]) {
                purity <- rbind(purity, c(-9, -9, -9))
            } else {
                purity <- rbind(purity,
                            matrix(get_purity(cs[[i]], X, Xcorr, squared,
                                      n_purity, use_rfast), 1, 3))
            }
        }
        purity <- as.data.frame(purity)
        if (squared){
            colnames(purity) <-
                c("min.sq.corr", "mean.sq.corr", "median.sq.corr")
        } else {
            colnames(purity) <-
                c("min.abs.corr", "mean.abs.corr", "mean.abs.corr")
        }
        threshold <- ifelse(squared, min_abs_corr^2, min_abs_corr)
        is_pure <- which(purity[, 1] >= threshold)
        if (length(is_pure) > 0) {
            cs        <- cs[is_pure]
            purity    <- purity[is_pure, ]
            row_names <- paste0("L", which(include_idx)[is_pure])
            names(cs) <- row_names
            rownames(purity) <- row_names

            # Re-order CS list and purity rows based on purity.
            ordering <- order(purity[, 1], decreasing = TRUE)

            return(list(cs = cs[ordering],
                        purity = purity[ordering, ],
                        cs_index = which(include_idx)[is_pure[ordering]],
                        coverage = claimed_coverage[ordering],
                        requested_coverage = coverage))
        } else {
           return(list(cs = NULL, coverage = NULL,
                       requested_coverage = coverage))
        }
    }
}

#' @title Get Correlations Between CSs, using Variable with Maximum PIP From Each CS
#'
#' @description This function evaluates the correlation between single effect
#'   CSs. It is not part of the SuSiE inference. Rather, it is designed as
#'   a diagnostic tool to assess how correlated the reported CS are.
#'
#' @param res A SuSiE fit, typically an output from
#'   \code{\link{susie}} or one of its variants.
#'
#' @param X n by p matrix of values of the p variables (covariates) in
#'   n samples. When provided, correlation between variables will be
#'   computed and used to remove CSs whose minimum correlation among
#'   variables is smaller than \code{min_abs_corr}.
#'
#' @param Xcorr p by p matrix of correlations between variables
#'   (covariates). When provided, it will be used to remove CSs whose
#'   minimum correlation among variables is smaller than
#'   \code{min_abs_corr}.
#'
#' @param max When \code{Max = FAFLSE}, return a matrix of CS
#'   correlations. When \code{Max = TRUE}, return only the maximum
#'   absolute correlation among all pairs of correlations.
#'
#' @return A matrix of correlations between CSs, or the maximum
#'   absolute correlation when \code{max = TRUE}.
#'
#' @export
#'
get_cs_correlation <- function(res, X = NULL, Xcorr = NULL, max = FALSE) {
    # if (is.null(res$sets$cs) || length(res$sets$cs) == 1)
    #     return(NA)
    # if (!is.null(X) && !is.null(Xcorr))
    #     stop("Only one of X or Xcorr should be specified")
    # if (!is.null(Xcorr) && !is.symmetric_matrix(Xcorr)) {
    #     warning_message("X is not symmetric; forcing Xcorr to be symmetric",
    #                     "by replacing Xcorr with (Xcorr + t(Xcorr))/2")
    #     Xcorr <- Xcorr + t(Xcorr)
    #     Xcorr <- Xcorr / 2
    # }
    # # Get index for the best PIP per CS
    # max_pip_idx <- sapply(res$sets$cs,
    #                       function(cs) cs[which.max(res$pip[cs])])
    # if (is.null(Xcorr)) {
    #     X_sub <- X[, max_pip_idx]
    #     cs_corr <- muffled_corr(as.matrix(X_sub))
    # } else {
    #     cs_corr <- Xcorr[max_pip_idx, max_pip_idx]
    # }
    # if (max) {
    #     cs_corr <- max(abs(cs_corr[upper.tri(cs_corr)]))
    # }
    # rownames(cs_corr) <- colnames(cs_corr) <- names(res$sets$cs)

  susieR::get_cs_correlation(res, X = X, Xcorr = Xcorr, max = max)
  return(cs_corr)
}

#' @rdname gsusie_get_methods
#'
#' @param prune_by_cs Whether or not to ignore single effects not in
#'   a reported CS when calculating PIP.
#'
#' @param prior_tol Filter out effects having estimated prior variance
#'   smaller than this threshold.
#'
#' @export
#'
gsusie_get_pip <- function(res, prune_by_cs = FALSE, prior_tol = 1e-9) {

    if (inherits(res, c("gsusie", "susie"))) {  # gsusie method applicable

        # Drop null weight columns
        if (!is.null(res$null_index) && res$null_index > 0) {
            res$alpha <- res$alpha[, -res$null_index, drop = FALSE]
        }

        # Drop the single-effects with estimated prior of zero
        if (is.numeric(res$V)) {
            include_idx <- which(res$V > prior_tol)
        } else {
            include_idx <- 1 : nrow(res$alpha)
        }

        # Now extract relevant rows from alpha matrix
        if (length(include_idx) > 0) {
            res <- res$alpha[include_idx, , drop = FALSE]
        } else {
            res <- matrix(0, 1, ncol(res$alpha))
        }

    }
    return(as.vector(1 - apply(1 - res, 2, prod)))
}

# Find how many variables in the CS
# x is a probability vector
#' @keywords internal
n_in_CS_x <- function(x, coverage = 0.9) {
    sum(cumsum(sort(x, decreasing = TRUE)) < coverage) + 1
}

# Return binary vector indicating if each point is in CS.
# x is a probability vector
#' @keywords internal
in_CS_x <- function(x, coverage = 0.9) {
    n <- n_in_CS_x(x, coverage)
    o <- order(x, decreasing = TRUE)
    result <- rep(0, length(x))
    result[o[1 : n]] <- 1
    return(result)
}


# Return an l-by-p binary matrix indicating which variables are
# in susie (gsusie) credible sets
#' @keywords internal
in_CS <- function(res, coverage = 0.9) {
    if (inherits(res, c("susie", "gsusie")))
        res <- res$alpha
    return(t(apply(res, 1, function(x) in_CS_x(x, coverage))))
}

#' @keywords internal
n_in_CS <- function(res, coverage = 0.9) {
    if (inherits(res, c("susie", "gsusie")))
        res <- res$alpha
    return(apply(res, 1, function(x) n_in_CS_x(x, coverage)))
}

# Subsample and compute min, mean, median and max abs corr.
#
#' @importFrom stats median
get_purity <- function(pos, X, Xcorr, squared = FALSE, n = 100,
                       use_rfast) {
    if (missing(use_rfast))
        use_rfast <- requireNamespace("Rfast", quietly = TRUE)
    if (use_rfast) {
        get_upper_tri <- Rfast::upper_tri
        get_median    <- Rfast::med
    } else {
        get_upper_tri <- function(R) R[upper.tri(R)]
        get_median    <- stats::median
    }
    if (length(pos) == 1)
        return(c(1, 1, 1))
    else{
        # Subsample the columns if necessary
        if (length(pos) > n)
            pos <- sample(pos, n)

        if (is.null(Xcorr)) {
            X_sub <- X[, pos]
            X_sub <- as.matrix(X_sub)
            value <- abs(get_upper_tri(muffled_corr(X_sub)))
        } else {
            value <- abs(get_upper_tri(Xcorr[pos, pos]))
        }
        if (squared)
            value <- value^2

        return(c(min(value),
                 sum(value) / length(value),
                 get_median(value)))
    }
}

# Correlation function with specified warning muffled.
#' @importFrom stats cor
#' @keywords internal
muffled_corr <- function(x) {
    withCallingHandlers(cor(x),
            warning = function(w) {
                if (grepl("the standard deviation is zero", w$message))
                    invokeRestart("muffleWarning")
            })
}

# cov2cor function with specified warning muffled
#
#' @importFrom stats cov2cor
#' @keywords internal
muffled_cov2cor <- function(x) {
    withCallingHandlers(cov2cor(x),
        warning = function(w) {
            if (grepl("had 0 or NA entries; non-finite result is doubtful",
                      w$message))
                invokeRestart("muffleWarning")
        })
}

# Check for symmetric matrix
#' @keywords internal
is.symmetric_matrix <- function(x) {
    if (requireNamespace("Rfast", quietly = TRUE))
        return(Rfast::is.symmetric(x))
    else
       return(isSymmetric(x))
}

# Compute standard error for MLE regression coefficients
# under each simple regression model
# Different from that on
# [https://github.com/stephenslab/susieR/blob/master/R/susie_utils.R#L557]
#' @keywords internal
calc_stderr <- function(X, weights) {
  XtWX <- colSums(sweep(X * X, 1, weights, "*"))
  shat2 <- 1 / XtWX
  return(sqrt(shat2))
}

# Return residuals (No such thing!)

# Slim the result of fitted gsusie model
# Different from that on
# [https://github.com/stephenslab/susieR/blob/master/R/susie_utils.R#L557]
#' @keywords internal
gsusie_slim <- function(res) {
    list(alpha = res$alpha, niter = res$niter, mu = res$mu, V = res$V)
}

# Prune single effects to given number L in susie/gsusie model subject.
#' @keywords internal
gsusie_prune_single_effects <- function(gs, L = 0, V = NULL) {
    num_effects <- nrow(gs$alpha)
    if (L == 0) {
        # Filtering will be based on non-zero elements in s$V.
        if (!is.null(gs$V))
            L <- length(which(gs$V > 0))
        else
            L <- num_effects
    }
    if (L == num_effects) {
        gs$sets <- NULL
        return(gs)
    }
    if (!is.null(gs$sets$cs_index))
        effects_rank <- c(gs$sets$cs_index,
                         setdiff(1:num_effects, gs$sets$cs_index))
    else
        effects_rank <- 1:num_effects

    if (L > num_effects) {
        message(paste("Specified number of effects L =",L,
                  "is greater the number of effects",num_effects,
                  "in input (G-)SuSiE model. The (G-)SuSiE model will be expanded",
                  "to have",L,"effects."))

        gs$alpha <- rbind(gs$alpha[effects_rank, ],
                          matrix(1/ncol(gs$alpha), L - num_effects,
                                 ncol(gs$alpha)))
        for (n in c("mu", "mu2", "lbf_variable")) {
            if (!is.null(gs[[n]]))
                gs[[n]] <- rbind(gs[[n]][effects_rank, ],
                                 matrix(0, L - num_effects, ncol(gs[[n]])))
        }
        # for (n in c("KL", "lbf")) {
        for (n in c("lbf")) {
            if (!is.null(gs[[n]]))
                gs[[n]] <- c(gs[[n]][effects_rank], rep(NA, L - num_effects))
        }
        if (!is.null(V)) {
            if (length(V) > 1)
                V[1 : num_effects] <- gs$V[effects_rank]
            else V <- rep(V, L)
        }
        gs$V <- V
    }
    gs$sets <- NULL
    return(gs)
}

#' @title Utility function to display warning messages as they occur
#' @param ... warning message
#' @param style either "warning" or "hint"
#' @importFrom crayon combine_styles
warning_message <- function(..., style = c("warning", "hint")) {
  style <- match.arg(style)
  if (style == "warning" && getOption("warn") >= 0) {
    alert <- combine_styles("bold", "underline", "red")
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- combine_styles("bold", "underline", "magenta")
    message(alert("HINT:"), " ", ...)
  }
}
