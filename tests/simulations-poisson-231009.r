
## source functions
# library(tidyverse)
library(glmnet)
library(coefplot)
library(readr)

`%&%` <- function(a, b) paste0(a, b)
if.needed <- function(.files, .code) {
  if(!all(file.exists(.files))){
    .code
  }
}

# getwd()
files_path <- "./R/"
r_file_names <- list.files(path = files_path)

for (i in 1 : length(r_file_names)) {
  source(files_path %&% r_file_names[i])
}

################

# get the selected variables
# res: an output from gsusie()
# included_X0: whether the original input data contains intercept X0.
# glmnet_top_n: select the top N variables with the largest absolute values
#  from cv.glmnet. If glmnet_top_n = NULL, it returns all variables with
#  non-zero estimated effect size.

get_selected_vars <- function(res,
                              gs_included_X0 = TRUE,
                              glmnet_top_n = NULL)
{
  if (class(res) == "gsusie") {
    out <- unlist(res$sets$cs)
    if (gs_included_X0) {  # drop X0 if it is selected.
      pp <- ncol(res$mu) - 1
      out <- out[out != (pp+1)]
    }
  }

  if (class(res) == "cv.glmnet") {

    require(coefplot)
    out <- extract.coef(res)[extract.coef(res)[,3] !=   "(Intercept)", ]

    if (!is.null(glmnet_top_n)) {
      out <- out[order(abs(out$Value), decreasing = T), ]
      out <- out[1 : glmnet_top_n, 3]
    } else {
      out <- out[, 3]
    }
   out <- parse_number(out)
  }

  return(out)
}

# generate filenames that store the results of simulations.
# .prefix can be used to add date/"sim"/"real"
gen_filenames <- function(n, p, h2, mod, .prefix = NULL) {
  .filename <- mod %&% "-n" %&% n %&% "-p" %&% p %&% "-h" %&% (10*h2) %&% ".rds"
  if (!is.null(.prefix)) {
    .filename <- .prefix %&% "_" %&% .filename
  }
  return(.filename)
}

################
## simulations

n_trials <- 10

nn <- 500
pp <- 1000
h2 <- 0.6

.save.path <- "./tests/"
.filenames <- gen_filenames(nn, pp, h2, "poisson")

if.needed(.save.path %&% .filenames, {
  set.seed(Sys.Date())

  verbose <- TRUE

  # Precision (PPV)
  PPV <- as.data.frame(matrix(nrow = n_trials, ncol = 6))
  # Recall (TPR)
  TPR <- as.data.frame(matrix(nrow = n_trials, ncol = 6))
  # Number of variables selcted
  NumVs <- as.data.frame(matrix(nrow = n_trials, ncol = 6))
  # Convergence Warnings
  convergeFlags <- as.data.frame(matrix(nrow = n_trials, ncol = 6))
  # ELBOs
  ELBO <- as.data.frame(matrix(nrow = n_trials, ncol = 6))
  ## Name the columns (methods):
  ## vn: vanilla (non-robust)
  ## hb: huber-weighting
  ## bs: bisquare-weighting
  ## fc1: fraction 0.01
  ## gn: glmnet
  ## gn10: glmnet with 10 variables with the largest absolute effect size
  colnames(PPV) <- colnames(TPR) <- colnames(NumVs) <-
    colnames(convergeFlags) <- colnames(ELBO) <-
    c(paste("gs", c("vn", "hb", "bs", "fc1"), sep = "_"), "gn", "gn10")

  for (tt in 1 : n_trials) {
    if (verbose) {
      if (tt %% 10 == 0)
        cat("Running the ", tt, "-th trial... \n")
    }

    #################################
    ## Data generation
    X <- matrix(rnorm(nn * pp), ncol = pp)
    n_effect_vars <- 10  ## number of effective variables
    effect_idx <- sort(sample(1 : pp, size = n_effect_vars))
    bb <- rnorm(n_effect_vars)
    # data.frame(variable = paste0("X", effect_idx), effect_size = bb)

    Eta <- sqrt(h2) * (X[ ,effect_idx, drop = F] %*% as.matrix(bb)) +
      sqrt(1-h2) * rnorm(nn)
    y1 <- rpois(nn, exp(Eta))
    y  <- exp(scale(log1p(y1)))  # scale the response
    X <- cbind(X, 1)  # Intercept is always put at last
    colnames(X) <- paste0("X", c(1:pp, 0))

    #################################
    ## Method comparison

    # Non-robust estimation (vanilla)
    print("gs-vanilla")
    res_gs_vn <- gsusie(X, y, family = "poisson",
                        max_iters = 200,
                        estimate_prior_method = "optim",
                        robust_estimation = F)

    # Huber weighting
    print("gs-huber")
    res_gs_hb <- gsusie(X, y, family = "poisson",
                        max_iters = 200,
                        estimate_prior_method = "optim",
                        robust_estimation = T,
                        robust_method = "huber")

    # Bisquare weighting
    print("gs-bisquare")
    res_gs_bs <- gsusie(X, y, family = "poisson",
                        max_iters = 200,
                        estimate_prior_method = "optim",
                        robust_estimation = T,
                        robust_method = "bisquare")

    # Weight dropped by fractions - 0.01
    print("gs-fraction")
    res_gs_fc1 <- gsusie(X, y, family = "poisson",
                         max_iters = 200,
                         estimate_prior_method = "optim",
                         robust_estimation = T,
                         robust_method = "simple",
                         simple_outlier_fraction = 0.01,
                         simple_outlier_thres = NULL,
                         verbose = T)

    # GLMNET
    res_glmnet <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson")

    ## convergeFlags
    convergeFlags$gs_vn[tt]  <- res_gs_vn$converged
    convergeFlags$gs_hb[tt]  <- res_gs_hb$converged
    convergeFlags$gs_bs[tt]  <- res_gs_bs$converged
    convergeFlags$gs_fc1[tt] <- res_gs_fc1$converged

    ## ELBO
    ELBO$gs_vn[tt]  <- res_gs_vn$elbo[length(res_gs_vn$elbo)]
    ELBO$gs_hb[tt]  <- res_gs_hb$elbo[length(res_gs_hb$elbo)]
    ELBO$gs_bs[tt]  <- res_gs_bs$elbo[length(res_gs_bs$elbo)]
    ELBO$gs_fc1[tt] <- res_gs_fc1$elbo[length(res_gs_fc1$elbo)]


    #####################################
    ## Evaluations - variable selection

    ### Number of variables selected
    NumVs$gs_vn[tt]  <- length(get_selected_vars(res_gs_vn))
    NumVs$gs_hb[tt]  <- length(get_selected_vars(res_gs_hb))
    NumVs$gs_bs[tt]  <- length(get_selected_vars(res_gs_bs))
    NumVs$gs_fc1[tt] <- length(get_selected_vars(res_gs_fc1))
    NumVs$gn[tt]     <- length(get_selected_vars(res_glmnet))
    NumVs$gn10[tt]   <- length(get_selected_vars(res_glmnet, glmnet_top_n = 10)) # 10

    ### Precision
    PPV$gs_vn[tt]  <- evaluate_ppv(get_selected_vars(res_gs_vn),  effect_idx)
    PPV$gs_hb[tt]  <- evaluate_ppv(get_selected_vars(res_gs_hb),  effect_idx)
    PPV$gs_bs[tt]  <- evaluate_ppv(get_selected_vars(res_gs_bs),  effect_idx)
    PPV$gs_fc1[tt] <- evaluate_ppv(get_selected_vars(res_gs_fc1), effect_idx)
    PPV$gn[tt]     <- evaluate_ppv(get_selected_vars(res_glmnet), effect_idx)
    PPV$gn10[tt] <- evaluate_ppv(get_selected_vars(res_glmnet, glmnet_top_n = 10),
                               effect_idx)

    ### Recall
    TPR$gs_vn[tt]  <- evaluate_tpr(get_selected_vars(res_gs_vn),  effect_idx)
    TPR$gs_hb[tt]  <- evaluate_tpr(get_selected_vars(res_gs_hb),  effect_idx)
    TPR$gs_bs[tt]  <- evaluate_tpr(get_selected_vars(res_gs_bs),  effect_idx)
    TPR$gs_fc1[tt] <- evaluate_tpr(get_selected_vars(res_gs_fc1), effect_idx)
    TPR$gn[tt]     <- evaluate_tpr(get_selected_vars(res_glmnet), effect_idx)
    TPR$gn10[tt] <- evaluate_tpr(get_selected_vars(res_glmnet, glmnet_top_n = 10),
                               effect_idx)

  }

  saveRDS(list(Precision = PPV,
               Recall    = TPR,
               NumVs     = NumVs,
               convergeFlags = convergeFlags,
               ELBO = ELBO),
          file = .save.path %&%.filenames)
})
