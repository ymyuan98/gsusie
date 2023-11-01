
## source functions
.libPaths("~/projects/Rlibs/")
library(glmnet)
# library(coefplot)
library(readr)
library(PRROC)

`%&%` <- function(a, b) paste0(a, b)
if.needed <- function(.files, .code) {
  if(!all(file.exists(.files))){
    .code
  }
}

# getwd()
files_path <- "~/projects/gsusie/R/"
# files_path <- "./R/"
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
  if (is.null(res)) {  # if an error is raised...
    return(NULL)
  }

  if (class(res) %in% c("gsusie", "susie")) {
    out <- unlist(res$sets$cs)
    if (gs_included_X0) {  # drop X0 if it is selected.
      pp <- ncol(res$mu) - 1
      out <- out[out != (pp+1)]
    }
  }

  if (class(res) == "cv.glmnet") {

    out <- coef(res, s = "lambda.min")[-1]

   if (is.null(glmnet_top_n)) glmnet_top_n <- length(out)
   out <- order(abs(out), decreasing = T)[1:glmnet_top_n]
  }

  return(out)
}

# generate filenames that store the results of simulations.
# .prefix can be used to add date/"sim"/"real"
gen_filenames <- function(n, p, h2, mod, .prefix = NULL) {
  if (h2 < 0.1 & h2 >= 0.01) {
    .filename <- mod %&% "-n" %&% n %&% "-p" %&% p %&% "-hz" %&% (100*h2) %&% ".rds"
  } else {
    .filename <- mod %&% "-n" %&% n %&% "-p" %&% p %&% "-h" %&% (10*h2) %&% ".rds"
  }

  if (!is.null(.prefix)) {
    .filename <- .prefix %&% "-" %&% .filename
  }
  return(.filename)
}

################
## simulations

n_trials <- 10

## read in arguments
args <- commandArgs(trailingOnly = TRUE)
nn <- as.numeric(args[1])
pp <- as.numeric(args[2])
h2 <- as.numeric(args[3])
# nn <- 500
# pp <- 1000
# h2 <- 0.5

model <- "poisson"
.result.dir <- "~/projects/gsusie/results/synthetic-" %&% Sys.Date() %&% "-" %&% model %&% "/"
# .result.dir <- "./results/synthetic-" %&% Sys.Date() %&% "-" %&% model %&% "/"
if (!dir.exists(.result.dir)) {
  dir.create(.result.dir, recursive=TRUE, showWarnings=FALSE)
}

.filenames <- gen_filenames(nn, pp, h2, model, Sys.Date())

if.needed(.result.dir %&% .filenames, {

  verbose <- TRUE

  # methods-names
  ## Name the columns (methods):
  ## vn: vanilla (non-robust)
  ## hb: huber-weighting
  ## bs: bisquare-weighting
  ## fc30: fraction 0.5
  ## fc5: fraction 0.05
  ## fc1: fraction 0.01
  ## la:   lasso
  ## la10: lasso with 10 variables with the largest absolute effect size
  ## rg:   ridge regression
  ## rg10: ridge with 10 variables ...
  ## en:   elastic net
  ## en10: elastic net with 10 variables
  methods_names <- c(paste("gs", c("vn", "hb_S", "hb_M", "hb_vSD", "hb_vMd",
                                   "bs", "fc30", "fc5", "fc1"),
                           sep = "_"),
                     paste0(rep(c("la", "rg", "en"), each=2),
                            rep(c("", "10"), times = 3)))

  # Precision (PPV)
  PPV <- as.data.frame(matrix(nrow = n_trials, ncol = length(methods_names)))
  # Recall (TPR)
  TPR <- as.data.frame(matrix(nrow = n_trials, ncol = length(methods_names)))
  # Number of variables selcted
  NumVs <- as.data.frame(matrix(nrow = n_trials, ncol = length(methods_names)))
  # AUPRC
  AUPRC <- as.data.frame(matrix(nrow = n_trials, ncol = length(methods_names)))
  # Convergence Warnings (gsusie)
  convergeFlags <- as.data.frame(matrix(nrow = n_trials,
                                ncol = length(startsWith(methods_names, "gs"))))
  # ELBOs (gsusie)
  ELBO <- as.data.frame(matrix(nrow = n_trials,
                               ncol = length(startsWith(methods_names, "gs"))))

  ##
  colnames(PPV) <- colnames(TPR) <- colnames(NumVs) <- colnames(AUPRC) <-
    methods_names
  colnames(convergeFlags) <- colnames(ELBO) <- methods_names[1:ncol(ELBO)]


  for (tt in 1 : n_trials) {
    set.seed(20231020+tt)

    if (verbose) {
      # if (tt %% 10 == 0)
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
    tryCatch({
      res_gs_vn <- gsusie(X, y, family = "poisson",
                          estimate_prior_method = "optim",
                          robust_estimation = F)
    }, error = function(err) {
      print(paste("vanilla: Error in the", tt, "trial: ", err))
      res_gs_vn <- NULL
    })

    # Huber weighting
    tryCatch({
      res_gs_hb_S <- gsusie(X, y, family = "poisson",
                            estimate_prior_method = "optim",
                            robust_estimation = T,
                            robust_method = "huber",
                            huber_tuning_k = "S")
    }, error = function(err) {
      print(paste("Huber-S: Error in the", tt, "trial: ", err))
      res_gs_hb_S <- NULL
    })

    tryCatch({
      res_gs_hb_M <- gsusie(X, y, family = "poisson",
                            estimate_prior_method = "optim",
                            robust_estimation = T,
                            robust_method = "huber",
                            huber_tuning_k = "M")
    }, error = function(err) {
      print(paste("Huber-M: Error in the", tt, "trial: ", err))
      res_gs_hb_M <- NULL
    })

    tryCatch({
      res_gs_hb_vSD <- gsusie(X, y, family = "poisson",
                              estimate_prior_method = "optim",
                              robust_estimation = T,
                              robust_method = "huber",
                              huber_tuning_k = "vSD")
    }, error = function(err) {
      print(paste("Huber-vSD: Error in the", tt, "trial: ", err))
      res_gs_hb_vSD <- NULL
    })

    tryCatch({
      res_gs_hb_vMd <- gsusie(X, y, family = "poisson",
                              estimate_prior_method = "optim",
                              robust_estimation = T,
                              robust_method = "huber",
                              huber_tuning_k = "vMd")
    }, error = function(err) {
      print(paste("Huber-vMd: Error in the", tt, "trial: ", err))
      res_gs_hb_vMd <- NULL
    })


    # Bisquare weighting
    tryCatch({
      res_gs_bs <- gsusie(X, y, family = "poisson",
                          estimate_prior_method = "optim",
                          robust_estimation = T,
                          robust_method = "bisquare")

    }, error = function(err) {
      print(paste("Bisquare-vMd: Error in the", tt, "trial: ", err))
      res_gs_bs <- NULL
    })

    # Weight dropped by fractions - 0.3
    tryCatch({
      res_gs_fc30 <- gsusie(X, y, family = "poisson",
                            estimate_prior_method = "optim",
                            robust_estimation = T,
                            robust_method = "simple",
                            simple_outlier_fraction = 0.3,
                            simple_outlier_thres = NULL)
    }, error = function(err) {
      print(paste("Simple-frac 0.3: Error in the", tt, "trial: ", err))
      res_gs_fc30 <- NULL
    })

    # Weight dropped by fractions - 0.05
    tryCatch({
      res_gs_fc5 <- gsusie(X, y, family = "poisson",
                           estimate_prior_method = "optim",
                           robust_estimation = T,
                           robust_method = "simple",
                           simple_outlier_fraction = 0.05,
                           simple_outlier_thres = NULL)
    }, error = function(err) {
      print(paste("Simple-frac 0.05: Error in the", tt, "trial: ", err))
      res_gs_fc5 <- NULL
    })


    # Weight dropped by fractions - 0.01
    tryCatch({
      res_gs_fc1 <- gsusie(X, y, family = "poisson",
                           estimate_prior_method = "optim",
                           robust_estimation = T,
                           robust_method = "simple",
                           simple_outlier_fraction = 0.01,
                           simple_outlier_thres = NULL)
    }, error = function(err) {
      print(paste("Simple-frac 0.01: Error in the", tt, "trial: ", err))
      res_gs_fc1 <- NULL
    })


    ### Comparison: GLMNET
    # Lasso
    res_la <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 1)
    # Ridge regression
    res_rg <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 0)
    # Elastic net: alpha=0.5
    res_en <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 0.5)


    #####################################
    ## convergeFlags
    convergeFlags$gs_vn[tt]  <- res_gs_vn$converged
    convergeFlags$gs_hb_S[tt]    <- res_gs_hb_S$converged
    convergeFlags$gs_hb_M[tt]    <- res_gs_hb_M$converged
    convergeFlags$gs_hb_vSD[tt]  <- res_gs_hb_vSD$converged
    convergeFlags$gs_hb_vMd[tt]  <- res_gs_hb_vMd$converged
    convergeFlags$gs_bs[tt]  <- res_gs_bs$converged
    convergeFlags$gs_fc30[tt] <- res_gs_fc30$converged
    convergeFlags$gs_fc5[tt] <- res_gs_fc5$converged
    convergeFlags$gs_fc1[tt] <- res_gs_fc1$converged

    ## ELBO
    ELBO$gs_vn[tt]  <- res_gs_vn$elbo[length(res_gs_vn$elbo)]
    ELBO$gs_hb_S[tt]    <- res_gs_hb_S$elbo[length(res_gs_hb_S$elbo)]
    ELBO$gs_hb_M[tt]    <- res_gs_hb_S$elbo[length(res_gs_hb_M$elbo)]
    ELBO$gs_hb_vSD[tt]  <- res_gs_hb_S$elbo[length(res_gs_hb_vSD$elbo)]
    ELBO$gs_hb_vMd[tt]  <- res_gs_hb_S$elbo[length(res_gs_hb_vMd$elbo)]
    ELBO$gs_bs[tt]  <- res_gs_bs$elbo[length(res_gs_bs$elbo)]
    ELBO$gs_fc30[tt] <- res_gs_fc30$elbo[length(res_gs_fc30$elbo)]
    ELBO$gs_fc5[tt] <- res_gs_fc5$elbo[length(res_gs_fc5$elbo)]
    ELBO$gs_fc1[tt] <- res_gs_fc1$elbo[length(res_gs_fc1$elbo)]


    #####################################
    ## Evaluations - variable selection
    vars_gs_vn   <- get_selected_vars(res_gs_vn)
    vars_gs_hb_S <- get_selected_vars(res_gs_hb_S)
    vars_gs_hb_M <- get_selected_vars(res_gs_hb_M)
    vars_gs_hb_vSD   <- get_selected_vars(res_gs_hb_vSD)
    vars_gs_hb_vMd   <- get_selected_vars(res_gs_hb_vMd)
    vars_gs_bs   <- get_selected_vars(res_gs_bs)
    vars_gs_fc30 <- get_selected_vars(res_gs_fc30)
    vars_gs_fc5  <- get_selected_vars(res_gs_fc5)
    vars_gs_fc1  <- get_selected_vars(res_gs_fc1)
    vars_la      <- get_selected_vars(res_la)
    vars_la10    <- get_selected_vars(res_la, glmnet_top_n = n_effect_vars)
    vars_rg      <- get_selected_vars(res_rg)
    vars_rg10    <- get_selected_vars(res_rg, glmnet_top_n = n_effect_vars)
    vars_en      <- get_selected_vars(res_en)
    vars_en10    <- get_selected_vars(res_en, glmnet_top_n = n_effect_vars)

    ### Number of variables selected
    NumVs$gs_vn[tt]   <- length(vars_gs_vn)
    NumVs$gs_hb_S[tt]   <- length(vars_gs_hb_S)
    NumVs$gs_hb_M[tt]   <- length(vars_gs_hb_M)
    NumVs$gs_hb_vSD[tt]   <- length(vars_gs_hb_vSD)
    NumVs$gs_hb_vMd[tt]   <- length(vars_gs_hb_vMd)
    NumVs$gs_bs[tt]   <- length(vars_gs_bs)
    NumVs$gs_fc30[tt] <- length(vars_gs_fc30)
    NumVs$gs_fc5[tt]  <- length(vars_gs_fc5)
    NumVs$gs_fc1[tt]  <- length(vars_gs_fc5)
    NumVs$la[tt]      <- length(vars_la)
    NumVs$la10[tt]    <- length(vars_la10) # 10
    NumVs$rg[tt]      <- length(vars_rg)
    NumVs$rg10[tt]    <- length(vars_rg10) # 10
    NumVs$en[tt]      <- length(vars_en)
    NumVs$en10[tt]    <- length(vars_en10) # 10

    ### Precision
    PPV$gs_vn[tt]   <- evaluate_ppv(vars_gs_vn,  effect_idx)
    PPV$gs_hb_S[tt]   <- evaluate_ppv(vars_gs_hb_S,  effect_idx)
    PPV$gs_hb_M[tt]   <- evaluate_ppv(vars_gs_hb_M,  effect_idx)
    PPV$gs_hb_vSD[tt]   <- evaluate_ppv(vars_gs_hb_vSD,  effect_idx)
    PPV$gs_hb_vMd[tt]   <- evaluate_ppv(vars_gs_hb_vMd,  effect_idx)
    PPV$gs_bs[tt]   <- evaluate_ppv(vars_gs_bs,  effect_idx)
    PPV$gs_fc30[tt] <- evaluate_ppv(vars_gs_fc30, effect_idx)
    PPV$gs_fc5[tt]  <- evaluate_ppv(vars_gs_fc5, effect_idx)
    PPV$gs_fc1[tt]  <- evaluate_ppv(vars_gs_fc5, effect_idx)
    PPV$la[tt]      <- evaluate_ppv(vars_la, effect_idx)
    PPV$la10[tt]    <- evaluate_ppv(vars_la10, effect_idx)
    PPV$rg[tt]      <- evaluate_ppv(vars_rg, effect_idx)
    PPV$rg10[tt]    <- evaluate_ppv(vars_rg10, effect_idx)
    PPV$en[tt]      <- evaluate_ppv(vars_en, effect_idx)
    PPV$en10[tt]    <- evaluate_ppv(vars_en10, effect_idx)

    ### Recall
    TPR$gs_vn[tt]  <- evaluate_tpr(vars_gs_vn,  effect_idx)
    TPR$gs_hb_S[tt]  <- evaluate_tpr(vars_gs_hb_S,  effect_idx)
    TPR$gs_hb_M[tt]  <- evaluate_tpr(vars_gs_hb_M,  effect_idx)
    TPR$gs_hb_vSD[tt]  <- evaluate_tpr(vars_gs_hb_vSD,  effect_idx)
    TPR$gs_hb_vMd[tt]  <- evaluate_tpr(vars_gs_hb_vMd,  effect_idx)
    TPR$gs_bs[tt]  <- evaluate_tpr(vars_gs_bs,  effect_idx)
    TPR$gs_fc30[tt]<- evaluate_tpr(vars_gs_fc30, effect_idx)
    TPR$gs_fc5[tt] <- evaluate_tpr(vars_gs_fc5, effect_idx)
    TPR$gs_fc1[tt] <- evaluate_tpr(vars_gs_fc5, effect_idx)
    TPR$la[tt]     <- evaluate_tpr(vars_la, effect_idx)
    TPR$la10[tt]   <- evaluate_tpr(vars_la10, effect_idx)
    TPR$rg[tt]     <- evaluate_tpr(vars_rg, effect_idx)
    TPR$rg10[tt]   <- evaluate_tpr(vars_rg10, effect_idx)
    TPR$en[tt]     <- evaluate_tpr(vars_en, effect_idx)
    TPR$en10[tt]   <- evaluate_tpr(vars_en10, effect_idx)


    ### AUPRC: remember to remove intercept!
    labs <- 1 * ((1 : pp) %in% effect_idx)  # drop the intercept term!
    AUPRC$gs_vn[tt] <- pr.curve(scores.class0 = res_gs_vn$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_hb_S[tt] <- pr.curve(scores.class0 = res_gs_hb_S$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_hb_M[tt] <- pr.curve(scores.class0 = res_gs_hb_M$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_hb_vSD[tt] <- pr.curve(scores.class0 = res_gs_hb_vSD$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_hb_vMd[tt] <- pr.curve(scores.class0 = res_gs_hb_vMd$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_bs[tt] <- pr.curve(scores.class0 = res_gs_bs$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_fc30[tt] <- pr.curve(scores.class0 = res_gs_fc30$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_fc5[tt] <- pr.curve(scores.class0 = res_gs_fc5$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    AUPRC$gs_fc1[tt] <- pr.curve(scores.class0 = res_gs_fc1$pip[-ncol(X)],
                                weights.class0 = labs)$auc.integral
    ## remember to remove the intercept in the fitted GLMNET result!
    la_coefs <- as.numeric(coef(res_la, s = "lambda.min"))[-1]
    AUPRC$la[tt] <- pr.curve(scores.class0 = abs(la_coefs),
                             weights.class0 = labs)$auc.integral
    rg_coefs <- as.numeric(coef(res_rg, s = "lambda.min"))[-1]
    AUPRC$rg[tt] <- pr.curve(scores.class0 = abs(rg_coefs),
                             weights.class0 = labs)$auc.integral
    en_coefs <- as.numeric(coef(res_en), s = "lambda.min")[-1]
    AUPRC$en[tt] <- pr.curve(scores.class0 = abs(en_coefs),
                             weights.class0 = labs)$auc.integral

    ###############
    # remove
    rm_ls <- ls()[startsWith(ls(), "vars_") | startsWith(ls(), "res_")]
    rm_ls <- append(rm_ls, ls()[endsWith(ls(), "_coefs")])
    rm_ls <- append(rm_ls, "labs")
    rm(list = rm_ls)

    ###############
    ## for clusters: save results every 50 runs
    if (tt %% 10 == 0) {  #
      saveRDS(list(Precision = PPV,
                   Recall    = TPR,
                   AUPRC     = AUPRC,
                   NumVs     = NumVs,
                   convergeFlags = convergeFlags,
                   ELBO = ELBO),
              file = .result.dir %&% .filenames)
    }

  }

  ###############
  ## for local computer: save results at the end
  saveRDS(list(Precision = PPV,
               Recall    = TPR,
               AUPRC     = AUPRC,
               NumVs     = NumVs,
               convergeFlags = convergeFlags,
               ELBO = ELBO),
          file = .result.dir %&% .filenames)
})
