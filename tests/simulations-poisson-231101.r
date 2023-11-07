
## source functions
.libPaths("~/projects/Rlibs/")
library(glmnet)
# library(readr)
# library(PRROC)

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

# generate filename that store the results of simulations.
# .prefix can be used to add date/"sim"/"real"
gen_filename <- function(n, p, h2, mod, seed,
                          .prefix = NULL) {
  if (h2 < 0.1 & h2 >= 0.01) {
    .filename <- mod %&% "-n" %&% n %&% "-p" %&% p %&% "-hz" %&% (100*h2) %&%
      "-s" %&% seed %&% ".rds"
  } else {
    .filename <- mod %&% "-n" %&% n %&% "-p" %&% p %&% "-h" %&% (10*h2) %&%
      "-s" %&% seed %&% ".rds"
  }

  if (!is.null(.prefix)) {
    .filename <- .prefix %&% "-" %&% .filename
  }
  return(.filename)
}

################
## simulations

## read in arguments
args <- commandArgs(trailingOnly = TRUE)
ID <- as.numeric(args[1])
seed <- as.numeric(args[2])
nn <- as.numeric(args[3])
pp <- as.numeric(args[4])
h2 <- as.numeric(args[5])
n_effect_vars <- as.numeric(args[6])  ## number of effective variables
# nn <- 500
# pp <- 1000
# h2 <- 0.5

model <- "poisson"
.result.dir <- "~/projects/gsusie/results/" %&%
  Sys.Date() %&% "-synthetic-" %&% model %&% "/"
# .result.dir <- "./results/synthetic-" %&% Sys.Date() %&% "-" %&% model %&% "/"
if (!dir.exists(.result.dir)) {
  dir.create(.result.dir, recursive=TRUE, showWarnings=FALSE)
}

.filename <- gen_filename(nn, pp, h2, model, seed, "ID"%&%ID)

if.needed(.result.dir %&% .filename, {

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
  methods_names <- c(paste("gs", c("vn", "hb_S", "hb_M", "bs_S", "bs_M",
                                   "fc30", "fc5", "fc1"),
                           sep = "_"),
                     c("la", "rr", "en"))


  #################################
  ## Data generation
  set.seed(seed)
  X <- matrix(rnorm(nn * pp), ncol = pp)
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
    print(paste("Vanilla. seed:", seed, "error:", err))
    res_gs_vn <- list()
    res_gs_vn$pip <- rep(NA, times = ncol(X))
  })

  # Huber weighting
  tryCatch({
    res_gs_hb_S <- gsusie(X, y, family = "poisson",
                          estimate_prior_method = "optim",
                          robust_estimation = T,
                          robust_method = "huber",
                          tuning_k = "S")
  }, error = function(err) {
    print(paste("Huber-S. seed:", seed, "error:", err))
    res_gs_hb_S <- list()
    res_gs_hb_S$pip <- rep(NA, times = ncol(X))
  })

  tryCatch({
    res_gs_hb_M <- gsusie(X, y, family = "poisson",
                          estimate_prior_method = "optim",
                          robust_estimation = T,
                          robust_method = "huber",
                          tuning_k = "M")
  }, error = function(err) {
    print(paste("Huber-M. seed:", seed, "error:", err))
    res_gs_hb_M <- list()
    res_gs_hb_M$pip <- rep(NA, times = ncol(X))
  })

  # Bisquare weighting
  tryCatch({
    res_gs_bs_S <- gsusie(X, y, family = "poisson",
                            estimate_prior_method = "optim",
                            robust_estimation = T,
                            robust_method = "bisquare",
                            tuning_k = "S")
  }, error = function(err) {
    print(paste("Bisquare-S. seed:", seed, "error:", err))
    res_gs_bs_S <- list()
    res_gs_bs_S$pip <- rep(NA, times = ncol(X))
  })

  tryCatch({
    res_gs_bs_M <- gsusie(X, y, family = "poisson",
                            estimate_prior_method = "optim",
                            robust_estimation = T,
                            robust_method = "bisquare",
                            tuning_k = "M")
  }, error = function(err) {
    print(paste("Bisquare-M. seed:", seed, "error:", err))
    res_gs_bs_M <- list()
    res_gs_bs_M$pip <- rep(NA, times = ncol(X))
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
    print(paste("Simple-frac-0.3. seed:", seed, "error:", err))
    res_gs_fc30 <- list()
    res_gs_fc30$pip <- rep(NA, times = ncol(X))
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
    print(paste("Simple-frac-0.05. seed:", seed, "error:", err))
    res_gs_fc5 <- list()
    res_gs_fc5$pip <- rep(NA, times = ncol(X))
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
    print(paste("Simple-frac-0.01. seed:", seed, "error:", err))
    res_gs_fc1 <- list()
    res_gs_fc1$pip <- rep(NA, times = ncol(X))
  })

  ### Comparison: GLMNET
  # Lasso
  res_la <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 1)
  # Ridge regression
  res_rr <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 0)
  # Elastic net: alpha=0.5
  res_en <- cv.glmnet(x = X[, -(pp+1)], y, family = "poisson", alpha = 0.5)

  #############################################
  output <- list()

  ## save setup
  output$setup <- list(
    seed = seed,
    n = nn,
    p = pp,
    n_effect_vars = n_effect_vars,
    effect_idx = effect_idx,
    h2 = h2
  )

  ## save results
  output$gs_pip <- data.frame(
    cbind(res_gs_vn$pip[c((pp+1), 1:pp)], # move intercept to the front
          res_gs_hb_M$pip[c((pp+1), 1:pp)],
          res_gs_hb_S$pip[c((pp+1), 1:pp)],
          res_gs_bs_M$pip[c((pp+1), 1:pp)],
          res_gs_bs_S$pip[c((pp+1), 1:pp)],
          res_gs_fc1$pip[c((pp+1), 1:pp)],
          res_gs_fc5$pip[c((pp+1), 1:pp)],
          res_gs_fc30$pip[c((pp+1), 1:pp)]),
    row.names = paste0("X", 0:pp)
  )
  colnames(output$gs_pip) <- methods_names[startsWith(methods_names, "gs")]

  output$gn_coef <- data.frame(
    cbind(as.numeric(coef(res_la, s = "lambda.min")),
          as.numeric(coef(res_rr, s = "lambda.min")),
          as.numeric(coef(res_en, s = "lambda.min"))),
    row.names = paste0("X", 0:pp)
  )
  colnames(output$gn_coef) <- c("la", "rr", "en")

  saveRDS(output, file = .result.dir %&% .filename)
})
