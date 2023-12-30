
## source functions
.libPaths("~/projects/Rlibs/")
library(glmnet)
# library(coefplot)
library(bigsnpr)
library(susieR)

`%&%` <- function(a, b) paste0(a, b)
if.needed <- function(.files, .code) {
  if(!all(file.exists(.files))){
    .code
  }
}

# getwd()
# files.dir <- "./R/"
files.dir <- "~/projects/gsusie/R/"
r_file_names <- list.files(path = files.dir)
for (i in 1 : length(r_file_names)) {
  source(files.dir %&% r_file_names[i])
}

# data path
# data.dir <- "~/Documents/research/data/genotype/"  ## local
## on clusters
data.dir <- "~/projects/data/genotype/"
# .subset.filename <- function(.index) {
#   return(paste0("snp_subset_", sprintf('%04d', .index), ".rds"))
# }

################

# generate filename that store the results of simulations.
# .prefix can be used to add date/"sim"/"snp"
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

model <- "binomial"
.result.dir <- "~/projects/gsusie/results/" %&%
  Sys.Date() %&% "-genotype-" %&% model %&% "/"
# .result.dir <- "./results/genotype-" %&% Sys.Date() %&% "-" %&% model %&% "/"
dir.create(.result.dir, recursive=TRUE, showWarnings=FALSE)

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
                     "su", "la", "rr", "en")


  expit <- function(eta) {
    res <- ifelse(eta > 0,
                  1 / (1 + exp(-eta)),
                  exp(eta) / (1 + exp(eta)))
    return(res)
  }
  #################################
  ## Data generation

  # preload_bed_data
  .bk.file <- data.dir %&% "1000G_phase3_common_norel.rds"
  if(!file.exists(.bk.file)){
    BED <- snp_readBed(.bed.file)
  }
  dat <- snp_attach(.bk.file)$genotypes

  set.seed(20231102)
  startpoint <- sample(1 : (ncol(dat)-10000), size = 9000)[seed]
  # size = the number of total trials we run

  set.seed(startpoint)
  ## sub-sample individuals
  if (nn < 2490) { ii.idx <- sample(1 : 2490, size = nn) }
  X <- dat[ii.idx, startpoint : (startpoint+pp-1)]
  rm(dat)

  effect_idx <- sort(sample(1 : pp, size = n_effect_vars))
  bb <- rnorm(n_effect_vars)
  # data.frame(variable = paste0("X", effect_idx), effect_size = bb)

  Eta <- sqrt(h2) * (scale(X[ ,effect_idx, drop = F]) %*% as.matrix(bb)) +
    sqrt(1-h2) * rnorm(nn)
  y <- rbinom(nn, 1, expit(Eta))  # binomial
  X <- cbind(X, 1)  # Intercept is always put at last
  colnames(X) <- paste0("X", c(1:pp, 0))

  #################################
  ## Method comparison

  # Non-robust estimation (vanilla)
  res_gs_vn <- tryCatch({
    gsusie(X, y, family = "binomial",
           estimate_prior_method = "optim",
           robust_estimation = F
           )
  }, error = function(err) {
    message("Vanilla. seed:", seed, ". error: ", err)
  })

  # # Huber weighting
  # res_gs_hb_S <- tryCatch({
  #   gsusie(X, y, family = "binomial",
  #          estimate_prior_method = "optim",
  #          robust_estimation = T,
  #          robust_method = "huber",
  #          robust_tuning_method = "S")
  # }, error = function(err) {
  #   message("Huber-S. seed: ", seed, " error: ", err)
  # })
  #
  #
  # res_gs_hb_M <- tryCatch({
  #   gsusie(X, y, family = "binomial",
  #          estimate_prior_method = "optim",
  #          robust_estimation = T,
  #          robust_method = "huber",
  #          robust_tuning_method = "M")
  # }, error = function(err) {
  #   message("Huber-M. seed: ", seed, " error: ", err)
  # })
  #
  # # Bisquare weighting
  # res_gs_bs_S <- tryCatch({
  #   gsusie(X, y, family = "binomial",
  #          estimate_prior_method = "optim",
  #          robust_estimation = T,
  #          robust_method = "bisquare",
  #          robust_tuning_method = "S")
  # }, error = function(err) {
  #   message("Bisquare-S. seed: ", seed, " error: ", err)
  # })
  #
  # res_gs_bs_M <- tryCatch({
  #   gsusie(X, y, family = "binomial",
  #          estimate_prior_method = "optim",
  #          robust_estimation = T,
  #          robust_method = "bisquare",
  #          robust_tuning_method = "M")
  # }, error = function(err) {
  #   message("Bisquare-M. seed: ", seed, " error: ", err)
  # })

  # Weight dropped by fractions - 0.3
  res_gs_fc30 <- tryCatch({
    gsusie(X, y, family = "binomial",
           estimate_prior_method = "optim",
           robust_estimation = T,
           robust_method = "simple",
           simple_outlier_fraction = 0.3,
           simple_outlier_thres = NULL)
  }, error = function(err) {
    message("Simple-frac-0.3. seed: ", seed, " error: ", err)
  })

  # Weight dropped by fractions - 0.05
  res_gs_fc5 <- tryCatch({
    gsusie(X, y, family = "binomial",
           estimate_prior_method = "optim",
           robust_estimation = T,
           robust_method = "simple",
           simple_outlier_fraction = 0.05,
           simple_outlier_thres = NULL)
  }, error = function(err) {
    message("Simple-frac-0.05. seed: ", seed, " error: ", err)
  })

  # Weight dropped by fractions - 0.01
  res_gs_fc1 <- tryCatch({
    gsusie(X, y, family = "binomial",
           estimate_prior_method = "optim",
           robust_estimation = T,
           robust_method = "simple",
           simple_outlier_fraction = 0.01,
           simple_outlier_thres = NULL)
  }, error = function(err) {
    message("Simple-frac-0.01. seed: ", seed, " error: ", err)
  })


  ### Comparison: SuSiE
  res_su <- susie(X = X[, -ncol(X)], y = y)

  ### Comparison: GLMNET
  # Lasso
  res_la <- cv.glmnet(x = X[, -ncol(X)], y, family = "binomial", alpha = 1)
  # Ridge regression
  # res_rr <- cv.glmnet(x = X[, -ncol(X)], y, family = "binomial", alpha = 0)
  # Elastic net: alpha=0.5
  res_en <- cv.glmnet(x = X[, -ncol(X)], y, family = "binomial", alpha = 0.5)


  ##############################################################################
  extract_pip <- function(gs_res) {
    if (is.null(gs_res)) {
      # gs is not fitted due to certain results
      pip <- rep(NA, times = pp)
    }
    else if (class(gs_res) %in% "gsusie"){
      pip <- gs_res$pip
      names(pip) <- names(gs_res$pip)
      pip <- pip[-length(pip)]  # remove the PIP of X0 (pip[-length(pip)])
    }
    return(pip)
  }

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
  output$pip <- data.frame(
    cbind(extract_pip(res_gs_vn), # move intercept to the front
          extract_pip(res_gs_hb_M),
          extract_pip(res_gs_hb_S),
          extract_pip(res_gs_bs_M),
          extract_pip(res_gs_bs_S),
          extract_pip(res_gs_fc1),
          extract_pip(res_gs_fc5),
          extract_pip(res_gs_fc30),
          res_su$pip),
    row.names = paste0("X", 1:pp)
  )
    colnames(output$pip) <-
      # c(methods_names[startsWith(methods_names, "gs")], "su")
      c("gs_vn", "gs_fc1", "gs_fc5", "gs_fc30", "su")

  output$gn_coef <- data.frame(
    cbind(as.numeric(coef(res_la, s = "lambda.min"))[-1],
          # as.numeric(coef(res_rr, s = "lambda.min"))[-1],
          as.numeric(coef(res_en, s = "lambda.min"))[-1]),
    row.names = paste0("X", 1:pp)
  )
  colnames(output$gn_coef) <- c("la", "en")

  saveRDS(output, file = .result.dir %&% .filename)
})
