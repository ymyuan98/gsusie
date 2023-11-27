set.seed(Sys.Date())

require(tidyverse)

expit <- function(x) {
  res <- ifelse(x > 0,
                1 / (1 + exp(-x)),      ## if yes
                exp(x) / (exp(x) + 1))  ## if no
  return(res)
}

## Synthesize data
n <- 500
p <- 10
h2 <- 0.6
X <- matrix(rnorm(n * p), ncol = p)
# b <- 2
# b <- 0.8
# cat("Effect size of X1:", b, "\n")  # Effect size of X1: 1.683509

# Eta <- sqrt(h2) * (2 * X[, 1] - 0.8 * X[, 6]) + sqrt(1 - h2) * rnorm(n)
Eta <- sqrt(h2) * (2 * X[, 1]) + sqrt(1 - h2) * rnorm(n)
y <- rbinom(n, 1, prob = expit(Eta))

XX <- cbind(1, X)
colnames(XX) <- paste0("X", 0:p)
res_gs <- gsusie(XX, y, family = "binomial", maxL = 5,
                 coef_prior_variance = 1)


res_gs$niter
res_gs$sets
summary.gsusie(res_gs)
coefficients.gsusie(res_gs)

round(res_gs$alpha, digits = 3)

plot(1:res_gs$niter, res_gs$elbo)
plot(1:res_gs$niter, res_gs$loglik_exact)
plot(1:res_gs$niter, res_gs$loglik_apprx)



library(glmnet)
library(coefplot)
res_glmnet <- cv.glmnet(X, y, family = "binomial", alpha = 1)
summary(res_glmnet)
coef(res_glmnet)
extract.coef(res_glmnet)  # No SE provided when n > p.

res_glm <- glm(y ~ X, family = binomial)
summary(res_glm)


library(susieR)

res_s <- susie(X, y)
summary(res_s)
coefficients.gsusie(res_s)


################################################################################
########################################
`%&%` <- function(a, b) paste0(a, b)

# getwd()
files_path <- "./R/"
filenames <- list.files(path = files_path)

# filenames <- filenames[filenames != "gsusie.r"]  # remove the "test.r"

for (i in 1 : length(filenames)) {
  source(files_path %&% filenames[i])
}


require(susieR)
########################################
########################################



##########################
set.seed(20230718)
library(MASS)

## Synthesize data
nn <- 500
pp <- 1000
h2 <- 0.5

X <- matrix(rnorm(nn * pp), ncol = pp)
n_effect_vars <- 10  ## number of effective variables
effect_idx <- sort(sample(1 : pp, size = n_effect_vars))
bb <- rnorm(n_effect_vars)
data.frame(variable = paste0("X", effect_idx), effect_size = bb)

Eta <- sqrt(h2) * (X[ ,effect_idx, drop = F] %*% as.matrix(bb)) +
  sqrt(1-h2) * rnorm(nn)
expEta <- exp(Eta)
y1 <- rpois(nn, expEta)
y <- exp(scale(log1p(y1)))
plot(y)



## our method
X <- cbind(X, 1)
colnames(X) <- paste0("X", c(1:pp, 0))

dim(X)
length(y)

## Non-robust-estimation
res_gs_vn <- gsusie(X, y, family = "poisson", maxL = 10,
                     max_iters = 100, tol = 1e-2,
                     coef_prior_variance = 1,
                     estimate_prior_method = "optim",
                     robust_estimation = F,
                     verbose = T)

# ## Threshold-removal
# res_gs_thres <- gsusie(X, y, family = "poisson", maxL = 10,
#                      max_iters = 500, tol = 1e-2,
#                      coef_prior_variance = 1,
#                      estimate_prior_method = "optim",
#                      robust_estimation = T,
#                      robust_method = "simple",
#                      simple_outlier_fraction = 0.01,
#                      simple_outlier_thres = NULL,
#                      verbose = T)

## Huber-weighting
res_gs_hb_M <- gsusie(X, y, family = "poisson", maxL = 10,
                     max_iters = 100, tol = 1e-2,
                     coef_prior_variance = 1,
                     estimate_prior_method = "optim",
                     robust_estimation = T,
                     robust_method = "huber",
                     tuning_k = "M")
res_gs_hb_S <- gsusie(X, y, family = "poisson", maxL = 10,
                       max_iters = 100, tol = 1e-2,
                       coef_prior_variance = 1,
                       estimate_prior_method = "optim",
                       robust_estimation = T,
                       robust_method = "huber",
                       tuning_k = "S")

## bisquare-weighting
res_gs_bs_M <- gsusie(X, y, family = "poisson", maxL = 10,
                       max_iters = 500, tol = 1e-2,
                       coef_prior_variance = 1,
                       estimate_prior_method = "optim",
                       robust_estimation = T,
                       robust_method = "bisquare",
                      tuning_k = "M",
                       verbose = T)

res_gs_bs_S <- gsusie(X, y, family = "poisson", maxL = 15,
                      max_iters = 500, tol = 1e-2,
                      coef_prior_variance = 1,
                      estimate_prior_method = "optim",
                      robust_estimation = T,
                      robust_method = "bisquare",
                      tuning_k = "S",
                      verbose = T)



# res_gsusie$sets
# round(res_gsusie$alpha, digits = 3)
plot(res_gsusie$pip, ylab = "PIP")
summary.gsusie(res_gsusie)
coefficients.gsusie(res_gsusie)
res_gsusie$elbo[res_gsusie$niter]
res_gsusie$V


indices <- 1:res_gsusie$niter
indices <- 300:500
plot(indices, res_gsusie$elbo[indices])
plot(indices, res_gsusie$loglik_exact[indices])
plot(indices, res_gsusie$elbo[indices] - res_gsusie$loglik_exact[indices])


res_susie <- susie(X, y)
summary(res_susie)
coefficients.gsusie(res_susie)
##############################################################################


nn <- 100
pp <- 20
XX <- matrix(rnorm(nn * pp), ncol = pp)
yy <- rpois(nn, exp(XX[,1]))
res1 <- glm(yy ~ XX[, 2] + XX[, 3], family = "poisson")
summary(res1)
res2 <- glm(yy ~ -1 + XX[, 9] , family = "poisson")
summary(res2)


##############################################################################
.result.dir <- "./tests/"
.filenames <- "poisson-n500-p1000-h5.rds"
RES <- readRDS(.result.dir %&% .filenames)

PPV <- RES$Precision
TPR <- RES$Recall
NumVs <- RES$NumVs
convergeFlags <- RES$convergeFlags
ELBO <- RES$ELBO



################################################################################
library(tidyverse)

set.seed(20231020)

h2 <- 0.5
## Generate data
nn <- 500
pp <- 1000
X <- matrix(rnorm(nn * pp), ncol = pp)
n_effect_vars <- 5  ## number of effective variables
effect_idx <- sort(sample(1 : pp, size = n_effect_vars))
bb <- rnorm(n_effect_vars)
# data.frame(variable = paste0("X", effect_idx), effect_size = bb)

Eta <- sqrt(h2) * (X[ ,effect_idx, drop = F] %*% as.matrix(bb)) +
  sqrt(1-h2) * rnorm(nn)
y1 <- rpois(nn, exp(Eta))
y  <- exp(scale(log1p(y1)))  # scale the response


X <- cbind(X, 1)  # Intercept is always put at last
colnames(X) <- paste0("X", c(1:pp, 0))

plot(y)

###############################################################################

res_gs_vn <- gsusie(X, y, family = "poisson",
                    robust_estimation = F)
summary(res_gs_vn)
gsusie_coefficients(res_gs_vn)

res_gs_hbS <- gsusie(X, y, family = "poisson",
                     robust_estimation = T,
                     robust_method = "huber",
                     huber_tuning_k = "S")
summary(res_gs_hbS)
gsusie_coefficients(res_gs_hbS)

res_gs_hbM <- gsusie(X, y, family = "poisson",
                     robust_estimation = T,
                     robust_method = "huber",
                     huber_tuning_k = "M")
summary(res_gs_hbM)
gsusie_coefficients(res_gs_hbM)


################################################################################













