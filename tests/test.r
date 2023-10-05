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

filenames <- filenames[filenames != "gsusie.r"]  # remove the "test.r"

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
h2 <- 0.8

X <- matrix(rnorm(nn * pp), ncol = pp)
Eta <- sqrt(h2) * (-2*X[,1] + 0.5*X[,6]) + sqrt(1-h2) * rnorm(nn)
expEta <- exp(Eta)
y1 <- rpois(nn, expEta)
# y <- y / sum(y) * 10
y <- exp(scale(log1p(y1)))
plot(y)

# res_susie <- susie(X, y, L = 10)
# summary(res_susie)s

## our method
X <- cbind(X, 1)
colnames(X) <- paste0("X", c(1:pp, 0))

dim(X)
length(y)


res_gsusie <- gsusie(X, y, family = "poisson", maxL = 10,
                     max_iters = 500, tol = 1e-2,
                     coef_prior_variance = 1,
                     estimate_prior_method = "optim",
                     robust_estimation = T,
                     robust_method = "simple",
                     simple_outlier_fraction = NULL,
                     simple_outlier_thres = 1000,
                     verbose = T)

# res_gsusie <- gsusie(X, y, family = "poisson", maxL = 10,
#                      max_iters = 500, tol = 1e-2,
#                      coef_prior_variance = 1,
#                      estimate_prior_method = "optim",
#                      robust_estimation = T,
#                      robust_method = "huber",
#                      verbose = T)

res_gsusie <- gsusie(X, y, family = "poisson", maxL = 10,
                     max_iters = 500, tol = 1e-2,
                     coef_prior_variance = 1,
                     estimate_prior_method = "optim",
                     robust_estimation = F,
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
