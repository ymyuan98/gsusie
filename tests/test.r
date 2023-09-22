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
p <- 2000
h2 <- 0.6
X <- matrix(rnorm(n * p), ncol = p)
# b <- 2
# b <- 0.8
# cat("Effect size of X1:", b, "\n")  # Effect size of X1: 1.683509

# Eta <- sqrt(h2) * (2 * X[, 1] - 0.8 * X[, 6]) + sqrt(1 - h2) * rnorm(n)
Eta <- sqrt(h2) * (2 * X[, 1]) + sqrt(1 - h2) * rnorm(n)
y <- rbinom(n, 1, prob = expit(Eta))


# res_gs <- gsusie(cbind(1, X), y, family = "binomial")  # default: grid_opt_prior_variance = F

res_gs <- gsusie(cbind(1, X), y, family = "binomial", maxL = 10,
                 coef_prior_variance = 1)  # default: grid_opt_prior_variance = F

# res_gs <- gsusie(X, y, family = "binomial", 
#                   grid_opt_prior_variance = T, 
#                   grid_prior_variance_value = c(0.1, 0.5, 1))


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
files_path <- "./gsusie-package/R/"
filenames <- list.files(path = files_path)

filenames <- filenames[filenames != "test.r"]  # remove the "test.r"

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
X <- matrix(rnorm(nn * pp), ncol = pp)

# Eta <- -2 * X[, 1, drop = F]
Eta <- -2 * X[,1] + 0.5 * X[,6]
expEta <- exp(Eta)
y1 <- rpois(nn, expEta) 
# y <- y / sum(y) * 10 
y <- exp(scale(log1p(y1)))

# res_susie <- susie(X, y, L = 10)
# summary(res_susie)


## our method
X <- cbind(1, X)
colnames(X) <- paste0("X", 0:pp)
res_gsusie <- gsusie(X, y, family = "poisson", maxL = 5,
                     max_iters = 500, tol = 1e-2, 
                     coef_prior_variance = 1,  
                     verbose = T)
# res_gsusie$sets
# round(res_gsusie$alpha, digits = 3)
plot(res_gsusie$pip)

summary.gsusie(res_gsusie)
coefficients.gsusie(res_gsusie)
res_gsusie$elbo[res_gsusie$niter]

indices <- 1:res_gsusie$niter
indices <- 40:80
plot(indices, res_gsusie$elbo[indices])
plot(indices, res_gsusie$loglik_exact[indices])
plot(indices, res_gsusie$elbo[indices] - res_gsusie$loglik_exact[indices])
# res_gsusie$elbo[indices]
# which.max(res_gsusie$elbo)

# bb <- res_gsusie$elbo[150 : 1000] - res_gsusie$elbo[149 : 999]
# plot(1:length(bb), bb)


res_susie <- susie(X, y)
summary(res_susie)
coefficients.gsusie(res_susie)
##############################################################################


nn <- 100
pp <- 20
XX <- matrix(rnorm(nn * pp), ncol = pp)
# yy <- rbinom(nn, 1, expit(XX[,1]))
yy <- rpois(nn, exp(XX[,1]))
res1 <- glm(yy ~ XX[, 2] + XX[, 3], family = "poisson")
summary(res1)
res2 <- glm(yy ~ -1 + XX[, 9] , family = "poisson")
summary(res2)
