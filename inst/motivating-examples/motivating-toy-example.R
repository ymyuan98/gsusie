library(tidyverse)
library(bigsnpr)
library(susieR)
library(glmnet)
library(MASS)
library(stringr)
library(patchwork)
library(gsusie)


`%&%` <- function(a, b) paste0(a,b)

source.dir <- "./R/"
r_file_names <- list.files(path = source.dir)
for (i in 1 : length(r_file_names)) {
  source(source.dir %&% r_file_names[i])
}


myPalette2 <- c("#E69F00", "#999999")
gsusie_pip_plot <- function(gs_res_pip, effect_idx) {
  tibble(pip = gs_res_pip) %>%
    mutate(id = 1 : n()) %>%
    mutate(true_label = if_else(id %in% effect_idx, "effect", "noneffect")) %>%
    arrange(desc(true_label)) %>%
    ggplot(aes(x = id, y = pip, color = as.factor(true_label))) +
    geom_point(size = 3) +
    scale_color_manual(values = myPalette2) +
    labs(x = "variable", color = "") +
    theme_classic()
}

###############################
## Toy example 1: independent predictors

set.seed(20231130)
## Synthesize data
nn <- 1000
pp <- 10
X <- matrix(rnorm(nn * pp), ncol = pp)

X[,1:2] <- MASS::mvrnorm(nn, mu = c(0,0),
                         Sigma = matrix(c(1, 0.8, 0.8, 1), nrow = 2))

X[,5:7] <- MASS::mvrnorm(nn, mu = rep(0, 3),
                         Sigma = matrix(c(1, 0.6, 0.9,
                                          0.6, 1, 0.75,
                                          0.9, 0.75, 1),
                                        nrow = 3, byrow = T))

effect_idx <- c(1, 6)
Eta <- scale(X[,effect_idx, drop=F] %*% as.matrix(c(-2, 0.2)))
y <- rpois(nn, exp(Eta))
plot(y)
plot(exp(scale(log1p(y))))

## G-SuSiE
res_gs_vn <- gsusie(cbind(X, 1), exp(scale(log1p(y))), family = "poisson")
summary(res_gs_vn)
gsusie_coefficients(res_gs_vn)

# res_gs_vn2 <- gsusie(cbind(X, 1), y, family = "poisson")
# summary(res_gs_vn2)
# gsusie_coefficients(res_gs_vn2)

res_gs_hb_M <- gsusie(cbind(X, 1), exp(scale(log1p(y))), family = "poisson",
                      robust_estimation = T,
                      robust_method = "huber",
                      robust_tuning_method = "M")
summary(res_gs_hb_M)
gsusie_coefficients(res_gs_hb_M)

# res_gs_hb_M2 <- gsusie(cbind(X, 1), y, family = "poisson",
#                        robust_estimation = T,
#                        robust_method = "huber",
#                        robust_tuning_method = "M")
# summary(res_gs_hb_M2)
# gsusie_coefficients(res_gs_hb_M2)

## SuSiE
res_su <- susie(X, exp(scale(log1p(y))))
summary(res_su)
gsusie_coefficients(res_su)

## Lasso
res_la <- glmnet(X, exp(scale(log1p(y))), family = "poisson", lambda = 1)
coef_la <- data.frame(coef = as.numeric(coefficients(res_la))[-1])
row.names(coef_la) <- paste0("X", 1 : pp)

## Elastic-Net
res_en <- glmnet(X, exp(scale(log1p(y))), family = "poisson", lambda = 0.5)
coef_en <- data.frame(coef = as.numeric(coefficients(res_en))[-1])
row.names(coef_en) <- paste0("X", 1 : pp)


p1 <- gsusie_pip_plot(res_gs_vn$pip[-(pp+1)], effect_idx) +
  labs(title = "GSuSiE (vanilla)", y = "PIP") +
  scale_x_continuous(name = "variable", breaks = seq(1, 10, by = 1)) +
  theme_classic()
p1
ggsave(filename = "independent-gsusie-vn.pdf",
       plot = p1, width = 9, height = 6,
       path = "./vignettes/motivating-example-independent/")

p2 <- gsusie_pip_plot(res_gs_hb_M$pip[-(pp+1)], effect_idx) +
  labs(title = "GSuSiE (Huber-M)", y = "PIP") +
  scale_x_continuous(name = "variable", breaks = seq(1, 10, by = 1)) +
  theme_classic()
p2

ggsave(filename = "independent-gsusie-hb-M.pdf",
       plot = p2, width = 9, height = 6,
       path = "./vignettes/motivating-example-independent/")

p3 <- gsusie_pip_plot(res_su$pip, effect_idx) +
  labs(title = "SuSiE", y = "PIP") +
  scale_x_continuous(name = "variable", breaks = seq(1, 10, by = 1))
p3
ggsave(filename = "independent-susie.pdf",
       plot = p3, width = 9, height = 6,
       path = "./vignettes/motivating-example-independent/")

p4 <- gsusie_pip_plot(abs(coef_la$coef), effect_idx) +
  labs(title = "Lasso", y = "|coefficient|") +
  scale_x_continuous(name = "variable", breaks = seq(1, 10, by = 1))
p4
ggsave(filename = "independent-lasso.pdf",
       plot = p4, width = 9, height = 6,
       path = "./vignettes/motivating-example-independent/")

p5 <- gsusie_pip_plot(abs(coef_en$coef), effect_idx) +
  labs(title = "Elastic Net", y = "|coefficient|") +
  scale_x_continuous(name = "variable", breaks = seq(1, 10, by = 1))
p5
ggsave(filename = "independent-elastic-net.pdf",
       plot = p5, width = 9, height = 6,
       path = "./vignettes/motivating-example-independent/")

wrap_plots(p2, p3, p4, p5, nrow = 2)

################################################################################
## Toy example 2: genotype data with n << p

# preload_bed_data
data.dir <- "~/Documents/research/data/genotype/"

.bed.file <- data.dir %&% "1000G_phase3_common_norel.bed"
if (!file.exists(.bed.file)) {
  bigsnpr::download_1000G(data.dir)
}

.bk.file <- data.dir %&% "1000G_phase3_common_norel.rds"
if (!file.exists(.bk.file)){
  BED <- snp_readBed(.bed.file)
}
dat <- snp_attach(.bk.file)$genotypes

nn <- 500
pp <- 2000
h2 <- 0.4
n_effect_vars <- 5

set.seed(20231129)
startpoint <- sample(1 : (ncol(dat)-10000), size = 9000)[7530]

set.seed(startpoint)
## sub-sample individuals
if (nn < 2490) { ii.idx <- sample(1 : 2490, size = nn) }
X <- dat[ii.idx, startpoint : (startpoint+pp-1)]

effect_idx <- sort(sample(1 : pp, size = n_effect_vars))
bb <- rnorm(n_effect_vars)
data.frame(variable = paste0("X", effect_idx), effect_size = bb)

Eta <- sqrt(h2) * (scale(X[ ,effect_idx, drop = F]) %*% as.matrix(bb)) +
  sqrt(1-h2) * rnorm(nn)
y1 <- rpois(nn, exp(Eta))
y  <- exp(scale(log1p(y1)))  # scale the response
plot(y)
X <- cbind(X, 1)  # Intercept is always put at last
colnames(X) <- paste0("X", c(1:pp, 0))

#################################
## Method comparison

# Non-robust estimation (vanilla)
res_gs_vn <- tryCatch({
  gsusie(X, y, family = "poisson",
         estimate_prior_method = "optim",
         robust_estimation = F)
}, error = function(err) {
  message("Vanilla. seed:", seed, ". error: ", err)
})
summary(res_gs_vn)
gsusie_coefficients(res_gs_vn)

### Huber-M
res_gs_hb_M <- tryCatch({
  gsusie(X, y, family = "poisson",
         estimate_prior_method = "optim",
         robust_estimation = T,
         robust_method = "huber",
         robust_tuning_method = "M")
}, error = function(err) {
  message("Huber-M. seed: ", seed, " error: ", err)
})
summary(res_gs_hb_M)
gsusie_coefficients(res_gs_hb_M)

### Comparison: SuSiE
res_su <- susie(X = X[, -ncol(X)], y = y)
summary(res_su)
gsusie_coefficients(res_su)

### Comparison: GLMNET
# Lasso
res_la <- cv.glmnet(x = X[, -ncol(X)], y, family = "poisson", alpha = 1)
coef_la <- data.frame(coef = as.numeric(coef(res_la, s = "lambda.min"))[-1])
row.names(coef_la) <- paste0("X", 1 : pp)
coef_la %>%
  filter(coef != 0) %>%
  arrange(desc(abs(coef))) %>%
  head(n = 10)


# Elastic net: alpha=0.5
res_en <- cv.glmnet(x = X[, -ncol(X)], y, family = "poisson", alpha = 0.5)
coef_en <- data.frame(coef = as.numeric(coef(res_en, s = "lambda.min"))[-1])
row.names(coef_en) <- paste0("X", 1 : pp)
coef_en %>%
  filter(coef != 0) %>%
  arrange(desc(abs(coef))) %>%
  head(n = 10)


p1 <- gsusie_pip_plot(res_gs_vn$pip[-(pp+1)], effect_idx) +
  labs(title = "GSuSiE (vanilla)", y = "PIP")
p1
ggsave(filename = "genotype-gsusie-vn.pdf",
       plot = p1, width = 9, height = 6,
       path = "./vignettes/motivating-example-genotype/")

p2 <- gsusie_pip_plot(res_gs_hb_M$pip[-(pp+1)], effect_idx) +
  labs(title = "GSuSiE (Huber-M)", y = "PIP")
p2
ggsave(filename = "genotype-gsusie-hb-M.pdf",
       plot = p2, width = 9, height = 6,
       path = "./vignettes/motivating-example-genotype/")

p3 <- gsusie_pip_plot(res_su$pip, effect_idx) +
  labs(title = "SuSiE", y = "PIP")
p3
ggsave(filename = "genotype-susie.pdf",
       plot = p3, width = 9, height = 6,
       path = "./vignettes/motivating-example-genotype/")

p4 <- gsusie_pip_plot(abs(coef_la$coef), effect_idx) +
  labs(title = "Lasso", y = "|coefficient|")
p4
ggsave(filename = "genotype-lasso.pdf",
       plot = p4, width = 9, height = 6,
       path = "./vignettes/motivating-example-genotype/")

p5 <- gsusie_pip_plot(abs(coef_en$coef), effect_idx) +
  labs(title = "Elastic Net", y = "|coefficient|")
p5
ggsave(filename = "genotype-elastic-net.pdf",
       plot = p5, width = 9, height = 6,
       path = "./vignettes/motivating-example-genotype/")

wrap_plots(p2, p3, p4, p5, nrow = 2)





gsusie_plot(res_gs_hb_M, y = "PIP")
gsusie_plot(res_su, y = "PIP")

# library(ggplotify)
# p6 <- as.ggplot(function() gsusie_plot(res_gs_hb_M, y = "PIP"))
