library(tidyverse)
library(data.table)


gsusie_pip_plot <- function(gs_res_pip, effect_idx) {
  tibble(pip = gs_res_pip) %>%
    mutate(id = 1 : n()) %>%
    mutate(true_label = if_else(id %in% effect_idx, "causal", "noncausal")) %>%
    arrange(desc(true_label)) %>%
    ggplot(aes(x = id, y = pip, color = as.factor(true_label))) +
    geom_point()
}


library("patchwork")
p1 <- gsusie_plot(res_gs_vn$pip, effect_idx) + ggtitle("vanilla")
p2 <- gsusie_plot(res_gs_hb_M$pip, effect_idx) + ggtitle("huber_M")
p3 <- gsusie_plot(res_gs_hb_S$pip, effect_idx) + ggtitle("huber_S")
p4 <- gsusie_plot(res_gs_bs_M$pip, effect_idx) + ggtitle("bisquare-M")
p5 <- gsusie_plot(res_gs_bs_S$pip, effect_idx) + ggtitle("bisquare_S")
p6 <- gsusie_plot(res_gs_fc30$pip, effect_idx) + ggtitle("fraction-0.3")
p7 <- gsusie_plot(res_gs_fc5$pip, effect_idx) + ggtitle("fraction-0.05")
p8 <- gsusie_plot(res_gs_fc1$pip, effect_idx) + ggtitle("fraction-0.01")
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8,
           nrow = 4)


res_su <- susie(X[, -(pp+1)], y)
p9 <- gsusie_pip_plot(res_su$pip, effect_idx)
p1| p2 | p3 | p9
p3

p1 | p9




################################################################################



## an example
## n = 500, p = 1000, h2 = 0.5, seed = 30000
effect_idx <- output$setup$effect_idx + 1
p1 <- gsusie_pip_plot(output$gs_pip[,"gs_vn"],   effect_idx) + ggtitle("vanilla")
p2 <- gsusie_pip_plot(output$gs_pip[,"gs_hb_M"], effect_idx) + ggtitle("huber_M")
p3 <- gsusie_pip_plot(output$gs_pip[,"gs_hb_S"], effect_idx) + ggtitle("huber_S")
p4 <- gsusie_pip_plot(output$gs_pip[,"gs_bs_M"], effect_idx) + ggtitle("bisquare-M")
p5 <- gsusie_pip_plot(output$gs_pip[,"gs_bs_S"], effect_idx) + ggtitle("bisquare_S")
p6 <- gsusie_pip_plot(output$gs_pip[,"gs_fc30"], effect_idx) + ggtitle("fraction-0.3")
p7 <- gsusie_pip_plot(output$gs_pip[,"gs_fc5"],  effect_idx) + ggtitle("fraction-0.05")
p8 <- gsusie_pip_plot(output$gs_pip[,"gs_fc1"],  effect_idx) + ggtitle("fraction-0.01")
wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9,
           nrow = 3)

p1 | p9


library(PRROC)
library(data.table)
results.dir <- "./results/2023-11-01-synthetic-poisson/"
results.filenames <- sort(list.files(path = results.dir, full.names = T))

output <- readRDS(results.filenames[1])

results.dt <- do.call(rbind, lapply(results.filenames, read.pip))

read.pip <- function(filename) {
  output <- readRDS(filename)
  PIPs <- setDT(output$gs_pip)
  PIPs <- PIPs[-1,]
  PIPs[, label:=0]
  PIPs[output$setup$effect_idx, label:=1]
  PIPs[, h2:=output$setup$h2]
  PIPs[, n:=output$setup$n]
  PIPs[, p:=output$setup$p]
  return(PIPs)
}


## PR curve for each h2, method, n, p




