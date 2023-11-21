library(tidyverse)
library(data.table)
library(PRROC)
library(patchwork)


results.dir <- "./results/2023-11-11-genotype-poisson/"
results.filenames <- sort(list.files(path = results.dir, full.names = T))
# output <- readRDS(results.filenames[1])

read.sim_results <- function(filename) {
  sim_results <- readRDS(filename)
  out1 <- setDT(sim_results$pip)
  out2 <- setDT(sim_results$gn_coef)
  out <- cbind(out1, out2)
  out[, label:=0]
  out[sim_results$setup$effect_idx, label:=1]
  out[, h2:=sim_results$setup$h2]
  out[, n:=sim_results$setup$n]
  out[, p:=sim_results$setup$p]
  out[, seed:=sim_results$setup$seed]
  return(out)
}

results.dt <- do.call(rbind, 
                      lapply(results.filenames[800:1800], 
                             read.sim_results))
# head(results.dt)
# 
# results.dt[label == 1]
############################
# Compute AUPRC for each trial

## Column `seed` defines a trial

# t1 <- pr.curve(scores.class0 = results.dt$gs_vn,
#                weights.class0 = results.dt$label,
#                curve = T)

compute.auprc <- function(res_col, label) {
    auc <- tryCatch({
      pr.curve(scores.class0 = res_col, 
               weights.class0 = label,
               curve = T)$auc.integral  
    }, error = function() {
      return(0)
    })
  
    return(auc)
}

# compute.auprc(results.dt$gs_hb_M, results.dt$label)

# t1 <- pr.curve(scores.class0 = results.dt[seed==4950, gs_hb_M], 
#                weights.class0 = results.dt[seed==4950, label],
#                curve = T)
# plot(t1)


auprc <- function(resultsDT) {
  resultsDT <- copy(resultsDT) 
  method_names <- colnames(resultsDT)[1:12]
  resultsDT[, `:=`(la = abs(la), rr = abs(rr), en = abs(en))]
  resultsDT[, lapply(.SD, compute.auprc, label = label), 
            by = .(seed, n, p, h2),
            .SDcols = method_names]
}

# AUPRC.dt <- auprc(results.dt)

AUPRC.dt.long <- 
  melt(AUPRC.dt, id.vars = c("seed", "n", "p", "h2"),
       variable.name = "method", value.name = "auprc")
  
dat <- AUPRC.dt.long[, .(m.auprc = mean(auprc), 
                         sd.auprc = sd(auprc), 
                         upper = min(mean(auprc) + sd(auprc), 1),
                         lower = max(mean(auprc) - sd(auprc), 0)), 
                     by = .(n, p, h2, method)]



lineplot.auprc <- function(auprc.dat, .n, .p, .method = NULL,
                           error_bar = TRUE) {
  if (is.null(.method)) {
    .method = unique(auprc.dat$method)
  }
  pt <- auprc.dat[n == .n & p == .p & (method %in% .method)] %>%
    ggplot(aes(x = h2, y = m.auprc, color = method)) + 
    geom_line(linewidth = 0.8, alpha = 0.8) + 
    geom_point(size = 2, alpha = 0.8) +
    labs(y = "average AU-PRC")
  if (error_bar) {
    pt <- pt + 
      geom_errorbar(aes(ymin = lower, ymax = upper, group = method),
                    width = 0.005, alpha = 0.8, linewidth = 0.6)
  }
  return(pt)
}

p1 <- lineplot.auprc(dat, .n = 1000, .p = 10000, 
                     .method = c("gs_vn", "gs_hb_M", "su", "la", "en"),
                     error_bar = T)  
p1

p2 <- lineplot.auprc(dat, .n = 500, .p = 10000, 
                     .method = c("gs_vn", "gs_hb_M", "su", "la", "en"),
                     error_bar = F)  
p2



# AUPRC.dt.long[, head(.SD, 2), by = .(n, p, h2, method)]


##############################
t1 <- pr.curve(scores.class0 = results.dt[seed==4950, gs_hb_M], 
               weights.class0 = results.dt[seed==4950, label],
               curve = T)
plot(t1$curve[, 1], t1$curve[, 2])
t2 <- pr.curve(scores.class0 = results.dt[seed==4951, gs_hb_M], 
               weights.class0 = results.dt[seed==4951, label],
               curve = T)
plot(t2)
plot(t2$curve[, 1], t2$curve[, 2])
t3 <- pr.curve(scores.class0 = results.dt[seed==4953, gs_hb_M], 
               weights.class0 = results.dt[seed==4953, label],
               curve = T)
plot(t3)
plot(t3$curve[, 1], t3$curve[, 2])

t4 <- rbind(t1$curve, t2$curve, t3$curve)
plot(t4[, 1], t4[, 2])


extract.prc.dt <- function(res_col, label) {
  prc.df <- tryCatch({
    setDT(pr.curve(res_col, label, curve = TRUE)$curve)
  }, error = function(){
    return(data.table(0, 0, 0))
  })
  return(prc.df)
}

prcdt <- function(resultsDT) {
  resultsDT <- copy(resultsDT)
  method_names <- colnames(resultsDT)[1:12]
  resultsDT[, `:=`(la = abs(la), rr = abs(rr), en = abs(en))]
  prc.dt <- resultsDT[, lapply(.SD, extract.prc.dt, label = label), 
                      by = .(seed, n, p, h2), 
                      .SDcols = method_names]
}

# auprc <- function(resultsDT) {
#   resultsDT <- copy(resultsDT) 
#   method_names <- colnames(resultsDT)[1:12]
#   resultsDT[, `:=`(la = abs(la), rr = abs(rr), en = abs(en))]
#   resultsDT[, lapply(.SD, compute.auprc, label = label), 
#             by = .(seed, n, p, h2),
#             .SDcols = method_names]
# }

