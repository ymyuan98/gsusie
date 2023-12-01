library(tidyverse)
# library(viridis)
# library(ggsci)
library(scales)
library(data.table)
library(stringr)
library(PRROC)
# library(patchwork)
setDTthreads(threads = 8)

########################################################################
## The colors are from pal_jco("default")(10)
## from `library(ggsci)`
## `scale_color_jco()`
myPalette <- c(
  "gs_vn"   = "#0073C2FF",
  "gs_hb_M" = "#EFC000FF",
  "su"      = "#868686FF",
  "la"      = "#CD534CFF",
  "en"      = "#7AA6DCFF",
  "rr"      = "#8F7700FF",
  "gs_hb_S" = "#003C67FF",
  "gs_fc30" = "#3B3B3BFF",
  "gs_fc5"  = "#A73030FF",
  "gs_fc1"  = "#4A6990FF"
)
########################################################################

`%&%` <- function(a, b) paste0(a, b)

read.sim_results <- function(filename) {
  sim_results <- readRDS(filename)
  out1 <- setDT(sim_results$pip)
  out2 <- setDT(sim_results$gn_coef)
  if (any(is.na(colnames(out2)))) {
    message("Warning: a column named NA")
    out2 <- out2[, .(la, en)]
  }
  out <- cbind(out1, out2)
  out[, label:=0]
  out[sim_results$setup$effect_idx, label:=1]
  out[, h2:=sim_results$setup$h2]
  out[, n:=sim_results$setup$n]
  out[, p:=sim_results$setup$p]
  out[, seed:=sim_results$setup$seed]
  return(out)
}

########################################################################

results.dir <- "./results/"
dir.name <- "2023-12-01-genotype-poisson"
results.filenames <- sort(list.files(path = results.dir%&%dir.name,
                                     full.names = T))

results.dt <- do.call(rbind, lapply(results.filenames, read.sim_results))
# head(results.dt)


########
## Specify methods for comparisons:

if (str_detect(dir.name, "binomial")) {
  print("Ha!!")
  results.dt <- results.dt[, c("gs_fc30", "gs_fc5", "gs_fc1"):= NULL]
  # "gs_hb_S", "gs_hb_M", "gs_bs_S", "gs_bs_M",
  # methods for comparisons
  methods_name <- c("gs_vn", "su", "la", "en")  # logistic
  # For logistic regression: remove the robust methods

} else if (str_detect(dir.name, "poisson")) {
  print("Ha Ha!!")
  # For poisson regression: remove gs_bs_S, gs_bs_M from comparison/evaluation
  # results.dt <- results.dt[, c("gs_bs_S", "gs_bs_M"):= NULL]
  ## methods for comparison
  methods_name <- c("gs_vn", "gs_hb_M", "su", "la", "en")  # poisson
}

########################################################################
# Compute AUPRC for each trial

compute.auprc <- function(res_col, label) {
  auc <- tryCatch({
    pr.curve(scores.class0 = res_col[label == 1],
             scores.class1 = res_col[label == 0])$auc.davis.goadrich
  }, error = function(err) {
    return(0)
  })

  return(auc)
}
## look at the cases where all pips are 0.

auprc <- function(resultsDT) {
  resultsDT <-
    melt(resultsDT, id.vars = c("seed", "n", "p", "h2", "label"),
         variable.name = "method", value.name = "pip") %>%
    filter(!is.na(pip)) %>%
    mutate(pip = abs(pip))
  resultsDT[, .(auprc=compute.auprc(pip, label)),
            by = .(method, seed, n, p, h2)]
}

########################################
## Function of lineplot of average AUPRC
lineplot.auprc <- function(auprcDT, .n, .p, .methods = NULL,
                           error_bar = TRUE) {
  if (is.null(.methods)) {
    .methods <- unique(auprcDT$method)
  }

  pt <- auprcDT[n == .n & p == .p & (method %in% .methods),
                .(m.auprc = mean(auprc),
                  se.auprc = sd(auprc)/sqrt(.N)),
                by = .(n, p, h2, method)][,
                `:=`(lower = max(m.auprc - se.auprc, 0),
                     upper = min(m.auprc + se.auprc, 1)),
                by = .(n, p, h2, method)] %>%
    mutate(method = fct_relevel(method, .methods)) %>%
    arrange(desc(method)) %>%
    ggplot(aes(x = h2, y = m.auprc, color = method))

  if (error_bar) {
    pt <- pt +
      geom_errorbar(aes(ymin = lower, ymax = upper, group = method),
                    width = 0.005, alpha = 0.8, linewidth = 0.6)
  }

  pt <- pt +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = myPalette) +
    labs(y = "Average AU-PRC", color = "Method", x = "h2") +
    theme_bw()

  return(pt)
}
##############################

file.name <- "AUPRC-plots/"%&%dir.name%&% "-DT-auprc3.txt.gz"
if (!file.exists(results.dir%&%file.name)) {
  AUPRC.dt <- auprc(results.dt)
  AUPRC.dt[, head(.SD, 2), by = .(method, h2)]
  fwrite(AUPRC.dt, file = results.dir%&%file.name)
} else {
  AUPRC.dt <- fread(results.dir %&% file.name)
}

## check set-up
sort(unique(AUPRC.dt[,seed-1] %/% 100)+1)
######################
## Plots

if (str_detect(dir.name, "genotype")) {
  print("Ha!!")
  ## genotype comb
  params_comb <- as.data.frame(
    expand.grid(n = c(500, 1000),
                p = c(2000, 5000, 10000)))
} else if (str_detect(dir.name, "synthetic")) {
  print("Ha Ha!!")
  ## synthetic comb
  params_comb <- as.data.frame(
    expand.grid(n = c(500, 800),
                p = c(1000, 3000, 5000)))
}

dir.create(results.dir%&%"AUPRC-plots/"%&%dir.name, recursive = T)
for (nr in 1 : nrow(params_comb)) {
  .n <- params_comb[nr, 1]
  .p <- params_comb[nr, 2]
  fig.title <- "n="%&% .n %&% ", p=" %&% .p
  fig.name <- "n"%&% .n %&%"-p"%&%.p%&%"-v3.pdf"

  fig <- lineplot.auprc(AUPRC.dt, .n = .n, .p = .p,
                        .method = methods_name,  ## remove for binomial cases
                        error_bar = T) +
    labs(title = fig.title)

  ggsave(filename = dir.name%&%"-"%&%fig.name,
         plot = fig,
         width = 9, height = 6,
         path = results.dir %&% "AUPRC-plots/"%&%dir.name)
}

rm(AUPRC.dt, fig)
gc()

################################################################################
################################################################################
## Aggregated precision-recall curves


aggregate.prcurve <- function(res_col, label) {
  auc <- tryCatch({
    pr.curve(scores.class0 = res_col[label == 1],
             scores.class1 = res_col[label == 0], curve = T)$curve
  }, error = function(err) {
    message(err)
    return(0)
  })
}

aggregate_prcurve_dt <- function(resultsDT, .n, .p, .h2, .methods = NULL) {

  resultsDT <- resultsDT[n == .n & p == .p & h2 == .h2] %>%
    melt(id.vars = c("seed", "n", "p", "h2", "label"),
         variable.name = "method", value.name = "pip")  %>%
    filter(!is.na(pip)) %>%
    mutate(pip = abs(pip))
  if (is.null(.methods)) {
    .methods <- unique(resultsDT[, method])
  }
  resultsDT <- resultsDT[method %in% .methods]

  outdt <- data.table()
  for (m in .methods) {
    temp <- as.data.table(
      aggregate.prcurve(resultsDT[method == m, pip],
                        resultsDT[method == m, label]))[, method := m]
    outdt <- rbind(outdt, temp)
    rm(temp)
    gc()
  }
  outdt <- outdt %>% setnames(., c("V1", "V2", "V3"),
                                 c("Recall", "Precision", "Threshold"))
  return(outdt)
}


dir.create(results.dir%&%"PRCurves-plots/"%&%dir.name%&%"-PRC-figs",
           recursive = T)
dir.create(results.dir%&%"PRCurves-plots/"%&%dir.name%&%"-PRC-DT",
           recursive = T)

hh2 <- c(0.01, 0.1, 0.2, 0.3, 0.4)
for (.h2 in hh2) {
  for (nr in 1 : nrow(params_comb)) {
    .n <- params_comb[nr, 1]
    .p <- params_comb[nr, 2]
    fig.title <- "n="%&%.n%&%", p="%&%.p%&%", h2="%&%.h2

    if (.h2 == 0.01) {
      fig.name <- "PRC-n"%&%.n%&%"-p"%&%.p%&%"-hz1"%&%"-v3"
    } else {
      fig.name <- "PRC-n"%&%.n%&%"-p"%&%.p%&%"-h"%&%(10*.h2)%&%"-v3"
    }

    if (!file.exists(results.dir %&% "PRCurves-plots/"%&%dir.name%&%"-PRC-DT/"%&%
                     fig.name%&%".txt.gz")) {
      prcurve_dt <- aggregate_prcurve_dt(results.dt, .n, .p, .h2,
                                         .methods = methods_name)
      fwrite(prcurve_dt,
             file = results.dir %&% "PRCurves-plots/"%&%dir.name%&%"-PRC-DT/"%&%
               fig.name%&%".txt.gz")
    } else {
      prcurve_dt <-
        fread(results.dir %&% "PRCurves-plots/"%&%dir.name%&%"-PRC-DT/"%&%
                fig.name%&%".txt.gz")
    }

    gc()

    fig <- prcurve_dt %>%
      mutate(method = fct_relevel(method, methods_name)) %>%
      arrange(desc(method)) %>%
      ggplot(aes(x = Recall, y = Precision, color = method)) +
      geom_line(linewidth = 1, alpha = 0.8) +
      scale_color_manual(values = myPalette) +
      xlim(0, 1) + ylim(0, 1) +
      labs(title = fig.title, color = "Method") +
      theme_bw()

    ggsave(filename = dir.name%&%"-"%&%fig.name%&%".pdf",
           plot = fig,
           width = 9, height = 6,
           path = results.dir %&% "PRCurves-plots/"%&%dir.name%&%"-PRC-figs")

    rm(fig, prcurve_dt)
    gc()

  }
}
rm(results.dt)
gc()

