####################################################################
## Parameters settings for synthetic (all-random) data

nn <- c(500, 1000)
pp <- c(1000, 3000, 5000)
h2 <- c(0.01, 0.1, 0.2, 0.4)
n_effect_vars <- c(10, 5)


param_comb <- as.data.frame(expand.grid(n = nn, p = pp, h2 = h2,
                                        n_effect_vars = n_effect_vars)) %>%
  mutate(ID = 1 : n()) %>%
  slice(rep(1 : n(), each = 100)) %>% # (repeat each simulation 100 times)
  mutate(seed = 1 : n()) %>%
  select(ID, seed, n, p, h2, n_effect_vars)


write.table(param_comb, file = "./tests/params.txt",
            row.names = F,col.names = F)
# save as `params.txt` in the remote host.


####################################################################
## Parameters settings for genomic data

nn <- c(500, 1000)
pp <- c(2000, 5000, 10000)
h2 <- c(0.01, 0.1, 0.3)
n_effect_vars <- c(1, 3, 5)

param_comb <- as.data.frame(expand.grid(n = nn, p = pp, h2 = h2,
                                        n_effect_vars = n_effect_vars)) %>%
  mutate(ID = 1 : n()) %>%
  slice(rep(1 : n(), each = 100)) %>% # (repeat each simulation 100 times)
  mutate(seed = 1 : n()) %>%
  select(ID, seed, n, p, h2, n_effect_vars)

param_comb %>%
  filter(ID == 54)

## 54 cases in total, repeat each simulation 100 times
## ID 1-18,  (row 1-1800),    n_effect_vars <- 1
## ID 19-36, (row 1801-3600), n_effect_vars <- 3
## ID 37-54, (row 3601-5400), n_effect_vars <- 5

write.table(param_comb, file = "./tests/params-genom.txt",
          row.names = F,col.names = F)
