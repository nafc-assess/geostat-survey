
library(SimSurvey)
library(furrr)
plan(multisession, workers = 10)

source("./pop_cod_fn.R")
source("./pop_yellowtail_fn.R")

## Mostly default SimSurvey settings -------------------------------------------

tic <- Sys.time()
default <- furrr::future_map_dfr(seq_len(1000), function(i) {
  sim <- sim_abundance() %>%
    sim_distribution(grid = make_grid(res = c(10, 10))) %>%
    sim_survey(set_den = 10/1000) %>%
    run_strat() %>%
    strat_error()
  sim$total_strat_error
},
.options = furrr::furrr_options(seed = 123, packages = "SimSurvey"))
toc <- Sys.time()
toc - tic

default_e <- default$error
hist(default_e, breaks = 200)
abline(v = 0, col = "red", lwd = 2)
mean(default_e)
##     610975.8    with     10 sims
##    -914137      with    100 sims
##    -400311.2    with   1000 sims
##    -5491.609    with   1000 sims when set_den = 10/1000

## Yellowtail ------------------------------------------------------------------

yt <- furrr::future_map_dfr(seq_len(1000), function(i) {
  sim <- population_yellowtail(i, n_sims = 1, set_den = 10/1000) %>%
    run_strat() %>%
    strat_error()
  sim$total_strat_error
  },
  .options = furrr::furrr_options(seed = 891, packages = "SimSurvey"))

yt_e <- yt$error
hist(yt_e, breaks = 200)
abline(v = 0, col = "red", lwd = 2)
mean(yt_e)
##     1179984       with     10 sims
##   -12064882       with    100 sims
##   -11867666       with   1000 sims
##   -15052764       with  10000 sims
##   -27799965       with  10000 sims with set_den = 5/1000
##       18799.67    with   1000 sims with set_den = 10/1000
##        3974.876   with  10000 sims with set_den = 10/1000


## Cod -------------------------------------------------------------------------

cod <- furrr::future_map_dfr(seq_len(10000), function(i) {
  sim <- population_cod(i, n_sims = 1, set_den = 10/1000) %>%
    run_strat() %>%
    strat_error()
  sim$total_strat_error
},
.options = furrr::furrr_options(seed = 438, packages = "SimSurvey"))

cod_e <- cod$error
hist(cod_e, breaks = 200)
abline(v = 0, col = "red", lwd = 2)
mean(cod_e)
##    1544895      with     10 sims
##    -388603.1    with    100 sims
##     -68757.93   with   1000 sims
##       1416.984  with   1000 sims and set_den = 10/1000
##       -253.761  with  10000 sims and set_den = 10/1000


