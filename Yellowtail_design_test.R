# Yellowtail-like species design and model based indices calculation

#############  Packages
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)
library(ggpubr)
plan(multisession)

### Load functions
source("./pop_yellowtail_fn.R")
source("./model_run_fn.R")
source("./bootstrapping_fn.R")
source("./data_prep_fn.R")

############               #############

# BASE SCENARIO

#############               #############


#############  Population simulations

set.seed(1)
survey_yellowtail <- furrr::future_map(seq_len(100), n_sims = 1, population_yellowtail, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

#save(survey_yellowtail, file = "./data/survey_yellowtail_base.Rdata")

############# True abundance

true_yellowtail <- map_df(seq_along(survey_yellowtail),function(i){
  tibble(year = unique(survey_yellowtail[[i]]$sp_N$year), true = as.numeric(colSums(survey_yellowtail[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like") |>
    group_by(pop)})

############# Design-based index

design_index_yellowtail <- map_df(seq_along(survey_yellowtail),function(i){
  strat <- SimSurvey::run_strat(survey_yellowtail[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Base") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()
