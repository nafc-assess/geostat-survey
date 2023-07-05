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



design_index_yellowtail <- merge(design_index_yellowtail, true_yellowtail, by=c("pop", "year", "species"))

design_r <-
  design_index_yellowtail |>
  filter(year > 10) |>
  ungroup() |>
  group_by(pop, type, scenario, species) |>
  summarise(MRE = mean((N - true) / true),
            ME = mean(N - true),
            difflog = mean(log(N) - log(true)),
            RMSE = sqrt(mean((log(N) - log(true))^2)))

hist(design_r$difflog)
