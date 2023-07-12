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
# plan(multisession)
plan(multisession, workers = 10L)

### Load functions
source("./pop_yellowtail_fn.R")
source("./model_run_fn.R")
source("./bootstrapping_fn.R")
source("./data_prep_fn.R")

# set to 1 so data.table doesn't spawn forks on forks
# usethis::edit_r_environ()
Sys.getenv("OMP_THREAD_LIMIT")
############               #############

            # BASE SCENARIO

#############               #############

message("BASE SCENARIO")

#############  Population simulations

set.seed(1)
survey_yellowtail <- furrr::future_map(seq_len(100L), n_sims = 1, population_yellowtail, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
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

############# Data prep: Separating populations

setdet_yellowtail <- map(survey_yellowtail, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail)){
  setdet_yellowtail[[i]]$pop <- as.numeric(i)
}

setdet_yellowtail <- lapply(setdet_yellowtail, function(x) split(x, x$sim)) |> flatten()

#save(setdet_yellowtail, file = "./data/setdet_yellowtail_base.Rdata")

############# Bootstrapped index

boot_index_yellowtail <- furrr::future_map_dfr(setdet_yellowtail, boot_wrapper, reps = 1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Base")

############# sdmTMB

### Data prep

sdm_data_yellowtail <- furrr::future_map(setdet_yellowtail, sdm_data_fn)

mesh_sdm_yellowtail <- purrr::map(sdm_data_yellowtail, mesh_sdm_fn)

sdm_newdata_yellowtail <- sdm_newdata_fn(survey_yellowtail[[1]], sdm_data_yellowtail[[1]]) ### since all populations has the same prediction area, newdata is same for all.

### Models

### IID + NB2
sdm_NB2_IID_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                       formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Base",
                                                       species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                             formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Base",
                                                             species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Base",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                            formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Base",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Base",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail <- furrr::future_map2_dfr(sdm_data_yellowtail, mesh_sdm_yellowtail,
                                                            formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Base",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)

#############                          #############

           # SET DENSITY REDUCTION SCENARIO

#############                          #############

message("SET DENSITY REDUCTION SCENARIO")
# Subsetting the sets by year
# 30 % effort reduction after Year 10

setdet_r30_fn <- function(x, year_cutoff) {
  a <- x |>
    filter(year >= year_cutoff) |>
    group_by(sim, year, strat) |>
    sample_frac(0.7) |>
    filter(unique(strat_sets) >= 2) ### additional rule set >=2
  b <- x|>
    filter(year < year_cutoff)
  ab <- rbind(a,b) |>
    as.data.table()
  return(ab)
}

setdet_yellowtail <- map(survey_yellowtail, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail)){
  setdet_yellowtail[[i]]$pop <- as.numeric(i)} ### sim level included

set.seed(45)
setdet_r30 <- map(setdet_yellowtail, setdet_r30_fn, year_cutoff = 11)

## updating the survey
survey_r30 <- survey_yellowtail
for(i in seq_along(survey_r30)) {
  survey_r30[[i]]$setdet <- setdet_r30[[i]]
}

############# Design-based index

## Running design-based calculation

design_index_r30 <- map_df(seq_along(survey_r30),function(i){
  strat <- SimSurvey::run_strat(survey_r30[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Set reduction") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_yellowtail_r30 <- map(survey_r30, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail_r30)){
  setdet_yellowtail_r30[[i]]$pop <- as.numeric(i)
}

setdet_yellowtail_r30 <- lapply(setdet_yellowtail_r30, function(x) split(x, x$sim)) |> flatten()

############# Bootstrapped index

boot_index_yellowtail_r30 <- furrr::future_map_dfr(setdet_yellowtail_r30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Set reduction")

#################  sdmTMB

### Data prep

sdm_data_yellowtail_r30 <- furrr::future_map(setdet_yellowtail_r30, sdm_data_fn)

mesh_sdm_yellowtail_r30 <- purrr::map(sdm_data_yellowtail_r30, mesh_sdm_fn)

#mesh_sdm_yellowtail_r30  <- map(sdm_data_yellowtail_r30, mesh_sdm_fn)

### IID + NB2
sdm_NB2_IID_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                       formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Set reduction",
                                                       species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                             formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Set reduction",
                                                             species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Set reduction",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                            formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Set reduction",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Set reduction",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_r30 <- furrr::future_map2_dfr(sdm_data_yellowtail_r30, mesh_sdm_yellowtail_r30,
                                                            formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Set reduction",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)

#############               #############

        # STRATA REMOVAL SCENARIO

#############               #############

message("STRATA REMOVAL SCENARIO")

setdet_yellowtail <- map(survey_yellowtail, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail)){
  setdet_yellowtail[[i]]$pop <- as.numeric(i)
}

# Getting all samples before year 10

samples10 <- map(setdet_yellowtail, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_yellowtail, function(x) {subset(x, year >= 11)})

blocked_samples <- map(seq_along(sample_a_y10), function(i) {
  blocked_samples <- sample_a_y10[[i]] |>
    filter(!strat %in% c(2,3,4,5,17,18,19,20))})

combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_SR <- survey_yellowtail
for(i in seq_along(survey_SR)){
  survey_SR[[i]]$setdet <- combined_samples_10[[i]]}

setdet_yellowtail_SR <- map(survey_SR, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_yellowtail_SR)){
  setdet_yellowtail_SR[[i]]$pop <- as.numeric(i)
}

################# Design-based index

design_index_SR <- map_df(seq_along(survey_SR),function(i){
  strat <- SimSurvey::run_strat(survey_SR[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Strata removal") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_yellowtail_SR <- map(survey_SR, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_yellowtail_SR)){
  setdet_yellowtail_SR[[i]]$pop <- as.numeric(i)
}

setdet_yellowtail_SR <- lapply(setdet_yellowtail_SR, function(x) split(x, x$sim)) |> flatten()

############# Bootstrapped index

boot_index_yellowtail_SR <- furrr::future_map_dfr(setdet_yellowtail_SR, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Strata removal")

################# sdmTMB

### Data prep

sdm_data_yellowtail_SR <- furrr::future_map(setdet_yellowtail_SR, sdm_data_fn)

mesh_sdm_yellowtail_SR  <- purrr::map(sdm_data_yellowtail_SR, mesh_sdm_fn)


### IID + NB2
sdm_NB2_IID_index_yellowtail_SR  <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                       formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Strata removal",
                                                       species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_SR  <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                             formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Strata removal",
                                                             species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_SR  <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Strata removal",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_SR  <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                            formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Strata removal",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_SR <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Strata removal",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)

#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_SR <- furrr::future_map2_dfr(sdm_data_yellowtail_SR, mesh_sdm_yellowtail_SR,
                                                            formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Strata removal",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)


#############                           #############

           # AREA BLOCKED REDUCTION SCENARIO

#############                           #############

message("AREA BLOCKED REDUCTION SCENARIO")

################# Blocking the setdets

setdet_yellowtail <- map(survey_yellowtail, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail)){
  setdet_yellowtail[[i]]$pop <- as.numeric(i)
}

# Getting all samples before year 10

samples10 <- map(setdet_yellowtail, function(x) {subset(x, year < 11)}) ### same for any perc

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_yellowtail, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples <- map(sample_a_y10, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_b30 <- survey_yellowtail
for( i in seq_along(survey_b30)){
  survey_b30[[i]]$setdet <- combined_samples[[i]]}

################# Design-based index

design_index_b30 <- map_df(seq_along(survey_b30),function(i){
  strat <- SimSurvey::run_strat(survey_b30[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Area blocked") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_yellowtail_b30 <- map(survey_b30, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_yellowtail_b30)){
  setdet_yellowtail_b30[[i]]$pop <- as.numeric(i)
}

setdet_yellowtail_b30 <- lapply(setdet_yellowtail_b30, function(x) split(x, x$sim)) |> flatten()

############# Bootstrapped index

boot_index_yellowtail_b30 <- furrr::future_map_dfr(setdet_yellowtail_b30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Area blocked")

################# sdmTMB

### Data prep

sdm_data_yellowtail_b30 <- furrr::future_map(setdet_yellowtail_b30, sdm_data_fn)

mesh_sdm_yellowtail_b30  <- map(sdm_data_yellowtail_b30, mesh_sdm_fn)

### IID + NB2
sdm_NB2_IID_index_yellowtail_b30  <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                            formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Area blocked",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_b30  <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                                  formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Area blocked",
                                                                  species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_b30  <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                             formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Area blocked",
                                                             species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)

### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_b30  <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                                   formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Area blocked",
                                                                   species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                            formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Area blocked",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)

#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                                  formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Area blocked",
                                                                  species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)


#############                           #############

# Scenario 5: RECOVERY

#############                           #############

message("RECOVERY SCENARIO")

#############  Population simulations

set.seed(1)
survey_yellowtail_rec <- furrr::future_map(seq_len(100), n_sims = 1, population_yellowtail_recovery, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

############# True abundance

true_yellowtail_rec <- map_df(seq_along(survey_yellowtail_rec),function(i){
  tibble(year = unique(survey_yellowtail_rec[[i]]$sp_N$year), true = as.numeric(colSums(survey_yellowtail_rec[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like") |>
    group_by(pop)})

############# Data prep: Separating populations

setdet_yellowtail_rec <- map(survey_yellowtail_rec, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail_rec)){
  setdet_yellowtail_rec[[i]]$pop  <- as.numeric(i)}

# Getting all samples before year 10

samples10_rec <- map(setdet_yellowtail_rec, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10_rec <- map(setdet_yellowtail_rec, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples_rec <- map(sample_a_y10_rec, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples_rec <- map(seq_along(blocked_samples_rec),function(i) {rbind(blocked_samples_rec[[i]], samples10_rec[[i]])})

################### Updating survey lists

survey_rec <- survey_yellowtail_rec
for( i in seq_along(survey_rec)){
  survey_rec[[i]]$setdet  <- combined_samples_rec[[i]]}

################# Design-based index

design_index_rec <- map_df(seq_along(survey_rec),function(i){
  strat  <- SimSurvey::run_strat(survey_rec[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Recovery") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_yellowtail_rec <- map(survey_rec, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_yellowtail_rec)){
  setdet_yellowtail_rec[[i]]$pop  <- as.numeric(i)
}

setdet_yellowtail_rec <- lapply(setdet_yellowtail_rec, function(x) split(x, x$sim)) |> flatten()

############# Bootstrapped index

boot_index_yellowtail_rec <- furrr::future_map_dfr(setdet_yellowtail_rec, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Recovery")
gc()

############# sdmTMB

### Data prep

sdm_data_yellowtail_rec <- furrr::future_map(setdet_yellowtail_rec, sdm_data_fn)

mesh_sdm_yellowtail_rec <- furrr::future_map(sdm_data_yellowtail_rec, mesh_sdm_fn)

### Models

### IID + NB2
sdm_NB2_IID_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                           formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Recovery",
                                                           species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                                 formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Recovery",
                                                                 species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                          formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Recovery",
                                                          species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                                formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Recovery",
                                                                species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                          formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Recovery",
                                                          species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_rec <- furrr::future_map2_dfr(sdm_data_yellowtail_rec, mesh_sdm_yellowtail_rec,
                                                                formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Recovery",
                                                                species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#############


#############                           #############

# Scenario 6: RECOVERY WITH SPILLOVER EFFECT

#############                           #############

message("RECOVERY WITH SPILLOVER EFFECT")

#############  Population simulations

set.seed(1)
survey_yellowtail_so <- furrr::future_map(seq_len(100), n_sims = 1, population_yellowtail_spillover, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

############# True abundance

true_yellowtail_so <- map_df(seq_along(survey_yellowtail_so),function(i){
  tibble(year = unique(survey_yellowtail_so[[i]]$sp_N$year), true = as.numeric(colSums(survey_yellowtail_so[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like") |>
    group_by(pop)})

############# Data prep: Separating populations

setdet_yellowtail_so <- map(survey_yellowtail_so, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail_so)){
  setdet_yellowtail_so[[i]]$pop  <- as.numeric(i)}

# Getting all samples before year 10

samples10_so <- map(setdet_yellowtail_so, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10_so <- map(setdet_yellowtail_so, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples_so <- map(sample_a_y10_so, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples_so <- map(seq_along(blocked_samples_so),function(i) {rbind(blocked_samples_so[[i]], samples10_so[[i]])})

################### Updating survey lists

survey_so <- survey_yellowtail_so
for( i in seq_along(survey_so)){
  survey_so[[i]]$setdet  <- combined_samples_so[[i]]}

################# Design-based index

design_index_so <- map_df(seq_along(survey_so),function(i){
  strat  <- SimSurvey::run_strat(survey_so[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Recovery + Spillover") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Yellowtail-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_yellowtail_so <- map(survey_so, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_yellowtail_so)){
  setdet_yellowtail_so[[i]]$pop  <- as.numeric(i)
}

setdet_yellowtail_so <- lapply(setdet_yellowtail_so, function(x) split(x, x$sim)) |> flatten()

############# Bootstrapped index

boot_index_yellowtail_so <- furrr::future_map_dfr(setdet_yellowtail_so, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Yellowtail-like", scenario = "Recovery + Spillover")
gc()

############# sdmTMB

### Data prep

sdm_data_yellowtail_so <- furrr::future_map(setdet_yellowtail_so, sdm_data_fn)

mesh_sdm_yellowtail_so <- furrr::future_map(sdm_data_yellowtail_so, mesh_sdm_fn)

### Models

### IID + NB2
sdm_NB2_IID_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                          formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Recovery + Spillover",
                                                          species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                                formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Recovery + Spillover",
                                                                species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                         formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Recovery + Spillover",
                                                         species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                               formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Recovery + Spillover",
                                                               species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                         formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Recovery + Spillover",
                                                         species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_so <- furrr::future_map2_dfr(sdm_data_yellowtail_so, mesh_sdm_yellowtail_so,
                                                               formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Recovery + Spillover",
                                                               species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)

#############               #############

# COMBINING ALL RESULTS

#############               #############

index_yellowtail_all_scenarios <- do.call(bind_rows, mget(ls(pattern = "index")))
index_yellowtail_all_scenarios <- merge(index_yellowtail_all_scenarios, true_yellowtail, by=c("pop", "year", "species"))
save(index_yellowtail_all_scenarios, file = "./data/index_yellowtail_all_scenarios.Rdata")
