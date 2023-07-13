# Cod-like species design and model based indices calculation

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
#plan(multisession)
plan(multisession, workers = 10L)

ITER <- 200L

### Load functions
source("./pop_cod_fn.R")
source("./model_run_fn.R")
source("./bootstrapping_fn.R")
source("./data_prep_fn.R")


# set to 1 so data.table doesn't spawn forks on forks
# usethis::edit_r_environ()
Sys.getenv("OMP_THREAD_LIMIT")

#############               #############

# Scenario 1: BASE SCENARIO #

#############               #############

message("BASE SCENARIO")

#############  Population simulations

set.seed(1)
survey_cod <- furrr::future_map(seq_len(ITER), n_sim = 1, population_cod, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

#save(survey_cod, file = "./data/survey_cod_base.Rdata")

############# True abundance

true_cod <- map_df(seq_along(survey_cod),function(i){
  tibble(year = unique(survey_cod[[i]]$sp_N$year), true = as.numeric(colSums(survey_cod[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Cod-like") |>
    group_by(pop)})

############# Design-based index

design_index_cod <- map_df(seq_along(survey_cod),function(i){
  strat <- SimSurvey::run_strat(survey_cod[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Base") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations

setdet_cod <- map(survey_cod, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod)){
  setdet_cod[[i]]$pop <- as.numeric(i)}

setdet_cod <- lapply(setdet_cod, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod, file = "./data/setdet_cod_base.Rdata")

############# Bootstrapped index

boot_index_cod <- furrr::future_map_dfr(setdet_cod, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Base")
gc()

############# sdmTMB

### Data prep

inla_mesh <- make_standard_mesh() # done once

sdm_data_cod <- furrr::future_map(setdet_cod, sdm_data_fn)

mesh_sdm_cod <- purrr::map(sdm_data_cod, mesh_sdm_fn, existing_mesh = inla_mesh)

sdm_newdata_cod <- sdm_newdata_fn(survey_cod[[1]], sdm_data_cod[[1]]) ### since all populations has the same prediction area, newdata is same for all.

#save(sdm_newdata_cod, file = "./data/sdm_newdata_cod.Rdata")

### Models

### IID + NB2
sdm_NB2_IID_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                       formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Base",
                                                       species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                             formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Base",
                                                             species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Base",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                            formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Base",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Base",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod <- furrr::future_map2_dfr(sdm_data_cod, mesh_sdm_cod,
                                                            formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Base",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############               #############

# Scenario 2: SET DENSITY REDUCTION #

#############               #############

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

setdet_cod <- map(survey_cod, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod)){
  setdet_cod[[i]]$pop <- as.numeric(i)}

set.seed(45)
setdet_r30 <- map(setdet_cod, setdet_r30_fn, year_cutoff = 11)

## updating the survey

survey_r30 <- survey_cod
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
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_r30 <- map(survey_r30, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod_r30)){
  setdet_cod_r30[[i]]$pop <- as.numeric(i)
}

setdet_cod_r30 <- lapply(setdet_cod_r30, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod_r30, file = "./data/setdet_cod_r30.Rdata")
############# Bootstrapped index

boot_index_cod_r30 <- furrr::future_map_dfr(setdet_cod_r30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Set reduction")

#################  sdmTMB

### Data prep

sdm_data_cod_r30 <- furrr::future_map(setdet_cod_r30, sdm_data_fn)

mesh_sdm_cod_r30  <- purrr::map(sdm_data_cod_r30, mesh_sdm_fn, existing_mesh = inla_mesh)

### IID + NB2
sdm_NB2_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                    formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Set reduction",
                                                    species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                          formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Set reduction",
                                                          species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Set reduction",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                         formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Set reduction",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Set reduction",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                         formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Set reduction",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############               #############

# Scenario 3: STRATA REMOVAL SCENARIO

#############               #############

message("STRATA REMOVAL SCENARIO")

setdet_cod <- map(survey_cod, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod)){
  setdet_cod[[i]]$pop <- as.numeric(i)
}

# Getting all samples before year 10

samples10 <- map(setdet_cod, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_cod, function(x) {subset(x, year >= 11)})

blocked_samples <- map(seq_along(sample_a_y10), function(i) {
  blocked_samples <- sample_a_y10[[i]] |>
    filter(!strat %in% c(2,3,4,5,17,18,19,20))})

combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_SR <- survey_cod
for(i in seq_along(survey_SR)){
  survey_SR[[i]]$setdet <- combined_samples_10[[i]]}

################# Design-based index

design_index_SR <- map_df(seq_along(survey_SR),function(i){
  strat <- SimSurvey::run_strat(survey_SR[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Strata removal") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_SR <- map(survey_SR, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_cod_SR)){
  setdet_cod_SR[[i]]$pop <- as.numeric(i)
}

setdet_cod_SR <- lapply(setdet_cod_SR, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod_SR, file = "./data/setdet_cod_SR.Rdata")

############# Bootstrapped index

boot_index_cod_SR <- furrr::future_map_dfr(setdet_cod_SR, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Strata removal")

################# sdmTMB

### Data prep

sdm_data_cod_SR <- furrr::future_map(setdet_cod_SR, sdm_data_fn)
#save(sdm_data_cod_SR, file = "./data/sdm_data_cod_SR.Rdata")

mesh_sdm_cod_SR  <- purrr::map(sdm_data_cod_SR, mesh_sdm_fn, existing_mesh = inla_mesh)
#save(mesh_sdm_cod_SR, file = "./data/mesh_sdm_cod_SR.Rdata")

### IID + NB2
sdm_NB2_IID_index_cod_SR  <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                    formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Strata removal",
                                                    species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)

### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_SR  <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                          formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Strata removal",
                                                          species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_SR  <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Strata removal",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth

sdm_TW_IID_depth_index_cod_SR  <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                         formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Strata removal",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG

sdm_DG_IID_index_cod_SR <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                  formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Strata removal",
                                                  species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#### IID + DG + depth

sdm_DG_IID_depth_index_cod_SR <- furrr::future_map2_dfr(sdm_data_cod_SR, mesh_sdm_cod_SR,
                                                        formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Strata removal",
                                                        species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############                           #############

# Scenario 4: AREA BLOCKED REDUCTION SCENARIO

#############                           #############

message("AREA BLOCKED REDUCTION SCENARIO")

################# Blocking the setdets

setdet_cod <- map(survey_cod, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod)){
  setdet_cod[[i]]$pop <- as.numeric(i)
}

# Getting all samples before year 10

samples10 <- map(setdet_cod, function(x) {subset(x, year < 11)}) ### samples before Year 11

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_cod, function(x) {subset(x, year >= 11)})### samples after Year 11

### Blocking the mpa coordinates

blocked_samples <- map(sample_a_y10, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_b30 <- survey_cod
for( i in seq_along(survey_b30)){
  survey_b30[[i]]$setdet <- combined_samples[[i]]}

################# Design-based index

design_index_b30 <- map_df(seq_along(survey_b30),function(i){
  strat <- SimSurvey::run_strat(survey_b30[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Area blocked") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_b30 <- map(survey_b30, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_cod_b30)){
  setdet_cod_b30[[i]]$pop <- as.numeric(i)
}

setdet_cod_b30 <- lapply(setdet_cod_b30, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod_b30, file = "./data/setdet_cod_b30.Rdata")

############# Bootstrapped index

boot_index_cod_b30 <- furrr::future_map_dfr(setdet_cod_b30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Area blocked")

################# sdmTMB

### Data prep

sdm_data_cod_b30 <- furrr::future_map(setdet_cod_b30, sdm_data_fn)

#save(sdm_data_cod_b30, file = "./data/sdm_data_cod_b30.Rdata")

mesh_sdm_cod_b30  <- purrr::map(sdm_data_cod_b30, mesh_sdm_fn, existing_mesh = inla_mesh)

### IID + NB2
sdm_NB2_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                       formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Area blocked",
                                                       species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                             formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Area blocked",
                                                             species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Area blocked",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                            formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Area blocked",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Area blocked",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                            formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Area blocked",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############                           #############

# Scenario 5: RECOVERY

#############                           #############

message("RECOVERY SCENARIO")

#############  Population simulations

set.seed(1)
survey_cod_rec <- furrr::future_map(seq_len(ITER), n_sims = 1, population_cod_recovery, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()
#save(survey_cod_rec, file = "./data/survey_cod_rec.Rdata")

############# True abundance

true_cod_rec <- map_df(seq_along(survey_cod_rec),function(i){
  tibble(year = unique(survey_cod_rec[[i]]$sp_N$year), true = as.numeric(colSums(survey_cod_rec[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Cod-like") |>
    group_by(pop)})

############# Data prep: Separating populations

setdet_cod_rec <- map(survey_cod_rec, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod_rec)){
  setdet_cod_rec[[i]]$pop  <- as.numeric(i)}

# Getting all samples before year 10

samples10_rec <- map(setdet_cod_rec, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10_rec <- map(setdet_cod_rec, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples_rec <- map(sample_a_y10_rec, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples_rec <- map(seq_along(blocked_samples_rec),function(i) {rbind(blocked_samples_rec[[i]], samples10_rec[[i]])})

################### Updating survey lists

survey_rec <- survey_cod_rec
for( i in seq_along(survey_rec)){
  survey_rec[[i]]$setdet  <- combined_samples_rec[[i]]}

################# Design-based index

design_index_rec <- map_df(seq_along(survey_rec),function(i){
  strat  <- SimSurvey::run_strat(survey_rec[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Recovery") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_rec <- map(survey_rec, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_cod_rec)){
  setdet_cod_rec[[i]]$pop  <- as.numeric(i)
}

setdet_cod_rec <- lapply(setdet_cod_rec, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod_rec, file = "./data/setdet_cod_rec.Rdata")

############# Bootstrapped index

boot_index_cod_rec <- furrr::future_map_dfr(setdet_cod_rec, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Recovery")
gc()

############# sdmTMB

### Data prep

sdm_data_cod_rec <- furrr::future_map(setdet_cod_rec, sdm_data_fn)

mesh_sdm_cod_rec <- purrr::map(sdm_data_cod_rec, mesh_sdm_fn, existing_mesh = inla_mesh)

sdm_newdata_cod <- sdm_newdata_fn(survey_cod_rec[[1]], sdm_data_cod_rec[[1]])


### Models

### IID + NB2
sdm_NB2_IID_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                    formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Recovery",
                                                    species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                          formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Recovery",
                                                          species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Recovery",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                         formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Recovery",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Recovery",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_rec <- furrr::future_map2_dfr(sdm_data_cod_rec, mesh_sdm_cod_rec,
                                                         formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Recovery",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)


#############                           #############

# Scenario 6: RECOVERY WITH SPILLOVER EFFECT

#############                           #############

message("RECOVERY WITH SPILLOVER EFFECT")

#############  Population simulations

set.seed(1)
survey_cod_so <- furrr::future_map(seq_len(ITER), n_sim = 1, population_cod_spillover, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

#save(survey_cod_so, file = "./data/survey_cod_so_05.Rdata")

############# True abundance

true_cod_so <- map_df(seq_along(survey_cod_so),function(i){
  tibble(year = unique(survey_cod_so[[i]]$sp_N$year), true = as.numeric(colSums(survey_cod_so[[i]]$I))) |>
    mutate(pop = as.numeric(i), species = "Cod-like") |>
    group_by(pop)})

############# Data prep: Separating populations

setdet_cod_so <- map(survey_cod_so, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod_so)){
  setdet_cod_so[[i]]$pop <- as.numeric(i)}

# Getting all samples before year 10

samples10_so <- map(setdet_cod_so, function(x) {subset(x, year < 11)})

# Removing the area from setdet, year >= 10

sample_a_y10_so <- map(setdet_cod_so, function(x) {subset(x, year >= 11)})

### Blocking the mpa coordinates

blocked_samples_so <- map(sample_a_y10_so, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples_so <- map(seq_along(blocked_samples_so),function(i) {rbind(blocked_samples_so[[i]], samples10_so[[i]])})

################### Updating survey lists

survey_so <- survey_cod_so
for( i in seq_along(survey_so)){
  survey_so[[i]]$setdet <- combined_samples_so[[i]]}

############# Design-based index

design_index_cod_so <- map_df(seq_along(survey_so),function(i){
  strat <- SimSurvey::run_strat(survey_so[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Recovery + Spillover") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations

setdet_cod_so <- map(survey_so, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod_so)){
  setdet_cod_so[[i]]$pop <- as.numeric(i)}

setdet_cod_so <- lapply(setdet_cod_so, function(x) split(x, x$sim)) |> flatten()

#save(setdet_cod_so, file = "./data/setdet_cod_so.Rdata")

############# Bootstrapped index

boot_index_cod_so <- furrr::future_map_dfr(setdet_cod_so, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Recovery + Spillover")
gc()

############# sdmTMB

### Data prep

sdm_data_cod_so <- furrr::future_map(setdet_cod_so, sdm_data_fn)

mesh_sdm_cod_so <- furrr::future_map(sdm_data_cod_so, mesh_sdm_fn, existing_mesh = inla_mesh)

sdm_newdata_cod <- sdm_newdata_fn(survey_so[[1]], sdm_data_cod_so[[1]])

### Models

### IID + NB2
sdm_NB2_IID_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                   formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "Recovery + Spillover",
                                                   species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                         formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "Recovery + Spillover",
                                                         species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                  formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "Recovery + Spillover",
                                                  species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                        formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "Recovery + Spillover",
                                                        species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                  formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "Recovery + Spillover",
                                                  species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_so <- furrr::future_map2_dfr(sdm_data_cod_so, mesh_sdm_cod_so,
                                                        formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "Recovery + Spillover",
                                                        species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############                           #############

# COMBINING ALL RESULTS

#############                           #############

index_cod_all_scenarios <- do.call(bind_rows, mget(ls(pattern = "index")))
index_cod_all_scenarios <- merge(index_cod_all_scenarios, true_cod, by=c("pop", "year", "species"))
save(index_cod_all_scenarios, file = "./data/index_cod_all_scenarios.Rdata")

Sys.time()
