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

#############               #############

            # BASE SCENARIO

#############               #############


#############  Population simulations

set.seed(1)
survey_yellowtail <- furrr::future_map(seq_len(100), n_sims = 1, population_yellowtail, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

#save(survey_yellowtail, file = "./Data/survey_yellowtail_base.Rdata")

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

#save(setdet_yellowtail, file = "./Data/setdet_yellowtail_base.Rdata")

############# Bootstrapped index

# boot_index_yellowtail <- furrr::future_map_dfr(setdet_yellowtail, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
#   mutate(species = "Yellowtail-like", scenario = "Base", .progress = TRUE)
# gc()

############# sdmTMB

### Data prep

sdm_data_fn <- function(x) {
  dat <- as_tibble(x) |>
    dplyr::select(x, y, set, sim, pop, year, count = n, tow_area, area = cell_area, depth) |>
    mutate(offset = log(tow_area), density = count / tow_area)
}

sdm_data_yellowtail <- furrr::future_map(setdet_yellowtail, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_fn <- function(sdm_data){
  mesh <- sdmTMB::make_mesh(sdm_data,
                            xy_cols = c("x", "y"),
                            cutoff = 45)
}

mesh_sdm_yellowtail <- furrr::future_map(sdm_data_yellowtail, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

sdm_newdata_fn <- function(survey, sdm_data){
  newdata <- as_tibble(dplyr::select(survey$grid_xy, x, y, depth)) |> distinct()
  newdata <- purrr::map_dfr(sort(unique(sdm_data$year)), ~ bind_cols(newdata, year = .)) |>
    mutate(offset = 0, area = sdm_data$area[1])
}

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

# boot_index_yellowtail_r30 <- furrr::future_map_dfr(setdet_yellowtail_r30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
#   mutate(species = "Yellowtail-like", scenario = "Set reduction")

#################  sdmTMB

### Data prep

sdm_data_yellowtail_r30 <- furrr::future_map(setdet_yellowtail_r30, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_yellowtail_r30  <- furrr::future_map(sdm_data_yellowtail_r30, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

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

#############                           #############

            # AREA BLOCKED REDUCTION SCENARIO

#############                           #############

################# Blocking the setdets

setdet_yellowtail <- map(survey_yellowtail, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_yellowtail)){
  setdet_yellowtail[[i]]$pop <- as.numeric(i)
}

# Converting the setdet to a SpatialPointsDataFrame and then to an sf object

for (i in 1:length(setdet_yellowtail)){
  coordinates(setdet_yellowtail[[i]]) = cbind(setdet_yellowtail[[i]]$x, setdet_yellowtail[[i]]$y)}

setdet_yellowtail_sf <- map(setdet_yellowtail, function(x) {st_as_sf(x)})

# Getting all samples before year 10

samples10 <- map(setdet_yellowtail_sf, function(x) {subset(x, year < 11) |> st_drop_geometry() |> data.table()}) ### same for any perc

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_yellowtail_sf, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples <- map(sample_a_y10, function(x) {subset(x, !(x < 50 & x > -100 & y > -70 & y < 110))})

combined_samples <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_b30 <- survey_yellowtail
for( i in seq_along(survey_b30)){
  survey_b30[[i]]$setdet <- combined_samples_10[[i]]}

################# Design-based index

design_index_b30 <- map_df(seq_along(survey_b30),function(i){
  strat <- SimSurvey::run_strat(survey_b30[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "Area blocked reduction") |>
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
  mutate(species = "Yellowtail-like", scenario = "Area blocked reduction")

################# sdmTMB

### Data prep

sdm_data_yellowtail_b30 <- furrr::future_map(setdet_yellowtail_b30, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_yellowtail_b30  <- furrr::future_map(sdm_data_yellowtail_b30, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

### IID + NB2
sdm_NB2_IID_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                       formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "NB", scenario = "Area blocked reduction",
                                                       species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                             formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "NB + Depth", scenario = "Area blocked reduction",
                                                             species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "TW", scenario = "Area blocked reduction",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                            formula = formula2, range_gt = 125, sigma_lt = 12.5, type = "TW + Depth", scenario = "Area blocked reduction",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                      formula = formula1, range_gt = 125, sigma_lt = 12.5, type = "DG", scenario = "Area blocked reduction",
                                                      species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_yellowtail_b30 <- furrr::future_map2_dfr(sdm_data_yellowtail_b30, mesh_sdm_yellowtail_b30,
                                                            formula = list(formula1, formula2), range_gt = 125, sigma_lt = 12.5, type = "DG + Depth", scenario = "Area blocked reduction",
                                                            species = "Yellowtail-like", newdata = sdm_newdata_yellowtail, model_run_DG, .id = "model", .progress = TRUE)
#############               #############

        # STRATA REMOVAL SCENARIO

#############               #############

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

# boot_index_yellowtail_SR <- furrr::future_map_dfr(setdet_yellowtail_SR, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
#   mutate(species = "Yellowtail-like", scenario = "Strata removal")

################# sdmTMB

### Data prep

sdm_data_yellowtail_SR <- furrr::future_map(setdet_yellowtail_SR, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_yellowtail_SR  <- furrr::future_map(sdm_data_yellowtail_SR, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))


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


#############               #############

# COMBINING ALL RESULTS

#############               #############

index_yellowtail_all_scenarios <- do.call(bind_rows, mget(ls(pattern = "index")))
index_yellowtail_all_scenarios <- merge(index_yellowtail_all_scenarios, true_yellowtail, by=c("pop", "year", "species"))
