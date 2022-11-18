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
plan(multisession)

### Load functions
source("./pop_cod_fn.R")
source("./model_run_fn.R")
source("./bootstrapping_fn.R")

#############               #############

            # BASE SCENARIO

#############               #############

#############  Population simulations

set.seed(1)
survey_cod <- furrr::future_map(seq_len(100), population_cod, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

save(survey_cod, file = "./Data/survey_cod_base.Rdata")

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
  setdet_cod[[i]]$pop <- as.numeric(i)
}

save(setdet_cod, file = "./Data/setdet_cod_base.Rdata")

############# Bootstrapped index

boot_index_cod <- furrr::future_map_dfr(setdet_cod, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Base", .progress = TRUE)
gc()

############# sdmTMB

### Data prep

sdm_data_fn <- function(x) {
  dat <- as_tibble(x) |>
    dplyr::select(x, y, set, sim, pop, year, count = n, tow_area, area = cell_area, depth) |>
    mutate(offset = log(tow_area), density = count / tow_area)
}

sdm_data_cod <- furrr::future_map(setdet_cod, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_fn <- function(sdm_data){
  mesh <- sdmTMB::make_mesh(sdm_data,
                            xy_cols = c("x", "y"),
                            cutoff = 45)
}

mesh_sdm_cod <- furrr::future_map(sdm_data_cod, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

sdm_newdata_fn <- function(survey, sdm_data){
  newdata <- as_tibble(dplyr::select(survey$grid_xy, x, y, depth)) |> distinct()
  newdata <- purrr::map_dfr(sort(unique(sdm_data$year)), ~ bind_cols(newdata, year = .)) |>
    mutate(offset = 0, area = sdm_data$area[1])
}

sdm_newdata_cod <- sdm_newdata_fn(survey_cod[[1]], sdm_data_cod[[1]]) ### since all populations has the same prediction area, newdata is same for all.

#save(sdm_newdata_cod, file = "./Data/sdm_newdata_cod.Rdata")

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
#############                           #############

# 30% SET DENSITY REDUCTION

#############                           #############

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
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "30% set reduction") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_r30 <- map(survey_r30, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod_r30)){
  setdet_cod_r30[[i]]$pop <- as.numeric(i)
}

############# Bootstrapped index

boot_index_cod_r30 <- furrr::future_map_dfr(setdet_cod_r30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "30% set reduction")

#################  sdmTMB

### Data prep

sdm_data_cod_r30 <- furrr::future_map(setdet_cod_r30, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_cod_r30  <- furrr::future_map(sdm_data_cod_r30, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

### IID + NB2
sdm_NB2_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                       formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "30% set reduction",
                                                       species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                             formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "30% set reduction",
                                                             species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "30% set reduction",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                            formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "30% set reduction",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "30% set reduction",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_r30 <- furrr::future_map2_dfr(sdm_data_cod_r30, mesh_sdm_cod_r30,
                                                            formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "30% set reduction",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############                           #############

# 30% AREA BLOCKED SCENARIO

#############                           #############

##  Survey scenarios

library(data.table)
library(raster)
library(sp)
library(sf)

set.seed(1)

# Making a grid with SimSurvey
grid <- make_grid(x_range = c(-150, 150),
                  y_range = c(-150, 150),
                  res = c(10, 10),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1000),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 40),
                  strat_splits = 2,
                  method = "spline")

# Converting from RasterBrick to SpatialPolygonsDataFrame

polylist_grid <- lapply(as.list(grid), rasterToPolygons)

poly_cell <- lapply(as.list(grid$cell), rasterToPolygons)

# Combining all attributes to a single SpatialPolygonsDataFrame and convert it to sf

grid_sf <- list(polylist_grid, makeUniqueIDs = T) |>
  flatten()

grid_sf <- do.call(cbind, grid_sf) |>
  st_as_sf()

# Creating a buffer area

grid_sf_buffered <- function(x, fraction) {
  width_half <- sqrt(st_area(st_union(x)) * fraction)/2
  buffer <- st_buffer(st_geometry(st_union(x)), -width_half)
  st_geometry(x)[!st_is_empty(buffer)] = buffer[!st_is_empty(buffer)]
  return(st_union(x))
}
grid_sf_buffered_perc30 <- grid_sf_buffered(grid_sf, 0.3)

### selecting a random point and creating the blocked area

set.seed(30)

block_poly_sf_fn <- function(x, fraction, buffer) {
  grid_centroid <- st_centroid(x)
  pointpool <- st_intersection(grid_centroid, buffer)

  selected <-  pointpool |>
    #group_by()|>
    slice_sample(n = 1) |>
    st_coordinates(selected)|>
    data.table()

  width_half <- sqrt(as.numeric(st_area(st_union(x))) * fraction)/2

  yplus <- selected$Y + width_half
  xplus <- selected$X + width_half
  yminus <- selected$Y - width_half
  xminus <- selected$X - width_half

  x1 <- c(xminus,xplus,xplus,xminus,xminus) # a(SW), b(SE), c(NE), d(NW), a(SW) - closing the vertices
  y1 <- c(yplus,yplus,yminus,yminus,yplus)

  block <- sp::Polygon(cbind(x1,y1))
  block <- sp::Polygons(list(block), ID = "1")
  block_poly <- sp::SpatialPolygons(list(block))
  block_poly_sf <- st_as_sf(block_poly)

  return(block_poly_sf)
}

block_poly_sf_30 <- block_poly_sf_fn(grid_sf, 0.3, grid_sf_buffered_perc30)

st_crs(grid_sf) = st_crs(block_poly_sf_30)

################# Blocking the setdets

setdet_cod <- map(survey_cod, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet_cod)){
  setdet_cod[[i]]$pop <- as.numeric(i)
}

# Converting the setdet to a SpatialPointsDataFrame and then to an sf object

for (i in 1:length(setdet_cod)){
  coordinates(setdet_cod[[i]]) = cbind(setdet_cod[[i]]$x, setdet_cod[[i]]$y)}

setdet_cod_sf <- map(setdet_cod, function(x) {st_as_sf(x)})

# Getting all samples before year 10

samples10 <- map(setdet_cod_sf, function(x) {subset(x, year < 11) |> st_drop_geometry() |> data.table()}) ### same for any perc

# Removing the area from setdet, year >= 10

sample_a_y10 <- map(setdet_cod_sf, function(x) {subset(x, year >= 11)})

### 30 percent setdet

blocked_samples <- map(seq_along(sample_a_y10), function(i) {
  blocked_samples <- st_difference(sample_a_y10[[i]], block_poly_sf_30)|>
    st_drop_geometry() |>
    data.table()})

combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

################### Updating survey lists

survey_b30 <- survey_cod
for( i in seq_along(survey_b30)){
  survey_b30[[i]]$setdet <- combined_samples_10[[i]]}

################# Design-based index

design_index_b30 <- map_df(seq_along(survey_b30),function(i){
  strat <- SimSurvey::run_strat(survey_b30[[i]])
  strat$total_strat |>
    mutate(N = total, upr = total_ucl, lwr = total_lcl, type = "Design-based", scenario = "30% area blocked") |>
    dplyr:: select(sim, year, N, type, lwr, upr, scenario) |>
    mutate(pop = as.numeric(i), species = "Cod-like")})
gc()

############# Data prep: Separating populations and iterations (survey simulations)

setdet_cod_b30 <- map(survey_b30, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_cod_b30)){
  setdet_cod_b30[[i]]$pop <- as.numeric(i)
}

############# Bootstrapped index

boot_index_cod_b30 <- furrr::future_map_dfr(setdet_cod_b30, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "30% area blocked")

################# sdmTMB

### Data prep

sdm_data_cod_b30 <- furrr::future_map(setdet_cod_b30, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_cod_b30  <- furrr::future_map(sdm_data_cod_b30, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

### IID + NB2
sdm_NB2_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                       formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "NB", scenario = "30% area blocked",
                                                       species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + NB2 + depth
sdm_NB2_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                             formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "NB + Depth", scenario = "30% area blocked",
                                                             species = "Cod-like", newdata = sdm_newdata_cod, model_run_NB, .id = "model", .progress = TRUE)
### IID + TW
sdm_TW_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "TW", scenario = "30% area blocked",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + TW + depth
sdm_TW_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                            formula = formula2, range_gt = 75, sigma_lt = 7.5, type = "TW + Depth", scenario = "30% area blocked",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_TW, .id = "model", .progress = TRUE)
### IID + DG
sdm_DG_IID_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                      formula = formula1, range_gt = 75, sigma_lt = 7.5, type = "DG", scenario = "30% area blocked",
                                                      species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)
#### IID + DG + depth
sdm_DG_IID_depth_index_cod_b30 <- furrr::future_map2_dfr(sdm_data_cod_b30, mesh_sdm_cod_b30,
                                                            formula = list(formula1, formula2), range_gt = 75, sigma_lt = 7.5, type = "DG + Depth", scenario = "30% area blocked",
                                                            species = "Cod-like", newdata = sdm_newdata_cod, model_run_DG, .id = "model", .progress = TRUE)

#############               #############

# STRATA REMOVAL SCENARIO

#############               #############

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

setdet_cod_SR <- map(survey_SR, function(x) {pluck(x, 'setdet')})

for(i in seq_along(setdet_cod_SR)){
  setdet_cod_SR[[i]]$pop <- as.numeric(i)
}

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

############# Bootstrapped index

boot_index_cod_SR <- furrr::future_map_dfr(setdet_cod_SR, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))|>
  mutate(species = "Cod-like", scenario = "Strata removal")

################# sdmTMB

### Data prep

sdm_data_cod_SR <- furrr::future_map(setdet_cod_SR, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))

mesh_sdm_cod_SR  <- furrr::future_map(sdm_data_cod_SR, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))

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

# COMBINING ALL RESULTS

#############                           #############

index_cod_all_scenarios <- do.call(bind_rows, mget(ls(pattern = "index")))
index_cod_all_scenarios <- merge(index_cod_all_scenarios, true_cod, by=c("pop", "year", "species"))

save(index_cod_all_scenarios, file = "./Data/index_cod_all_scenarios.Rdata")
