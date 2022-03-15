#############  Packages

# library(SimSurvey)
# library(sdmTMB)
# library(tidyr)
# library(future)
# library(purrr)
# library(tictoc)
# library(data.table)
# library(ggplot2)
# library(dplyr)
#
# plan(multisession, workers = floor(availableCores()/2))
#
# # Population
# set.seed(1)
#
# population <- function(iter) {
#   set.seed(iter * 10)
#   pop <- sim_abundance(ages = 1:20,
#                        years = 1:20,
#                        Z = sim_Z(log_mean = log(0.63),
#                                  log_sd = 0.3,
#                                  phi_age = 0.9,
#                                  phi_year = 0.5,
#                                  plot = FALSE),
#                        R = sim_R(log_mean = log(75e+06),
#                                  log_sd = 0.5,
#                                  random_walk = TRUE,
#                                  plot = FALSE),
#                        N0 = sim_N0(N0 = "exp",
#                                    plot = FALSE),
#                        growth = sim_vonB(Linf = 120,
#                                          L0 = 5,
#                                          K = 0.11,
#                                          log_sd = 0.15,
#                                          plot = FALSE,
#                                          length_group = 3))
# }
#
# # tic()
# # purrr::map(seq_len(6), population)
# # toc()
#
# tic()
# pop <- furrr::future_map(seq_len(6), population, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
# toc()
#
# # Population distribution
#
# distribution <- function(x) {
#   sim_distribution(x,
#                    grid = make_grid(x_range = c(-150, 150),
#                                     y_range = c(-150, 150),
#                                     res = c(10, 10),
#                                     shelf_depth = 200,
#                                     shelf_width = 100,
#                                     depth_range = c(0, 1000),
#                                     n_div = 1,
#                                     strat_breaks = seq(0, 1000, by = 40),
#                                     strat_splits = 2,
#                                     method = "spline"),
#                    ays_covar = sim_ays_covar(sd = 2.5,
#                                              range = 200,
#                                              phi_age = 0.5,
#                                              phi_year = 0.9,
#                                              group_ages = 12:20),
#                    depth_par = sim_parabola(mu = log(80),
#                                             sigma = 0.25, plot=FALSE, log_space = TRUE))
# }
#
# #a <- distribution(pop[[1]])
#
# # tic()
# # purrr::map(pop, distribution)
# # toc()
#
# tic()
# sim <- furrr::future_map(pop, distribution, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
# toc()

# Survey scenarios

library(sp)
library(sf)
library(data.table)
library(raster)

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
grid_sf <- list(polylist_grid, makeUniqueIDs = T) %>%
  flatten() %>%
  do.call(cbind, .) %>%
  st_as_sf() %>%
  mutate(terrain = case_when(
    depth < 200 ~ "shallow",
    depth == 200 ~ "shelf",
    depth > 200 ~ "slope"))

# Creating a buffer area

grid_sf_buffered <- function(x, fraction) {
  width_half <- sqrt(st_area(st_union(x)) * fraction)/2
  buffer <- st_buffer(st_geometry(st_union(x)), -width_half)
  st_geometry(x)[!st_is_empty(buffer)] = buffer[!st_is_empty(buffer)]
  return(st_union(x))
}
grid_sf_buffered_perc10 <- grid_sf_buffered(grid_sf, 0.1)
grid_sf_buffered_perc20 <- grid_sf_buffered(grid_sf, 0.2)


### selecting a random point and creating the blocked area

set.seed(30)

block_poly_sf_fn <- function(x, fraction, buffer) {
  grid_centroid <- st_centroid(x)
  pointpool <- st_intersection(grid_centroid, buffer)

  selected <-  pointpool %>%
    #group_by()%>%
    slice_sample(n = 1) %>%
    st_coordinates(selected)%>%
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

block_poly_sf_10 <- block_poly_sf_fn(grid_sf, 0.1, grid_sf_buffered_perc10) ### block area is 10 % of the survey area

block_poly_sf_20 <- block_poly_sf_fn(grid_sf, 0.2, grid_sf_buffered_perc20) ### block area is 20 % of the survey area

# 10 %

st_crs(grid_sf) = st_crs(block_poly_sf_10)

drop_cells_sc3_10 <-st_intersection(grid_sf, block_poly_sf_10) %>%
  dplyr::select(cell)%>%
  st_drop_geometry() %>%
  data.table()


set_rule_sc3_10_fn <- function(x) {

  standard_sets <- sim_sets(x, year < 10,  min_sets = 2, set_den = 2 / 1000)

  reduced_sets <- sim_sets(x, year >= 10 & !cell %in% drop_cells_sc3_10$cell,  min_sets = 2, set_den = 2/1000)

  sets <- rbind(standard_sets, reduced_sets)

  sets$set <- seq(nrow(sets)) # Important - make sure set has a unique ID.

  return(sets)
}

set_rule_sc3_10 <- furrr::future_map(sim, set_rule_sc3_10_fn, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))


survey_fn <- function(x, y) {survey <- SimSurvey::sim_survey(sim = x, custom_sets = y)}


survey_3M <- furrr::future_map2(sim, set_rule_sc3_10, survey_fn, .options = furrr::furrr_options(seed = TRUE))


## 20 %

st_crs(grid_sf) = st_crs(block_poly_sf_20)

drop_cells_sc3_20 <-st_intersection(grid_sf, block_poly_sf_20) %>%
  dplyr::select(cell)%>%
  st_drop_geometry() %>%
  data.table()


set_rule_sc3_20_fn <- function(x) {

  standard_sets <- sim_sets(x, year < 10,  min_sets = 2, set_den = 2 / 1000)

  reduced_sets <- sim_sets(x, year >= 10 & !cell %in% drop_cells_sc3_20$cell,  min_sets = 2, set_den = 2/1000)

  sets <- rbind(standard_sets, reduced_sets)

  sets$set <- seq(nrow(sets)) # Important - make sure set has a unique ID.

  return(sets)
}

set_rule_sc3_20 <- furrr::future_map(sim, set_rule_sc3_20_fn, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))


survey_fn <- function(x, y) {survey <- SimSurvey::sim_survey(sim = x, custom_sets = y)}


survey_3H <- furrr::future_map2(sim, set_rule_sc3_20, survey_fn, .options = furrr::furrr_options(seed = TRUE))


############# True abundance for graphs

true_index_3M <- map(seq_along(survey_3M),function(i){
  tibble(year = unique(survey_3M[[i]]$sp_N$year), N = as.numeric(colSums(survey_3M[[i]]$I)))
})

for( i in seq_along(true_index_3M)){
  true_index_3M[[i]]$iter <- as.numeric(i)
}

true_index_3M_df <- as.data.frame(do.call(rbind, true_index_3M))


################# Design-based index

design <- function(x) {
  SimSurvey::run_strat(x) }

### 10 %

tic()
design_sc3_10 <- furrr::future_map(survey_3M, design)
toc()

design_index_sc3_10 <- map(seq_along(design_sc3_10), function(i){
  design_sc3_10[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    dplyr:: select(year, N, type, lwr, upr)})

true_index_3M <- map(seq_along(survey_3M),function(i){
  tibble(year = unique(survey_3M[[i]]$sp_N$year), N = as.numeric(colSums(survey_3M[[i]]$I)))
})


for( i in seq_along(true_index_3M)){
  true_index_3M[[i]]$iter <- as.numeric(i)
  true_index_3M[[i]]$type <- "true"
}

for( i in seq_along(design_index_sc3_10)){
  design_index_sc3_10[[i]]$iter <- as.numeric(i)
  design_index_sc3_10[[i]]$true <- true_index_3M[[i]]$N
}

design_index_sc3_10_b <- as.data.frame(do.call(rbind, design_index_sc3_10))
design_index_sc3_10_b$scenario <- "3M"

### 20 %

tic()
design_sc3_20 <- furrr::future_map(survey_3H, design)
toc()


design_index_sc3_20 <- map(seq_along(design_sc3_20), function(i){
  design_sc3_20[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    dplyr:: select(year, N, type, lwr, upr)})

true_index_3H <- map(seq_along(survey_3H),function(i){
  tibble(year = unique(survey_3H[[i]]$sp_N$year), N = as.numeric(colSums(survey_3H[[i]]$I)))
})

for( i in seq_along(true_index_3H)){
  true_index_3H[[i]]$iter <- as.numeric(i)
  true_index_3H[[i]]$type <- "true"
}


for( i in seq_along(design_index_sc3_20)){
  design_index_sc3_20[[i]]$iter <- as.numeric(i)
  design_index_sc3_20[[i]]$true <- true_index_3H[[i]]$N
}

design_index_sc3_20_b <- as.data.frame(do.call(rbind, design_index_sc3_20))
design_index_sc3_20_b$scenario <- "3H"


################# Bootstrapped

sumYst <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>%
    ### stratum level
    group_by(year, strat, strat_area) %>%
    summarise(meanYh = mean(n), .groups = "drop_last") %>%
    mutate(Nh = strat_area/(1.5 * 0.02)) %>%
    group_by(year) %>%
    mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)%>%
    ### year level
    summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") %>%
    pull(sumYst)
}

boot_one_year <- function(x, reps) {
  b <- boot::boot(x, statistic = sumYst, strata = x$strat, R = reps)
  suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
  tibble::tibble(
    total = sumYst(x),
    mean_boot = mean(b$t),
    median_boot = median(b$t),
    lwr = bci$perc[[4]],
    upr = bci$perc[[5]],
    cv = sd(b$t) / mean(b$t),
    type = "Bootstrapped"
  )
}

boot_wrapper <- function(dat, reps) {
  out <- dat %>%
    split(dat$year) %>%
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}

### 10 %

setdet_3M <- map(survey_3M, function(x) {pluck(x, 'setdet')})

tic()
boot_index_sc3_10 <- furrr::future_map(setdet_3M, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc3_10 )){
  boot_index_sc3_10 [[i]]$iter <- as.numeric(i)
  boot_index_sc3_10[[i]]$true <- true_index_3M[[i]]$N
}

boot_index_sc3_10_b <- as.data.frame(do.call(rbind, boot_index_sc3_10 ))
boot_index_sc3_10_b$scenario <- "3M"

boot_index_sc3_10_b <- boot_index_sc3_10_b %>%
  dplyr:: select(year, mean_boot, lwr, upr, type, iter, scenario, true) %>%
  rename(N=mean_boot)

### 20 %

setdet_3H <- map(survey_3H, function(x) {pluck(x, 'setdet')})

tic()
boot_index_sc3_20 <- furrr::future_map(setdet_3H, boot_wrapper, reps=1000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc3_20 )){
  boot_index_sc3_20 [[i]]$iter <- as.numeric(i)
  boot_index_sc3_20[[i]]$true <- true_index_3H[[i]]$N
}

boot_index_sc3_20_b <- as.data.frame(do.call(rbind, boot_index_sc3_20 ))
boot_index_sc3_20_b$scenario <- "3H"

boot_index_sc3_20_b <- boot_index_sc3_20_b %>%
  dplyr:: select(year, mean_boot, lwr, upr, type, iter, scenario, true) %>%
  rename(N=mean_boot)

################# sdmTMB

### Functions
sdm_data_fn <- function(x) {

  xy <- as_tibble(x$grid_xy)
  dat <- as_tibble(x$setdet) %>%
    dplyr::select(x, y, set, year, N = n, tow_area) %>%
    left_join(., xy, by = c("x", "y")) %>%
    mutate(offset = log(tow_area))
}

mesh_sdm_fn <- function(sdm_data){
  mesh <- sdmTMB::make_mesh(
    sdm_data,
    xy_cols = c("x", "y"),
    cutoff = 40)
}

sdm_IID_fn <- function(x, y){

  fit_IID <- sdmTMB(N ~ 0 + as.factor(year) + offset,
                    data = x,
                    mesh = y,
                    time = "year",
                    family = nbinom2(link = "log"),
                    spatial = TRUE,
                    spatiotemporal = "IID")
}

sdm_AR1_fn <- function(x, y) {

  fit_AR1 <- sdmTMB(N ~ 0 + as.factor(year) + offset,
                    data = x,
                    mesh = y,
                    time = "year",
                    family = nbinom2(link = "log"),
                    spatial = TRUE,
                    spatiotemporal = c("AR1"))
}

sdm_newdata_fn <- function(survey, sdm_data) {

  grid_dat <- as_tibble(dplyr::select(survey$grid_xy, x, y, depth)) %>% distinct()
  grid_dat <- purrr::map_dfr(sort(unique(sdm_data$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(sdm_data$offset)
  grid_dat$area <-survey$setdet$cell_area[1]/survey$setdet$tow_area[1]
  return(grid_dat)
}

sdm_prediction_IID_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE,
                      area = y$area)
}

sdm_prediction_AR1_fn <- function(x, y){
  pred_AR1 <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE,
                      area = y$area)
}

sdm_index_IID_fn <- function(x){
  index_IID <- get_index(x) %>%
    mutate(type = "IID", N = est)
}

sdm_index_AR1_fn <- function(x){
  index_AR1 <- get_index(x) %>%
    mutate(type = "AR1", N = est)
}


### 10 %

tic()
sdm_data_sc3_10 <- furrr::future_map(survey_3M, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
mesh_sdm_sc3_10 <- furrr::future_map(sdm_data_sc3_10, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
sdm_IID_sc3_10 <- furrr::future_map2(sdm_data_sc3_10, mesh_sdm_sc3_10, sdm_IID_fn)
toc()

tic()
sdm_AR1_sc3_10 <- furrr::future_map2(sdm_data_sc3_10, mesh_sdm_sc3_10, sdm_AR1_fn)
toc()

tic()
sdm_newdata_sc3_10 <- furrr::future_map2(survey_3M, sdm_data_sc3_10, sdm_newdata_fn)
toc()

tic()
sdm_prediction_IID_sc3_10 <- furrr::future_map2(sdm_IID_sc3_10, sdm_newdata_sc3_10, sdm_prediction_IID_fn)
toc()

tic()
sdm_prediction_AR1_sc3_10 <- furrr::future_map2(sdm_AR1_sc3_10, sdm_newdata_sc3_10, sdm_prediction_AR1_fn)
toc()

tic()
sdm_index_IID_sc3_10 <- furrr::future_map(sdm_prediction_IID_sc3_10, sdm_index_IID_fn)
toc()

for( i in seq_along(sdm_index_IID_sc3_10)){
  sdm_index_IID_sc3_10[[i]]$iter <- as.numeric(i)
  sdm_index_IID_sc3_10[[i]]$true <- true_index_3M[[i]]$N
}
sdm_index_IID_2_sc3_10 <- as.data.frame(do.call(rbind, sdm_index_IID_sc3_10))
sdm_index_IID_2_sc3_10$scenario <- "3M"

tic()
sdm_index_AR1_sc3_10 <- furrr::future_map(sdm_prediction_AR1_sc3_10, sdm_index_AR1_fn)
toc()

for( i in seq_along(sdm_index_AR1_sc3_10)){
  sdm_index_AR1_sc3_10[[i]]$iter <- as.numeric(i)
  sdm_index_AR1_sc3_10[[i]]$true <- true_index_3M[[i]]$N
}
sdm_index_AR1_2_sc3_10 <- as.data.frame(do.call(rbind, sdm_index_AR1_sc3_10))
sdm_index_AR1_2_sc3_10$scenario <- "3M"

### 20 %

tic()
sdm_data_sc3_20 <- furrr::future_map(survey_3H, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
mesh_sdm_sc3_20 <- furrr::future_map(sdm_data_sc3_20, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
sdm_IID_sc3_20 <- furrr::future_map2(sdm_data_sc3_20, mesh_sdm_sc3_20, sdm_IID_fn)
toc()

tic()
sdm_AR1_sc3_20 <- furrr::future_map2(sdm_data_sc3_20, mesh_sdm_sc3_20, sdm_AR1_fn)
toc()

tic()
sdm_newdata_sc3_20 <- furrr::future_map2(survey_3H, sdm_data_sc3_20, sdm_newdata_fn)
toc()

tic()
sdm_prediction_IID_sc3_20 <- furrr::future_map2(sdm_IID_sc3_20, sdm_newdata_sc3_20, sdm_prediction_IID_fn)
toc()

tic()
sdm_prediction_AR1_sc3_20 <- furrr::future_map2(sdm_AR1_sc3_20, sdm_newdata_sc3_20, sdm_prediction_AR1_fn)
toc()

tic()
sdm_index_IID_sc3_20 <- furrr::future_map(sdm_prediction_IID_sc3_20, sdm_index_IID_fn)
toc()

for( i in seq_along(sdm_index_IID_sc3_20)){
  sdm_index_IID_sc3_20[[i]]$iter <- as.numeric(i)
  sdm_index_IID_sc3_20[[i]]$true <- true_index_3H[[i]]$N
}
sdm_index_IID_2_sc3_20 <- as.data.frame(do.call(rbind, sdm_index_IID_sc3_20))
sdm_index_IID_2_sc3_20$scenario <- "3H"

tic()
sdm_index_AR1_sc3_20 <- furrr::future_map(sdm_prediction_AR1_sc3_20, sdm_index_AR1_fn)
toc()

for( i in seq_along(sdm_index_AR1_sc3_20)){
  sdm_index_AR1_sc3_20[[i]]$iter <- as.numeric(i)
  sdm_index_AR1_sc3_20[[i]]$true <- true_index_3H[[i]]$N
}
sdm_index_AR1_2_sc3_20 <- as.data.frame(do.call(rbind, sdm_index_AR1_sc3_20))
sdm_index_AR1_2_sc3_20$scenario <- "3H"

################# Combining results

str(true_index_df)

str(design_index_sc3_10_b)
str(design_index_sc3_20_b)
str(boot_index_sc3_10_b)
str(boot_index_sc3_20_b)
str(sdm_index_IID_2_sc3_10)
str(sdm_index_IID_2_sc3_20)
str(sdm_index_AR1_2_sc3_10)
str(sdm_index_AR1_2_sc3_20)

results_scenario3 <- bind_rows(design_index_sc3_10_b,
                        design_index_sc3_20_b,
                        boot_index_sc3_10_b,
                        boot_index_sc3_20_b,
                        sdm_index_IID_2_sc3_10,
                        sdm_index_IID_2_sc3_20,
                        sdm_index_AR1_2_sc3_10,
                        sdm_index_AR1_2_sc3_20)

### Medium intensity

results_scenario3 %>%
  filter(scenario == "3M")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~factor(type, levels=c('Design-based','Bootstrapped','IID','AR1')),  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_index_3M_df, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  theme(legend.position = "none")

### High intensity

results_scenario3 %>%
  filter(scenario == "3H")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~factor(type, levels=c('Design-based','Bootstrapped','IID','AR1')),  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_index_3M_df, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  theme(legend.position = "none")
