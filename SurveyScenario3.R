library(dplyr)
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(purrr)
library(tictoc)

plan(multisession, workers = floor(availableCores()/2))

# Population
set.seed(1)

population <- function(iter) {
  set.seed(iter * 10)
  pop <- sim_abundance(ages = 1:20,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.63),
                                 log_sd = 0.3,
                                 phi_age = 0.9,
                                 phi_year = 0.5,
                                 plot = FALSE),
                       R = sim_R(log_mean = log(75e+06),
                                 log_sd = 0.5,
                                 random_walk = TRUE,
                                 plot = FALSE),
                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),
                       growth = sim_vonB(Linf = 120,
                                         L0 = 5,
                                         K = 0.11,
                                         log_sd = 0.15,
                                         plot = FALSE,
                                         length_group = 3))
}

# tic()
# purrr::map(seq_len(6), population)
# toc()

tic()
pop <- furrr::future_map(seq_len(6), population, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()

# Population distribution

distribution <- function(x) {
  sim_distribution(x,
                   grid = make_grid(x_range = c(-150, 150),
                                    y_range = c(-150, 150),
                                    res = c(10, 10),
                                    shelf_depth = 200,
                                    shelf_width = 100,
                                    depth_range = c(0, 1000),
                                    n_div = 1,
                                    strat_breaks = seq(0, 1000, by = 40),
                                    strat_splits = 2,
                                    method = "spline"),
                   ays_covar = sim_ays_covar(sd = 2.5,
                                             range = 200,
                                             phi_age = 0.5,
                                             phi_year = 0.9,
                                             group_ages = 12:20),
                   depth_par = sim_parabola(mu = log(80),
                                            sigma = 0.25, plot=FALSE, log_space = TRUE))
}

#a <- distribution(pop[[1]])

# tic()
# purrr::map(pop, distribution)
# toc()

tic()
sim <- furrr::future_map(pop, distribution, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()


# Survey

survey <-  function(x) {
  sim_survey(x,
             n_sims = 1,
             q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
             trawl_dim = c(1.5, 0.02),
             resample_cells = FALSE,
             binom_error = TRUE,
             min_sets = 2,
             set_den = 2/1000)
}

# tic()
# purrr::map(sim, survey)
# toc()

tic()
survey <- furrr::future_map(sim, survey, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()

library(sp)
library(sf)
library(data.table)
library(raster)


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

block_poly_sf <- block_poly_sf_fn(grid_sf, 0.1, grid_sf_buffered_perc10) ### block area is 10 % of the survey area

block_poly_sf2 <- block_poly_sf_fn(grid_sf, 0.2, grid_sf_buffered_perc20) ### block area is 20 % of the survey area


################# Blocking the setdets

setdet <- map(survey, function(x) {pluck(x, 'setdet')})

setdet <- map(setdet, function(x) {x %>%
    mutate(terrain = case_when(
      depth < 200 ~ "shallow",
      depth ==200 ~ "shelf",
      depth > 200 ~ "slope"
    ))})

# Converting the setdet to a SpatialPointsDataFrame and then to an sf object

for (i in 1:length(setdet)){
  coordinates(setdet[[i]]) = cbind(setdet[[i]]$x, setdet[[i]]$y)}

setdet_sf <- map(setdet, function(x) {st_as_sf(x)})

# Getting all samples before year 10

samples10 <- map (setdet_sf, function(x) {subset(x, year < 10) %>% st_drop_geometry() %>% data.table()}) ### same for any perc

# Removing the area from setdet, year >= 10

sample_a_y10 <- map (setdet_sf, function(x) {subset(x, year >= 10)})

### 10 percent setdet

blocked_samples <- map(seq_along(sample_a_y10), function(i) {
  blocked_samples <- st_difference(sample_a_y10[[i]], block_poly_sf) %>%
    st_drop_geometry() %>%
    data.table()
})

combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})


### 20 percent setdet

blocked_samples2 <- map(seq_along(sample_a_y10), function(i) {
  blocked_samples <- st_difference(sample_a_y10[[i]], block_poly_sf2) %>%
    st_drop_geometry() %>%
    data.table()
})

combined_samples_20 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples2[[i]], samples10[[i]])})


# plot(st_union(grid_sf), graticule = TRUE, axes = TRUE)
# plot(grid_sf_buffered_perc10, add=TRUE, col="yellow")
# plot(block_poly_sf, add= TRUE, col="green")
#
# plot(st_union(grid_sf), axes = TRUE)
# plot(grid_sf_buffered_perc20, add=TRUE, col="yellow")
# plot(block_poly_sf2, add= TRUE, col="green")


################# Design-based index

design_index_sc3_10  <- map(seq_along(survey), function(i){
  survey[[i]]$setdet <- combined_samples_10[[i]]
  survey_index <- survey[[i]] %>% run_strat()})

for( i in seq_along(design_index_sc3_10)){
  design_index_sc3_10[[i]]$iter <- as.numeric(i)
}

design_index_sc3_10_b <- as.data.frame(do.call(rbind, design_index_sc3_10))
design_index_sc3_10_b$scenario <- "3L"


design_index_sc3_20  <- map(seq_along(survey), function(i){
  survey[[i]]$setdet <- combined_samples_20[[i]]
  survey_index <- survey[[i]] %>% run_strat()})

for( i in seq_along(design_index_sc3_20)){
  design_index_sc3_20[[i]]$iter <- as.numeric(i)
}

design_index_sc3_20_b <- as.data.frame(do.call(rbind, design_index_sc3_20))
design_index_sc3_20_b$scenario <- "3H"


################# Bootstrapped

tic()
boot_index_sc3_10 <- furrr::future_map(combined_samples_10, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc3_10 )){
  boot_index_sc3_10 [[i]]$iter <- as.numeric(i)
}

boot_index_sc3_10_b <- as.data.frame(do.call(rbind, boot_index_sc3_10 ))
boot_index_sc3_10_b$scebario <- "3L"



tic()
boot_index_sc3_20 <- furrr::future_map(combined_samples_20, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc3_20 )){
  boot_index_sc3_20 [[i]]$iter <- as.numeric(i)
}

boot_index_sc3_20_b <- as.data.frame(do.call(rbind, boot_index_sc3_20 ))
boot_index_sc3_20_b$scebario <- "3H"

################### Updating survey lists

survey_sc3_10 <- survey
for( i in seq_along(survey_sc3_10)){
  survey_sc3_10[[i]]$setdet <- combined_samples_10[[i]]}

survey_sc3_20 <- survey
for( i in seq_along(survey_sc3_20)){
  survey_sc3_20[[i]]$setdet <- combined_samples_20[[i]]}


################# Model
library(data.table)

tic()
model_index_sc3_10 <- furrr::future_map(survey_sc3_10, sdm, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(model_index_sc3_10)){
  model_index_sc3_10[[i]]$iter <- as.numeric(i)
}
model_index_sc3_10_b <- as.data.frame(do.call(rbind, model_index_sc3_10))
model_index_sc3_10_b$scenario <- "3L"

tic()
model_index_sc3_20 <- furrr::future_map(survey_sc3_20, sdm, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(model_index_sc3_20)){
  model_index_sc3_20[[i]]$iter <- as.numeric(i)
}
model_index_sc3_20_b <- as.data.frame(do.call(rbind, model_index_sc3_20))
model_index_sc3_20_b$scenario <- "3H"


################# Ensamble

str(true_index2)
str(design_index_sc3_10_b)
str(design_index_sc3_20_b)
str(boot_index_sc3_10_b)
str(boot_index_sc3_20_b)
str(model_index_sc3_10_b)
str(model_index_sc3_20_b)

ensamb_sc3 <- bind_rows(true_index2,
                        design_index_sc3_10_b,
                        design_index_sc3_20_b,
                        boot_index_sc3_10_b,
                        boot_index_sc3_20_b,
                        model_index_sc3_10_b,
                        model_index_sc3_20_b)
