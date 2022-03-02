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

# Plucking samples

setdet <- map(survey, function(x) {pluck(x, 'setdet')})

# classifying the terrain group

setdet <- map(setdet, function(x) {x %>%
      mutate(terrain = case_when(
      depth < 200 ~ "Shallow",
      depth ==200 ~ "Shelf",
      depth > 200 ~ "Slope"
    ))})

# Converting the setdet to a SpatialPointsDataFrame and then to an sf object

for (i in 1:length(setdet)){
  coordinates(setdet[[i]]) = cbind(setdet[[i]]$x, setdet[[i]]$y)}

setdet_sf <- map(setdet, function(x) {st_as_sf(x)})

# Selecting a random point withing a terrain class

set.seed(5)
point_shelf <- spsample(setdet[[1]][which("terrain"=="Shelf")],n=1,"random")
point_shelf@coords

# Defining the width of the polygon

get_width <- function(total_area, fraction) {sqrt(total_area * fraction)}

width <- get_width(300*300, 0.1)
width20 <- get_width(300*300, 0.2)

# Define the edges

yplus <- point_shelf$coords.x1+width/2
xplus <- point_shelf$coords.x2+width/2
yminus <- point_shelf$coords.x1-width/2
xminus <- point_shelf$coords.x2-width/2

# Define the vertices

x1 <- c(xminus,xplus,xplus,xminus,xminus) # a(SW), b(SE), c(NE), d(NW), a(SW) - closing the vertices
y1 <- c(yplus,yplus,yminus,yminus,yplus)

# Assign the vertices to a polygon and then creating SpatialPolygons and an sf object

block_shelf <- sp::Polygon(cbind(x1,y1))
block_shelf <- sp::Polygons(list(block_shelf), ID = "shelf")
str(block_shelf,1)

block_shelfpoly <- sp::SpatialPolygons(list(block_shelf))
block_shelfpoly

block_shelfpoly_sf <- st_as_sf(block_shelfpoly)
block_shelfpoly_sf

# Plotting the block on the survey area

plot(setdet[[1]])
plot(block_shelfpoly_sf, add=TRUE, col = "red")

# Removing the area from setdet, year = 10

blocked_samples10 <- map (setdet_sf, function(x) {subset(x, year >= 10)})

blocked_samples <- map(seq_along(blocked_samples10), function(i) {
                  blocked_samples <-st_difference(blocked_samples10[[i]], block_shelfpoly_sf) %>%
                    st_drop_geometry() %>%
                    data.table()
                  })

# Getting all samples before year 10
samples10 <- map (setdet_sf, function(x) {subset(x, year < 10) %>% st_drop_geometry() %>% data.table()})

# combining all samples
combined_samples_10 <- map(seq_along(blocked_samples),function(i) {rbind(blocked_samples[[i]], samples10[[i]])})

# calculating index
block_10_index <- map(seq_along(survey), function(i){
  survey[[i]]$setdet <- combined_samples_10[[i]]
  survey_index <- survey[[i]] %>% run_strat()})



















