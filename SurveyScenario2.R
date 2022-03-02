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

### Samples only

setdet <- map(survey, function(x) {pluck(x, 'setdet')})

set.seed(2)

### Removing 20% of samples after year 10

#### Defining the new sample sets

perc20 <- function(x, year_cutoff) {
  a <- x %>%
    filter(.$year >= year_cutoff) %>%
    sample_frac(0.8)

  b <- x%>%
    filter(.$year < year_cutoff)

  ab <- rbind(a,b)
}

perc20 <- map(setdet, perc20, year_cutoff=10)


#### Calculating new index

survey20f <- function(i) {
  survey[[i]]$setdet <- perc20[[i]]
  survey_index <- survey[[i]] %>% run_strat()
}

survey20 <- map(seq_along(survey), survey20f)


### Removing 20% of samples after year 10

#### Defining the new sample sets

perc50 <- function(x, year_cutoff) {
  a <- x %>%
    filter(.$year >= year_cutoff) %>%
    sample_frac(0.5)

  b <- x%>%
    filter(.$year < year_cutoff)

  ab <- rbind(a,b)
}

perc50 <- map(setdet, perc50, year_cutoff=10)

#### Calculating new index

survey50f <- function(i) {
  survey[[i]]$setdet <- perc50[[i]]
  survey_index <- survey[[i]] %>% run_strat()
}

survey50 <- map(seq_along(survey), survey50f)
