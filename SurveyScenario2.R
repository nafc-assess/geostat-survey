library(dplyr)
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(purrr)
library(tictoc)
library(ggplot2)

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

############# True abundance

true_index <- map(seq_along(survey), function(i){
  survey[[i]]$setdet %>%
    group_by(year) %>%
    summarise(N = sum(N)) %>%
    mutate(type= "true")})

for( i in seq_along(true_index)){
  true_index[[i]]$iter <- as.numeric(i)
}

true_index2 <- as.data.frame(do.call(rbind, true_index))

true_index2$type <- NULL

############# Scenarios

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


### Removing 50% of samples after year 10

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

############# Design-based index

#### 20 %

survey20f <- function(i) {
  survey[[i]]$setdet <- perc20[[i]]
  survey_index <- survey[[i]] %>% run_strat()
}

design_index_sc2_20 <- map(seq_along(survey), survey20f)

design_index_sc2_20 <- map(seq_along(design_index_sc2_20), function(i){
  design_index_sc2_20[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    select(year, N, type, lwr, upr)})

for( i in seq_along(design_index_sc2_20)){
  design_index_sc2_20[[i]]$iter <- as.numeric(i)
}

design_index_sc2_20_b <- as.data.frame(do.call(rbind, design_index_sc2_20))
design_index_sc2_20_b$scenario <- "2M"


#### 50 %

survey50f <- function(i) {
  survey[[i]]$setdet <- perc50[[i]]
  survey_index <- survey[[i]] %>% run_strat()
}

design_index_sc2_50 <- map(seq_along(survey), survey50f)

design_index_sc2_50 <- map(seq_along(design_index_sc2_50), function(i){
  design_index_sc2_50[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    select(year, N, type, lwr, upr)})

for( i in seq_along(design_index_sc2_50)){
  design_index_sc2_50[[i]]$iter <- as.numeric(i)
}

design_index_sc2_50_b <- as.data.frame(do.call(rbind, design_index_sc2_50))
design_index_sc2_50_b$scenario <- "2H"

#### Bootstrapped

#### 20 %
tic()
boot_index_sc2_20 <- furrr::future_map(perc20, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc2_20 )){
  boot_index_sc2_20 [[i]]$iter <- as.numeric(i)
}

boot_index_sc2_20_b <- as.data.frame(do.call(rbind, boot_index_sc2_20 ))
boot_index_sc2_20_b$scenario <- "2M"

boot_index_sc2_20_b <- boot_index_sc2_20_b %>% select(year, mean_boot, lwr, upr, type, iter, scenario) %>%
  rename(N=mean_boot)

#### 50 %
tic()
boot_index_sc2_50 <- furrr::future_map(perc50, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(boot_index_sc2_50 )){
  boot_index_sc2_50 [[i]]$iter <- as.numeric(i)
}

boot_index_sc2_50_b <- as.data.frame(do.call(rbind, boot_index_sc2_50 ))
boot_index_sc2_50_b$scenario <- "2H"

boot_index_sc2_50_b <- boot_index_sc2_50_b %>% select(year, mean_boot, lwr, upr, type, iter, scenario) %>%
  rename(N  = mean_boot)

#################  sdmTMB

#### 20 %
survey_sc2_20 <- survey
for( i in seq_along(survey_sc2_20)){
  survey_sc2_20[[i]]$setdet <- perc20[[i]]}

tic()
sdm_data_sc2_20 <- furrr::future_map(survey_sc2_20, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
mesh_sdm_sc2_20 <- furrr::future_map(sdm_data_sc2_20, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
sdm_IID_sc2_20 <- furrr::future_map2(sdm_data_sc2_20, mesh_sdm_sc2_20, sdm_IID_fn)
toc()

tic()
sdm_AR1_sc2_20 <- furrr::future_map2(sdm_data_sc2_20, mesh_sdm_sc2_20, sdm_AR1_fn)
toc()

tic()
sdm_newdata_sc2_20 <- furrr::future_map2(survey, sdm_data_sc2_20, sdm_newdata_fn)
toc()

tic()
sdm_prediction_IID_sc2_20 <- furrr::future_map2(sdm_IID_sc2_20, sdm_newdata_sc2_20, sdm_prediction_IID_fn)
toc()

tic()
sdm_prediction_AR1_sc2_20 <- furrr::future_map2(sdm_AR1_sc2_20, sdm_newdata_sc2_20, sdm_prediction_AR1_fn)
toc()

tic()
sdm_index_IID_sc2_20 <- furrr::future_map(sdm_prediction_IID_sc2_20, sdm_index_IID_fn)
toc()

for( i in seq_along(sdm_index_IID_sc2_20)){
  sdm_index_IID_sc2_20[[i]]$iter <- as.numeric(i)
  sdm_index_IID_sc2_20[[i]]$true <- true_index[[i]]$N
}
sdm_index_IID_2_sc2_20 <- as.data.frame(do.call(rbind, sdm_index_IID_sc2_20))
sdm_index_IID_2_sc2_20$scenario <- "2M"

tic()
sdm_index_AR1_sc2_20 <- furrr::future_map(sdm_prediction_AR1_sc2_20, sdm_index_AR1_fn)
toc()

for( i in seq_along(sdm_index_AR1_sc2_20)){
  sdm_index_AR1_sc2_20[[i]]$iter <- as.numeric(i)
  sdm_index_AR1_sc2_20[[i]]$true <- true_index[[i]]$N
}
sdm_index_AR1_2_sc2_20 <- as.data.frame(do.call(rbind, sdm_index_AR1_sc2_20))
sdm_index_AR1_2_sc2_20$scenario <- "2M"

#### 50 %

survey_sc2_50 <- survey
for( i in seq_along(survey_sc2_50)){
  survey_sc2_50[[i]]$setdet <- perc50[[i]]}

tic()
sdm_data_sc2_50 <- furrr::future_map(survey_sc2_50, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
mesh_sdm_sc2_50 <- furrr::future_map(sdm_data_sc2_50, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
sdm_IID_sc2_50 <- furrr::future_map2(sdm_data_sc2_50, mesh_sdm_sc2_50, sdm_IID_fn)
toc()

tic()
sdm_AR1_sc2_50 <- furrr::future_map2(sdm_data_sc2_50, mesh_sdm_sc2_50, sdm_AR1_fn)
toc()

tic()
sdm_newdata_sc2_50 <- furrr::future_map2(survey, sdm_data_sc2_50, sdm_newdata_fn)
toc()

tic()
sdm_prediction_IID_sc2_50 <- furrr::future_map2(sdm_IID_sc2_50, sdm_newdata_sc2_50, sdm_prediction_IID_fn)
toc()

tic()
sdm_prediction_AR1_sc2_50 <- furrr::future_map2(sdm_AR1_sc2_50, sdm_newdata_sc2_50, sdm_prediction_AR1_fn)
toc()

tic()
sdm_index_IID_sc2_50 <- furrr::future_map(sdm_prediction_IID_sc2_50, sdm_index_IID_fn)
toc()

for( i in seq_along(sdm_index_IID_sc2_50)){
  sdm_index_IID_sc2_50[[i]]$iter <- as.numeric(i)
  sdm_index_IID_sc2_50[[i]]$true <- true_index[[i]]$N
}
sdm_index_IID_2_sc2_50 <- as.data.frame(do.call(rbind, sdm_index_IID_sc2_50))
sdm_index_IID_2_sc2_50$scenario <- "2H"

tic()
sdm_index_AR1_sc2_50 <- furrr::future_map(sdm_prediction_AR1_sc2_50, sdm_index_AR1_fn)
toc()

for( i in seq_along(sdm_index_AR1_sc2_50)){
  sdm_index_AR1_sc2_50[[i]]$iter <- as.numeric(i)
  sdm_index_AR1_sc2_50[[i]]$true <- true_index[[i]]$N
}
sdm_index_AR1_2_sc2_50 <- as.data.frame(do.call(rbind, sdm_index_AR1_sc2_50))
sdm_index_AR1_2_sc2_50$scenario <- "2H"



#############
str(true_index2)

str(design_index_sc2_20_b)
str(design_index_sc2_50_b)
str(boot_index_sc2_20_b)
str(boot_index_sc2_50_b)
str(sdm_index_IID_2_sc2_20)
str(sdm_index_IID_2_sc2_50)
str(sdm_index_AR1_2_sc2_20)
str(sdm_index_AR1_2_sc2_50)

result_scenario2 <- bind_rows(design_index_sc2_20_b, design_index_sc2_50_b, boot_index_sc2_20_b, boot_index_sc2_50_b,
                              sdm_index_IID_2_sc2_20, sdm_index_IID_2_sc2_50,sdm_index_AR1_2_sc2_20, sdm_index_AR1_2_sc2_50)

### Medium intensity

result_scenario2 %>%
  filter(scenario == "2M")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~factor(type, levels=c('Design-based','Bootstrapped','IID','AR1')),  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_index2, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  theme(legend.position = "none")


 ### Medium intensity

result_scenario2 %>%
  filter(scenario == "2H")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~factor(type, levels=c('Design-based','Bootstrapped','IID','AR1')),  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_index2, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  theme(legend.position = "none")
