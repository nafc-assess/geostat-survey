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

design_index_sc2_20 <- map(seq_along(survey), survey20f)


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

#### Calculating new index

survey50f <- function(i) {
  survey[[i]]$setdet <- perc50[[i]]
  survey_index <- survey[[i]] %>% run_strat()
}

design_index_sc2_50 <- map(seq_along(survey), survey50f)

#### Bootstrapped

tic()
boot_index_sc2_20 <- furrr::future_map(perc20, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
boot_index_sc2_50 <- furrr::future_map(perc50, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()

#################  Model

survey_sc2_20 <- survey
for( i in seq_along(survey_sc2_20)){
  survey_sc2_20[[i]]$setdet <- perc20[[i]]}

tic()
model_index_sc2_20 <- furrr::future_map(survey_sc2_20, sdm, .options = furrr::furrr_options(seed = TRUE))
toc()

survey_sc2_50 <- survey
for( i in seq_along(survey_sc2_50)){
  survey_sc2_50[[i]]$setdet <- perc20[[i]]}

tic()
model_index_sc2_50 <- furrr::future_map(survey_sc2_50, sdm, .options = furrr::furrr_options(seed = TRUE))
toc()


#################  True abundance

true_index <- map(seq_along(survey), function(i){
  survey[[i]]$setdet %>%
    group_by(year) %>%
    summarise(N = sum(N)) %>%
    mutate(type= "True")})

for( i in seq_along(true_index)){
  true_index[[i]]$iter <- as.numeric(i)
}

true_index2 <- as.data.frame(do.call(rbind, true_index))
true_index2$scenario <- "true"

#################  Design based

design_index_sc2_20 <- map(seq_along(design_index_sc2_20), function(i){
  design_index_sc2_20[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    select(year, N, type, lwr, upr)})

for( i in seq_along(design_index_sc2_20)){
  design_index_sc2_20[[i]]$iter <- as.numeric(i)
}

design_index_sc2_20_b <- as.data.frame(do.call(rbind, design_index_sc2_20))
design_index_sc2_20_b$scenario <- "2L"


design_index_sc2_50 <- map(seq_along(design_index_sc2_50), function(i){
  design_index_sc2_50[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    select(year, N, type, lwr, upr)})

for( i in seq_along(design_index_sc2_50)){
  design_index_sc2_50[[i]]$iter <- as.numeric(i)
}

design_index_sc2_50_b <- as.data.frame(do.call(rbind, design_index_sc2_50))
design_index_sc2_50_b$scenario <- "2H"

#################  bootstapped

for( i in seq_along(boot_index_sc2_20 )){
  boot_index_sc2_20 [[i]]$iter <- as.numeric(i)
}

boot_index_sc2_20_b <- as.data.frame(do.call(rbind, boot_index_sc2_20 ))
boot_index_sc2_20_b$scenario <- "2L"

boot_index_sc2_20_b <- boot_index_sc2_20_b %>% select(year, mean_boot, lwr, upr, type, iter, scenario) %>%
  rename(N=mean_boot)


for( i in seq_along(boot_index_sc2_50 )){
  boot_index_sc2_50 [[i]]$iter <- as.numeric(i)
}

boot_index_sc2_50_b <- as.data.frame(do.call(rbind, boot_index_sc2_50 ))
boot_index_sc2_50_b$scenario <- "2H"

boot_index_sc2_50_b <- boot_index_sc2_50_b %>% select(year, mean_boot, lwr, upr, type, iter, scenario) %>%
  rename(N  = mean_boot)


################# model index

for( i in seq_along(model_index_sc2_20)){
  model_index_sc2_20[[i]]$iter <- as.numeric(i)
}
model_index_sc2_20_b <- as.data.frame(do.call(rbind, model_index_sc2_20))
model_index_sc2_20_b$scenario <- "2L"

for( i in seq_along(model_index_sc2_50)){
  model_index_sc2_50[[i]]$iter <- as.numeric(i)
}
model_index_sc2_50_b <- as.data.frame(do.call(rbind, model_index_sc2_50))
model_index_sc2_50_b$scenario <- "2H"

#############
str(true_index2)
str(design_index_sc2_20_b)
str(design_index_sc2_50_b)
str(boot_index_sc2_20_b)
str(boot_index_sc2_50_b)
str(model_index_sc2_20_b)
str(model_index_sc2_50_b)

ensamb_sc2 <- bind_rows(true_index2, design_index_sc2_20_b, design_index_sc2_50_b, boot_index_sc2_20_b, boot_index_sc2_50_b,
                        model_index_sc2_20_b, model_index_sc2_50_b)


df_2L <- ensamb_sc2 %>%
  filter(scenario == "2L"  | scenario == "true")
df_2L


df_2L %>%
  filter(type!="True")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~type,  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_db, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()


df_2H <- ensamb_sc2 %>%
  filter(scenario == "2H"  | scenario == "true")
df_2H


df_2H %>%
  filter(type!="True")%>%
  ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~type,  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_db, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()
