
#############  Packages

library(dplyr)
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(purrr)
library(tictoc)
library(ggplot2)

plan(multisession, workers = floor(availableCores()/2))

#############  Population

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

############# Population distribution

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

# tic()
# purrr::map(pop, distribution)
# toc()

tic()
sim <- furrr::future_map(pop, distribution, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()

############# Survey

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

############# Design-based index

design <- function(x) {
  run_strat(x) }

# tic()
# purrr::map(survey, design)
# toc()

tic()
design <- furrr::future_map(survey, design, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()

design_index <- map(seq_along(design), function(i){
  design[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    select(year, N, type, lwr, upr)})

for( i in seq_along(design_index)){
  design_index[[i]]$iter <- as.numeric(i)
}

design_index2 <- as.data.frame(do.call(rbind, design_index))
design_index2$scenario <- "base"

############# Bootstrapped index

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

setdet <- map(survey, function(x) {pluck(x, 'setdet')})

# tic()
# purrr::map(setdet, boot_wrapper, reps=100)
# toc()

tic()
boot_index <- furrr::future_map(setdet, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
toc()


for( i in seq_along(boot_index)){
  boot_index[[i]]$iter <- as.numeric(i)
}

boot_index2 <- as.data.frame(do.call(rbind, boot_index))
boot_index2$scenario <- "base"

boot_index2 <- boot_index2 %>% select(year, mean_boot, lwr, upr, type, iter, scenario) %>%
  rename(N  = mean_boot)


############# sdmTMB

sdm <- function(x) {

xy <- as_tibble(x$grid_xy)
dat <- as_tibble(x$setdet) %>%
  dplyr::select(x, y, set, year, N = n, tow_area) %>%
  left_join(., xy, by = c("x", "y")) %>%
  mutate(offset = log(tow_area))

mesh <- sdmTMB::make_mesh(
  dat,
  xy_cols = c("x", "y"),
  cutoff = 20)

fit_IID <- sdmTMB(N ~ 0 + as.factor(year) + offset,
              data = dat,
              mesh = mesh,
              time = "year",
              family = nbinom2(link = "log"),
              spatial = TRUE,
              spatiotemporal = "IID")

fit_AR1 <- sdmTMB(N ~ 0 + as.factor(year) + offset,
               data = dat,
               mesh = mesh,
               time = "year",
               family = nbinom2(link = "log"),
               spatial = TRUE,
               spatiotemporal = c("AR1"))

grid_dat <- as_tibble(dplyr::select(x$grid_xy, x, y, depth)) %>% distinct()
grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
grid_dat$offset <- mean(dat$offset)
grid_dat$area <-x$setdet$cell_area[1]/x$setdet$tow_area[1]

pred_IID <- predict(fit_IID,
                newdata = grid_dat,
                return_tmb_object = TRUE,
                area = grid_dat$area)

pred_AR1 <- predict(fit_AR1,
                newdata = grid_dat,
                return_tmb_object = TRUE,
                area = grid_dat$area)


index_IID <- get_index(pred_IID) %>%
mutate(type = "IID", N = est)

index_AR1 <- get_index(pred_AR1) %>%
  mutate(type = "AR1", N = est)

index_sdm <- index_IID %>%
  bind_rows(index_AR1)

return(index_sdm)

}

# tic()
# purrr::map(survey, sdm)
# toc()

tic()
model_index <- furrr::future_map(survey, sdm, .options = furrr::furrr_options(seed = TRUE))
toc()

for( i in seq_along(model_index)){
  model_index[[i]]$iter <- as.numeric(i)
}
model_index2 <- as.data.frame(do.call(rbind, model_index))
model_index2$scenario <- "base"

############# True abundance

true_index <- map(seq_along(survey), function(i){
  survey[[i]]$setdet %>%
    group_by(year) %>%
    summarise(N = sum(N)) %>%
    mutate(type= "True")})

for( i in seq_along(true_index)){
  true_index[[i]]$iter <- as.numeric(i)
}

true_index2 <- as.data.frame(do.call(rbind, true_index))

true_index2_no_type$type <- NULL

############# Combining results

str(true_index2)
str(design_index2)
str(boot_index2)
str(model_index2)

result_base <- bind_rows(true_index2, design_index2, boot_index2, model_index2)

result_base %>%
  filter(type!="True")%>%
ggplot(aes(year, N, group = type)) +
  geom_line(aes(colour = type), size=1) +
  facet_grid(iter~type,  scales = "free_y")+
  geom_line(aes(year, N), size=1, data= true_index2_no_type, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3)+
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  theme(legend.title= "Method")
