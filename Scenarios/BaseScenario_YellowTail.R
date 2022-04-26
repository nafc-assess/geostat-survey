#############  Packages
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(purrr)
library(tictoc)
library(ggplot2)
library(data.table)
library(dplyr)


plan(multisession, workers = floor(availableCores()/2))
#plan(multisession, workers = 2)
#options(future.globals.onReference = "error")

#sessionInfo()


#############  Population

set.seed(1)

population <- function(iter) {
  set.seed(iter * 17)
  pop <- sim_abundance(ages = 1:10,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.64),
                                 log_sd = 0.1,
                                 phi_age = 0.2,
                                 phi_year = 0.4,
                                 plot = FALSE),
                       R = sim_R(log_mean = log(12e+09),
                                 log_sd = 0.6,
                                 random_walk = TRUE,
                                 plot = FALSE),
                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),
                       growth = sim_vonB(Linf = 56,
                                         L0 = 0,
                                         K = 0.13,
                                         log_sd = 0.12,
                                         plot = FALSE,
                                         length_group = 2))
}

# tic()
# purrr::map(seq_len(6), population)
# toc()

tic()
pop <- furrr::future_map(seq_len(12), population, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
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
                   ays_covar = sim_ays_covar(sd = 1.4,
                                             range = 300,
                                             phi_age = 0.8,
                                             phi_year = 0.7,
                                             group_ages = 5:9),
                   depth_par = sim_parabola(mu = log(100),
                                            sigma = 0.15, plot=FALSE, log_space = TRUE))
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
             n_sims = 5,
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

true_index <- map(seq_along(survey),function(i){
  tibble(year = unique(survey[[i]]$sp_N$year), N = as.numeric(colSums(survey[[i]]$I)))
})

for( i in seq_along(true_index)){
  true_index[[i]]$pop <- as.numeric(i)
}

true_index_yellowtail <- as.data.frame(do.call(rbind, true_index)) |> rename(true = N)
true_index_yellowtail$type <- "True index"


#save(true_index_yellowtail, file = "./Yellowtail/true_index_yellowtail.Rdata")


############# Design-based index


design <- function(x) {
  SimSurvey::run_strat(x) }

# tic()
# purrr::map(survey, design)
# toc()

tic()
design <- furrr::future_map(survey, design, .options = furrr::furrr_options(seed = TRUE))
toc()

# library(parallel)
# detectCores()
# #plan(multicore)
# tic()
# design <- parallel::mclapply(survey, design)
# toc()


design_index <- map(seq_along(design), function(i){
  design[[i]]$total_strat %>%
    mutate(N = total, upr =total_ucl, lwr = total_lcl, type = "Design-based") %>%
    dplyr:: select(sim, year, N, type, lwr, upr)})

for( i in seq_along(design_index)){
  design_index[[i]]$pop <- as.numeric(i)
  design_index[[i]]$species <- "yellowtail"
}

design_index <- as.data.frame(do.call(rbind, design_index))
design_index$scenario <- "base"


ggplot(design_index, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Design-based index")


### Separating populations and iterations (survey simulations)


setdet <- map(survey, function(x) {pluck(x, 'setdet')})

for( i in seq_along(setdet)){
  setdet[[i]]$pop <- as.numeric(i)
}

setdet_sim <- NULL
for( i in seq_along(setdet)){
  setdet_sim[[i]] <- setdet[[i]] |> group_split(sim)
}

setdet_sim <- unlist(setdet_sim,recursive=FALSE)

############# Bootstrapped index

# sumYst <- function(dat, i = seq_len(nrow(dat))) {
#   dat[i, ] %>%
#     ### stratum level
#     group_by(year, strat, strat_area) %>%
#     summarise(meanYh = mean(n), .groups = "drop_last") %>%
#     mutate(Nh = strat_area/(1.5 * 0.02)) %>%
#     group_by(year) %>%
#     mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)%>%
#     ### year level
#     summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") %>%
#     pull(sumYst)
# }
#
# boot_one_year <- function(x, reps) {
#   b <- boot::boot(x, statistic = sumYst, strata = x$strat, R = reps)
#   suppressWarnings(bci <- boot::boot.ci(b, type = "bca"))
#   tibble::tibble(
#     pop = mean(x$pop),
#     sim = mean(x$sim),
#     total = sumYst(x),
#     mean_boot = mean(b$t),
#     median_boot = median(b$t),
#     lwr = bci$bca[[4]],
#     upr = bci$bca[[5]],
#     cv = sd(b$t) / mean(b$t),
#     type = "Bootstrapped"
#   )
# }
#
# boot_wrapper <- function(dat, reps) {
#   out <- dat %>%
#   split(dat$year) %>%
#   purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
#   out$year <- as.numeric(out$year)
#   out
# }
#
# tic()
# boot_index <- furrr::future_map(setdet_sim, boot_wrapper, reps=2000, .options = furrr::furrr_options(seed = TRUE))
# toc()
#
#
# for( i in seq_along(boot_index)){
#   #boot_index[[i]]$true <- true_index[[i]]$N
#   boot_index[[i]]$species <- "yellowtail"
# }
#
# boot_index2 <- as.data.frame(do.call(rbind, boot_index))
# boot_index2$scenario <- "base"
# boot_index2$N <- boot_index2$mean_boot
#
#
# ggplot(boot_index2, aes(year, mean_boot, group = sim)) +
#   facet_wrap(~pop, scales = "free_y")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
#   geom_line(aes(colour = type), size=1) +
#   geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
#   labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
#   #facet_wrap(~iter)+
#   theme_minimal()+
#   scale_y_log10() +
#   theme(legend.position = "none") +
#   ggtitle("Bootstrapped index")

############# sdmTMB


sdm_data_fn <- function(x) {

    dat <- as_tibble(x) %>%
    dplyr::select(x, y, set, sim, pop, year, count = n, tow_area, area = cell_area, depth) |>
    mutate(offset = log(tow_area), density = count / tow_area)
}


tic()
sdm_data <- furrr::future_map(setdet_sim, sdm_data_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

tic()
sdm_data <- parallel::mclapply(setdet_sim, sdm_data_fn)
toc()

#save(sdm_data, file = "sdm_data.Rdata")

mesh_sdm_fn <- function(sdm_data){
  mesh <- sdmTMB::make_mesh(sdm_data,
    xy_cols = c("x", "y"),
    cutoff = 40)
}

tic()
mesh_sdm <- furrr::future_map(sdm_data, mesh_sdm_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

# tic()
# mesh_sdm <- parallel::mclapply(sdm_data, mesh_sdm_fn)
# toc()

#save(mesh_sdm, file = "mesh_sdm.Rdata")

# for prediction

sdm_newdata <- as_tibble(dplyr::select(survey[[1]]$grid_xy, x, y, depth)) %>% distinct()
sdm_newdata <- purrr::map_dfr(sort(unique(sdm_data[[1]]$year)), ~ bind_cols(sdm_newdata, year = .))
sdm_newdata$offset <- 0
sdm_newdata$area <-sdm_data[[1]]$area[1]

#save(sdm_newdata, file = "sdm_newdata.Rdata")

### IID + NB2

sdm_IID_fn <- function(x, y){

  fit_IID <- sdmTMB(count ~ 0 + as.factor(year),
                    offset = x$offset,
                    data = x,
                    mesh = y,
                    time = "year",
                    family = nbinom2(link = "log"),
                    spatial = TRUE,
                    spatiotemporal = "IID",
                    #priors = sdmTMBpriors(
                      #matern_s = pc_matern(range_gt = 40, sigma_lt = 5),
                      #matern_st = pc_matern(range_gt = 40, sigma_lt = 5)),
                    #share_range = FALSE
                    )
}

tic()
sdm_IID <- furrr::future_map2(sdm_data, mesh_sdm, sdm_IID_fn)
toc()

# tic()
# sdm_IID <- parallel::mcmapply(sdm_data, mesh_sdm, FUN=sdm_IID_fn)
# toc()


sdm_prediction_IID_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID <- furrr::future_map(sdm_IID, y=sdm_newdata, sdm_prediction_IID_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_fn <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + NB2", N = est)
}

tic()
sdm_index_IID <- furrr::future_map(sdm_prediction_IID, sdm_index_IID_fn)
toc()

for( i in seq_along(sdm_index_IID)){
  sdm_index_IID[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID[[i]]$species <- "yellowtail"
}
sdm_index_IID_2 <- as.data.frame(do.call(rbind, sdm_index_IID))
sdm_index_IID_2$scenario <- "base"

ggplot(sdm_index_IID_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - IID + NB2")

### IID + NB2 + depth

sdm_IID_depth_fn <- function(x, y){

  fit_IID_depth <- sdmTMB(count ~  0 + as.factor(year) + s(log(depth), k= 4),
                    data = x,
                    offset = x$offset,
                    mesh = y,
                    time = "year",
                    family = nbinom2(link = "log"),
                    spatial = TRUE,
                    spatiotemporal = "IID")
}

tic()
sdm_IID_depth <- furrr::future_map2(sdm_data, mesh_sdm, sdm_IID_depth_fn)
toc()

sdm_prediction_IID_depth_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID_depth <- furrr::future_map(sdm_IID_depth, y=sdm_newdata, sdm_prediction_IID_depth_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_depth_fn <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + NB2 + depth", N = est)
}

tic()
sdm_index_IID_depth <- furrr::future_map(sdm_prediction_IID_depth, sdm_index_IID_depth_fn)
toc()

for( i in seq_along(sdm_index_IID_depth)){
  sdm_index_IID_depth[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID_depth[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID[[i]]$species <- "yellowtail"
}
sdm_index_IID_depth_2 <- as.data.frame(do.call(rbind, sdm_index_IID_depth))
sdm_index_IID_depth_2$scenario <- "base"

ggplot(sdm_index_IID_depth_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - IID + NB2 + depth")

### IID + TW

sdm_IID_tw_fn <- function(x, y){

  fit_IID_tw <- sdmTMB(density ~ 0 + as.factor(year),
                    data = x,
                    mesh = y,
                    time = "year",
                    family = tweedie(),
                    spatial = TRUE,
                    spatiotemporal = "IID",
                    #priors = sdmTMBpriors(
                      #matern_s = pc_matern(range_gt = 40, sigma_lt = 5),
                      #matern_st = pc_matern(range_gt = 40, sigma_lt = 5)),
                    #share_range = FALSE
                    )
}

tic()
sdm_IID_tw <- furrr::future_map2(sdm_data, mesh_sdm, sdm_IID_tw_fn)
toc()

sdm_prediction_IID_tw_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID_tw <- furrr::future_map(sdm_IID_tw, y=sdm_newdata, sdm_prediction_IID_tw_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_tw_fn <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + TW", N = est)
}

tic()
sdm_index_IID_tw <- furrr::future_map(sdm_prediction_IID_tw, sdm_index_IID_tw_fn)
toc()

for( i in seq_along(sdm_index_IID_tw)){
  sdm_index_IID_tw[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID_tw[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID_tw[[i]]$species <- "yellowtail"
}
sdm_index_IID_tw_2 <- as.data.frame(do.call(rbind, sdm_index_IID_tw))
sdm_index_IID_tw_2$scenario <- "base"

ggplot(sdm_index_IID_tw_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - Tweedie + IID")


### IID + TW + depth

fit_IID_tw_depth <- function(x, y){

  fit_IID_tw_depth <- sdmTMB(density ~ 0 + as.factor(year) + s(log(depth), k= 4),
                       data = x,
                       mesh = y,
                       time = "year",
                       family = tweedie(),
                       spatial = TRUE,
                       spatiotemporal = "IID")
}

tic()
sdm_IID_tw_depth <- furrr::future_map2(sdm_data, mesh_sdm, fit_IID_tw_depth)
toc()

sdm_prediction_IID_tw_fn_depth <- function(x, y){

  pred_IID_depth <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID_tw_depth <- furrr::future_map(sdm_IID_tw_depth, y=sdm_newdata, sdm_prediction_IID_tw_fn_depth, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_tw_fn_depth <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + TW + depth", N = est)
}

tic()
sdm_index_IID_tw_depth <- furrr::future_map(sdm_prediction_IID_tw_depth, sdm_index_IID_tw_fn_depth)
toc()

for( i in seq_along(sdm_index_IID_tw_depth)){
  sdm_index_IID_tw_depth[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID_tw_depth[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID_tw_depth[[i]]$species <- "yellowtail"
}
sdm_index_IID_tw_depth_2 <- as.data.frame(do.call(rbind, sdm_index_IID_tw_depth))
sdm_index_IID_tw_depth_2$scenario <- "base"


ggplot(sdm_index_IID_tw_depth_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - Tweedie + IID + depth")


### IID + DG

sdm_IID_dg_fn <- function(x, y){

  fit_IID_dg <- sdmTMB(density ~ 0 + as.factor(year),
                       data = x,
                       mesh = y,
                       #offset = x$offset,
                       time = "year",
                       family = delta_gamma(),
                       spatial = TRUE,
                       spatiotemporal = "IID",
                       #priors = sdmTMBpriors(
                         #matern_s = pc_matern(range_gt = 40, sigma_lt = 5),
                         #matern_st = pc_matern(range_gt = 40, sigma_lt = 5)),
                       #share_range = FALSE
                       )
}

tic()
sdm_IID_dg <- furrr::future_map2(sdm_data, mesh_sdm, sdm_IID_dg_fn)
toc()

#save(sdm_IID_dg, file = "sdm_IID_dg.Rdata")

sdm_prediction_IID_dg_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID_dg <- furrr::future_map(sdm_IID_dg, y=sdm_newdata, sdm_prediction_IID_dg_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_dg_fn <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + DG", N = est)
}

tic()
sdm_index_IID_dg <- furrr::future_map(sdm_prediction_IID_dg, sdm_index_IID_dg_fn)
toc()

for( i in seq_along(sdm_index_IID_dg)){
  sdm_index_IID_dg[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID_dg[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID_dg[[i]]$species <- "yellowtail"
}
sdm_index_IID_dg_2 <- as.data.frame(do.call(rbind, sdm_index_IID_dg))
sdm_index_IID_dg_2$scenario <- "base"

ggplot(sdm_index_IID_dg_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - IID + DG")

#### IID + DG + depth

sdm_IID_dg_depth_fn <- function(x, y){

  fit_IID_dg_depth <- sdmTMB(density ~ 0 + as.factor(year) + s(log(depth), k= 4),
                       data = x,
                       #offset = x$offset,
                       mesh = y,
                       time = "year",
                       family = delta_gamma(),
                       spatial = TRUE,
                       spatiotemporal = "IID",
                       #priors = sdmTMBpriors(
                         #matern_s = pc_matern(range_gt = 40, sigma_lt = 5),
                         #matern_st = pc_matern(range_gt = 40, sigma_lt = 5)),
                       #share_range = FALSE
                       )
}

tic()
sdm_IID_dg_depth <- furrr::future_map2(sdm_data, mesh_sdm, sdm_IID_dg_depth_fn)
toc()

sdm_prediction_IID_dg_depth_fn <- function(x, y){

  pred_IID <- predict(x,
                      newdata = y,
                      return_tmb_object = TRUE)
}

tic()
sdm_prediction_IID_dg_depth <- furrr::future_map(sdm_IID_dg_depth, y=sdm_newdata, sdm_prediction_IID_dg_depth_fn, .options = furrr::furrr_options(seed = TRUE))
toc()

sdm_index_IID_dg_depth_fn <- function(x){
  index_IID <- get_index(x, area = x$data$area, bias_correct = TRUE) %>%
    mutate(type = "IID + DG + depth", N = est)
}

tic()
sdm_index_IID_dg_depth <- furrr::future_map(sdm_prediction_IID_dg_depth, sdm_index_IID_dg_depth_fn)
toc()

for( i in seq_along(sdm_index_IID_dg_depth)){
  sdm_index_IID_dg_depth[[i]]$sim <- unique(sdm_data[[i]]$sim)
  sdm_index_IID_dg_depth[[i]]$pop <- unique(sdm_data[[i]]$pop)
  sdm_index_IID_dg_depth[[i]]$species <- "yellowtail"
}
sdm_index_IID_dg_depth_2 <- as.data.frame(do.call(rbind, sdm_index_IID_dg_depth))
sdm_index_IID_dg_depth_2$scenario <- "base"


ggplot(sdm_index_IID_dg_depth_2, aes(year, N, group = sim)) +
  facet_wrap(~pop, scales = "free_y")+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+
  geom_line(aes(colour = type), size=1) +
  geom_line(aes(year, true), size=1, data= true_index_yellowtail, inherit.aes = FALSE) +
  labs(x = "Year", y = "N", colour = "type", fill = "type", size = "type")+
  #facet_wrap(~iter)+
  theme_minimal()+
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle("Model-based index - IID + DG + depth")


############### index

index_yellowtail <- bind_rows(
  design_index,
  #boot_index2,
  sdm_index_IID_2,
  sdm_index_IID_depth_2,
  sdm_index_IID_tw_2,
  sdm_index_IID_tw_depth,
  sdm_index_IID_dg_2,
  sdm_index_IID_dg_depth_2)


save(index_yellowtail, file = "index_yellowtail.Rdata")


index_all_models_gm <- index_all_models|>
  group_by(pop, sim) |>
  mutate(N_gm = N / exp(mean(log(N)))) |>
  as.data.frame()

true_index_yellowtail_gm <- true_index_yellowtail |>
  group_by(pop, sim) |>
  mutate(N_gm = true / exp(mean(log(true))))


index_all_models <- merge(index_all_models_gm, true_index_yellowtail_gm, by=c("year", "pop"))

index_all_models$error <- index_all_models$N_gm.x - index_all_models$N_gm.y
index_all_models$Rerror <- (index_all_models$N_gm.x - index_all_models$N_gm.y) / index_all_models$N_gm.y
