# Experiment with various fits to SimSurvey output
# Fix residuals?
# 2022-01-24 SA

library(SimSurvey)
library(ggplot2)
theme_set(theme_minimal())
library(dplyr)
library(sdmTMB)

set.seed(1)
pop <- sim_abundance(
  ages = 1:20, # fish age vector
  years = 1:10, # year vector
  Z = sim_Z(
    log_mean = log(0.6), # Total mortality matrix -- rnorm
    log_sd = 0.1,
    phi_age = 0.9,
    phi_year = 0.6,
    plot = FALSE
  ),
  R = sim_R(
    log_mean = log(3e+07), # Recruitment vector -- rnorm
    log_sd = 0.5,
    random_walk = TRUE,
    plot = FALSE
  ),
  N0 = sim_N0(
    N0 = "exp", # Starting abundance vector -- exponential Z; N0 at age 1 is R0
    plot = FALSE
  ),
  growth = sim_vonB(
    Linf = 100, # abundance-at-age matrix
    L0 = 5,
    K = 0.2,
    log_sd = 0.05,
    plot = FALSE,
    length_group = 3
  )
) %>%
  sim_distribution(
    grid = make_grid(
      x_range = c(-150, 150), # km
      y_range = c(-150, 150), # km
      res = c(10, 10), # km
      shelf_depth = 200, # Approximate depth of the shelf in m
      shelf_width = 50, # Approximate width of the shelf in km
      depth_range = c(10, 1000), # m
      n_div = 1, # number of division
      strat_breaks = seq(0, 1000, by = 40), # strata depth breaks, meter
      strat_splits = 2, # horizontally split strat
      method = "spline"
    ), # spline, loess, linear interpolation options
    ays_covar = sim_ays_covar(
      range = 300, # adding noise correlated across space, year, and age dimensions -- covariance matrix, higher value higher spatial correlation
      phi_age = 0.8, # strong correlation across ages
      phi_year = 0.1
    ), # low correlation across years
    depth_par = sim_parabola(
      mu = 200, # defining relationship between abundance and depth, i.e., adding a parabolic depth 'preference', fish prefers occurring around 200 meters
      sigma = 70, plot = FALSE
    )
  )

survey <- sim_survey(
  pop, # sim_distribution
  n_sims = 1, # number of survey
  q = sim_logistic(k = 2, x0 = 3, plot = FALSE), # function of catchability at age (k= The steepness of the curve, x0= x-value of the sigmoid's midpoint, curve max value=1)
  # trawl_dim = c(1.5, 0.02), #Trawl width and distance (default=same units as grid)
  # resample_cells = FALSE, #Allow resampling of sampling units?
  binom_error = TRUE, # Impose binomial error?
  # min_sets = 2, #Minimum number of sets per strat
  set_den = 2 / 1000, # set per strata ~ strata area with a set density rule, sample size = sum(set numbers per strat) * years
  # lengths_cap = 500, #Maximum number of lengths measured per set
  # ages_cap = 10,
  # age_sampling = "stratified",
  # age_length_group = 1,
  # age_space_group = "division",
  # light = TRUE
) %>%
  SimSurvey::run_strat() |>
  SimSurvey::strat_error()

xy <- as_tibble(survey$grid_xy)
dat <- as_tibble(survey$setdet) %>%
  select(x, y, set, year, N = n, tow_area) %>%
  left_join(., xy, by = c("x", "y")) %>%
  mutate(offset = log(tow_area)) # setting up an "offset" variable to represent the tow area

mesh <- sdmTMB::make_mesh(
  dat,
  xy_cols = c("x", "y"),
  cutoff = 15
)

ctl <- sdmTMBcontrol(start = list(
  ln_tau_O = 2.5, ln_tau_E = 2,
  ln_kappa = rep(-4, 2)
))

m0 <- sdmTMB(N ~ 0 + as.factor(year) + offset,
  data = dat,
  mesh = mesh,
  time = "year",
  family = nbinom2(),
  control = ctl,
  silent = FALSE
)
simulate(m0, nsim = 300) |> dharma_residuals(m0)
simulate(m0, re_form = NA, nsim = 300) |> dharma_residuals(m0)

# mcmc residuals ----------------------------------------------------

fitted_ln_kappa <- unname(m0$model$par[names(m0$model$par) == "ln_kappa"])
m0_fixed <- sdmTMB(N ~ 0 + as.factor(year) + offset,
  data = dat,
  mesh = mesh,
  time = "year",
  family = nbinom2(),
  control = sdmTMBcontrol(
    start = list(ln_kappa = rep(fitted_ln_kappa, 2)), # fix kappa for speed
    map = list(ln_kappa = rep(factor(NA), 2))
  ),
  silent = FALSE
)
# warning: a bit slow!
stan_fit <- tmbstan::tmbstan(m0_fixed$tmb_obj,
  iter = 104L,
  warmup = 100L, chains = 1L
)
stan_eta <- predict(m0_fixed, tmbstan_model = stan_fit)
stan_mu <- m0_fixed$family$linkinv(stan_eta)
quant_res <- function(object, y, mu) {
  pars <- object$model$par
  if ("b_disp_k" %in% names(pars)) {
    phi <- exp(object$model$par[["b_disp_k"]]) # if on dispformula branch 2022-01-24
  } else {
    phi <- exp(object$model$par[["ln_phi"]]) # if on main branch 2022-01-24
  }
  a <- pnbinom(y - 1, size = phi, mu = mu)
  b <- pnbinom(y, size = phi, mu = mu)
  u <- runif(n = length(y), min = a, max = b)
  qnorm(u)
}
par(mfrow = c(2, 2), cex = 0.8)
for (i in seq_len(4)) {
  mcmc_res <- quant_res(m0_fixed, y = dat$N, mu = stan_mu[, i])
  qqnorm(mcmc_res)
  qqline(mcmc_res)
}

# delta-gamma -------------------------------------------------------
dat$density <- dat$N / dat$tow_area
dat$present <- ifelse(dat$density > 0, 1L, 0L)

m1 <- sdmTMB(present ~ 0 + as.factor(year),
  data = dat,
  mesh = mesh,
  time = "year",
  family = binomial(),
  control = ctl,
  silent = FALSE
)
simulate(m1, nsim = 300) |> dharma_residuals(m1)

dat_present <- filter(dat, present == 1L)
mesh_present <- mesh <- sdmTMB::make_mesh(
  dat_present,
  xy_cols = c("x", "y"),
  cutoff = 15,
  mesh = mesh$mesh
)
m2 <- sdmTMB(density ~ 0 + as.factor(year),
  data = dat_present,
  mesh = mesh,
  time = "year",
  family = Gamma(link = "log"),
  control = ctl,
  silent = FALSE
)
simulate(m2, nsim = 300) |> dharma_residuals(m2)

p2 <- predict(m2, newdata = NULL)
ggplot(p2, aes(x, y, colour = exp(est))) +
  geom_point() +
  scale_color_viridis_c(trans = "log10") +
  facet_wrap(~year)
ggplot(dat_present, aes(x, y, colour = density)) +
  geom_point() +
  scale_color_viridis_c(trans = "log10") +
  facet_wrap(~year)
p2$residual <- log(dat_present$density) - p2$est
ggplot(p2, aes(x, y, colour = residual)) +
  geom_point() +
  scale_color_gradient2() +
  facet_wrap(~year)

# stan_fit <- tmbstan::tmbstan(m2$tmb_obj, iter = 100, warmup = 99, chains = 1)
# r <- residuals(m2)

# m3 <- sdmTMB(density ~ 0 + as.factor(year),
#   data = dat_present,
#   mesh = mesh,
#   time = "year",
#   control = ctl,
#   family = lognormal(link = "log"),
#   silent = FALSE
# )
# simulate(m3, nsim = 300) |> dharma_residuals(m3)
#
# dat_present$depth_scaled <- as.numeric(scale(dat_present$depth))
# m4 <- sdmTMB(
#   density ~ 0 + as.factor(year),
#   dispformula = ~ 0 + as.factor(year) + depth_scaled,
#   data = dat_present,
#   mesh = mesh,
#   time = "year",
#   control = ctl,
#   family = lognormal(link = "log"),
#   silent = FALSE
# )
# m4
# simulate(m4, nsim = 300) |> dharma_residuals(m4)

# Get index and plot ------------------------------------------------

xy_sim <- as_tibble(pop$grid_xy)
df_sim <- as_tibble(pop$sp_N)
df_sim <- left_join(df_sim, xy_sim, by = "cell")
grid_dat <- tidyr::expand_grid(x = sort(unique(df_sim$x)), y = sort(unique(df_sim$y)))
grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
grid_dat$offset <- mean(dat$offset)

tow_area <- survey$setdet$tow_area[1] # cell and tow area is consistent in the simulation
cell_area <- survey$setdet$cell_area[1]
pred_pos <- predict(m2, newdata = grid_dat, sims = 200L)
pred_bin <- predict(m1, newdata = grid_dat, sims = 200L)
pred_tot <- log(exp(pred_pos) * plogis(pred_bin))

ind_delta_gamma <- get_index_sims(pred_tot, area = rep(1, nrow(grid_dat)))
ggplot(ind_delta, aes(year, ymin = lwr, ymax = upr, y = est)) +
  geom_ribbon() +
  geom_line()

pred_nb2 <- predict(m0, newdata = grid_dat, sims = 200L)
ind_nbinom2 <- get_index_sims(pred_nb2, area = rep(cell_area / tow_area, nrow(grid_dat)))
ggplot(ind_nbinom2, aes(year, ymin = lwr, ymax = upr, y = est)) +
  geom_ribbon() +
  geom_line()

true_abund <- tibble(
  year = unique(df_sim$year),
  N = as.numeric(colSums(survey$I))
) %>%
  mutate(type = "True")

strat_abund <- tibble::as_tibble(survey$total_strat) %>%
  mutate(N = total, type = "Design-based") %>%
  select(year, N, type)

index <-
  mutate(ind_nbinom2, type = "Model-based NB2", N = est) %>%
  bind_rows(mutate(ind_delta_gamma, type = "Model-based delta-Gamma", N = est)) |>
  bind_rows(true_abund) %>%
  bind_rows(strat_abund)

index |>
  group_by(type) |>
  mutate(lwr = lwr / median(N)) |> # quick hack, haven't got area-swept right
  mutate(upr = upr / median(N)) |>
  mutate(N = N / median(N)) |>
  ggplot(aes(year, N, group = type, fill = type)) +
  geom_line(aes(colour = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2)
