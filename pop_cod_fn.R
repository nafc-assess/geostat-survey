#' Function of cod-like species population simulation based on SimSurvey package
#' The parameters of population dynamics and spatial patterns are given in Supplementary Table 1
#'
#' @param iter: # of simulation
#'
population_cod <- function(iter, n_sims) {
  set.seed(iter * 35)
  pop <- sim_abundance(ages = 1:20,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.3), # changed this from 0.4 to 0.3
                                 log_sd = 0.05,
                                 phi_age = 0.9,
                                 phi_year = 0.5,
                                 plot = FALSE),

                       R = sim_R(log_mean = log(60e+06),
                                 log_sd = 0.4,
                                 random_walk = FALSE,
                                 plot = FALSE),

                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),

                       growth = sim_vonB(Linf = 120,
                                         L0 = 5,
                                         K = 0.11,
                                         log_sd = 0.15,
                                         plot = FALSE,
                                         length_group = 3))|>

    SimSurvey::sim_distribution(grid = make_grid(x_range = c(-150, 150),
                                                 y_range = c(-150, 150),
                                                 res = c(10, 10),
                                                 shelf_depth = 200,
                                                 shelf_width = 100,
                                                 depth_range = c(0, 1000),
                                                 n_div = 1,
                                                 strat_breaks = seq(0, 1000, by = 40),
                                                 strat_splits = 2,
                                                 method = "spline"),
                                ays_covar = sim_ays_covar(sd = 2.8,
                                                          range = 300,
                                                          phi_age = 0.4,
                                                          phi_year = 0.8,
                                                          group_ages = 5:20),

                                depth_par = sim_parabola(mu = log(160),
                                                         sigma = 0.5,
                                                         log_space = TRUE, plot=FALSE)) |>
    SimSurvey::sim_survey(n_sims = n_sims,
                          q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = set_den = 2/1000)
}


#### Recovery Scenario

grid_with_mpa_r <- make_grid(x_range = c(-150, 150),
                             y_range = c(-150, 150),
                             res = c(10, 10),
                             shelf_depth = 200,
                             shelf_width = 100,
                             depth_range = c(0, 1000),
                             n_div = 1,
                             strat_breaks = seq(0, 1000, by = 40),
                             strat_splits = 2,
                             method = "spline")

grid_xy <- as.data.frame(grid_with_mpa_r)
grid_xy$mpa <- 0
grid_xy$mpa[grid_xy$x < 50 & grid_xy$x > -100 & grid_xy$y > -70 & grid_xy$y < 110] <- 1
grid_with_mpa_r$mpa <- grid_xy$mpa

depth_mpa_par_recovery <- sim_nlf(formula = ~ alpha + ifelse(year > 10, (beta * mpa * (year - 10) / 10), 0) - ((log(depth) - mu) ^ 2) / (2 * sigma ^ 2),
                                  coeff = list(alpha = 0,
                                               mu = log(160),
                                               sigma = 0.5,
                                               beta = 5))

population_cod_recovery <- function(iter, n_sims) {
  set.seed(iter * 35)
  pop <- sim_abundance(ages = 1:20,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.3), # changed this from 0.4 to 0.3
                                 log_sd = 0.05,
                                 phi_age = 0.9,
                                 phi_year = 0.5,
                                 plot = FALSE),

                       R = sim_R(log_mean = log(60e+06),
                                 log_sd = 0.4,
                                 random_walk = FALSE,
                                 plot = FALSE),

                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),

                       growth = sim_vonB(Linf = 120,
                                         L0 = 5,
                                         K = 0.11,
                                         log_sd = 0.15,
                                         plot = FALSE,
                                         length_group = 3))|>

    SimSurvey::sim_distribution(grid = grid_with_mpa_r,
                                ays_covar = sim_ays_covar(sd = 2.8,
                                                          range = 300,
                                                          phi_age = 0.4,
                                                          phi_year = 0.8,
                                                          group_ages = 5:20),
                                depth_par = depth_mpa_par_recovery) |>

    SimSurvey::sim_survey(n_sims = n_sims,
                          q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}


#### Spillover effect
# Create grid with mpa column
grid_with_mpa <- make_grid(x_range = c(-150, 150),
                           y_range = c(-150, 150),
                           res = c(10, 10),
                           shelf_depth = 200,
                           shelf_width = 100,
                           depth_range = c(0, 1000),
                           n_div = 1,
                           strat_breaks = seq(0, 1000, by = 40),
                           strat_splits = 2,
                           method = "spline")

grid_xy <- as.data.frame(grid_with_mpa)
grid_xy$mpa <- 0
grid_xy$mpa[grid_xy$x < 60 & grid_xy$x > -110 & grid_xy$y > -80 & grid_xy$y < 120] <- 0.7
grid_xy$mpa[grid_xy$x < 50 & grid_xy$x > -100 & grid_xy$y > -70 & grid_xy$y < 110] <- 1
grid_with_mpa$mpa <- grid_xy$mpa

depth_mpa_par <- sim_nlf(formula = ~ alpha + ifelse(year > 10, (beta * mpa * (year - 10) / 10), 0) - ((log(depth) - mu) ^ 2) / (2 * sigma ^ 2),
        coeff = list(alpha = 0,
                     mu = log(160),
                     sigma = 0.5,
                     beta = 5))

population_cod_spillover <- function(iter, n_sims) {
  set.seed(iter * 35)
  pop <- sim_abundance(ages = 1:20,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.3), # changed this from 0.4 to 0.3
                                 log_sd = 0.05,
                                 phi_age = 0.9,
                                 phi_year = 0.5,
                                 plot = FALSE),

                       R = sim_R(log_mean = log(60e+06),
                                 log_sd = 0.4,
                                 random_walk = FALSE,
                                 plot = FALSE),

                       N0 = sim_N0(N0 = "exp",
                                   plot = FALSE),

                       growth = sim_vonB(Linf = 120,
                                         L0 = 5,
                                         K = 0.11,
                                         log_sd = 0.15,
                                         plot = FALSE,
                                         length_group = 3))|>

    SimSurvey::sim_distribution(grid = grid_with_mpa,
                                ays_covar = sim_ays_covar(sd = 2.8,
                                                          range = 300,
                                                          phi_age = 0.4,
                                                          phi_year = 0.8,
                                                          group_ages = 5:20),
                                depth_par = depth_mpa_par) |>

                                  SimSurvey::sim_survey(n_sims = n_sims,
                                                        q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                                                        trawl_dim = c(3, 0.02),
                                                        resample_cells = FALSE,
                                                        binom_error = TRUE,
                                                        min_sets = 2,
                                                        set_den = 2/1000)
}
