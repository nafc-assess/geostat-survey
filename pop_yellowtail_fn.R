#' Function of yellowtail-like species population simulation based on SimSurvey package
#' The parameters of population dynamics and spatial patterns are given in Supplementary Table 1
#'
#' @param iter: # of simulation
#'
population_yellowtail <- function(iter, n_sims) {
  set.seed(iter * 1)
  pop <- SimSurvey::sim_abundance(ages = 1:10,
                                  years = 1:20,

                                  Z = sim_Z(log_mean = log(0.5),
                                            log_sd = 0.2,
                                            phi_age = 0.3,
                                            phi_year = 0.4),

                                  R = sim_R(log_mean = log(8e+09),
                                            log_sd = 0.2,
                                            random_walk = FALSE),

                                  N0 = sim_N0(N0 = "exp", plot = FALSE),
                                  growth = sim_vonB(Linf = 56, L0 = 0,
                                                    K = 0.13, log_sd = 0.12,
                                                    length_group = 2, digits = 0)) |>
    sim_distribution(make_grid(x_range = c(-150, 150),
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
                                               range = 500,
                                               phi_age = 0.4,
                                               phi_year = 0.8,
                                               group_ages = 5:9),

                     depth_par = sim_parabola(mu = log(90),
                                              sigma = 0.3,
                                              log_space = TRUE))
  set.seed(35)
  survey <- SimSurvey::sim_survey(pop, n_sims = n_sims,
                          q = sim_logistic(k = 1.6, x0 = 5.5, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
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
                                               mu = log(90),
                                               sigma = 0.3,
                                               beta = 5))

population_yellowtail_recovery <- function(iter, n_sims) {
  set.seed(iter * 1)
  pop <- SimSurvey::sim_abundance(ages = 1:10,
                                  years = 1:20,

                                  Z = sim_Z(log_mean = log(0.5),
                                            log_sd = 0.2,
                                            phi_age = 0.3,
                                            phi_year = 0.4),

                                  R = sim_R(log_mean = log(8e+09),
                                            log_sd = 0.2,
                                            random_walk = FALSE),

                                  N0 = sim_N0(N0 = "exp", plot = FALSE),
                                  growth = sim_vonB(Linf = 56, L0 = 0,
                                                    K = 0.13, log_sd = 0.12,
                                                    length_group = 2, digits = 0))|>

    SimSurvey::sim_distribution(grid = grid_with_mpa_r,
                                sim_ays_covar(sd = 2.8,
                                              range = 500,
                                              phi_age = 0.4,
                                              phi_year = 0.8,
                                              group_ages = 5:9),
                                depth_par = depth_mpa_par_recovery)
  set.seed(35)
  survey <- SimSurvey::sim_survey(pop, n_sims = n_sims,
                                  q = sim_logistic(k = 1.6, x0 = 5.5, plot = FALSE),
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
grid_xy$mpa[grid_xy$x < 60 & grid_xy$x > -110 & grid_xy$y > -80 & grid_xy$y < 120] <- 0.5
grid_xy$mpa[grid_xy$x < 50 & grid_xy$x > -100 & grid_xy$y > -70 & grid_xy$y < 110] <- 1
grid_with_mpa$mpa <- grid_xy$mpa

depth_mpa_par <- sim_nlf(formula = ~ alpha + ifelse(year > 10, (beta * mpa * (year - 10) / 10), 0) - ((log(depth) - mu) ^ 2) / (2 * sigma ^ 2),
                         coeff = list(alpha = 0,
                                      mu = log(90),
                                      sigma = 0.3,
                                      beta = 5))

population_yellowtail_spillover <- function(iter, n_sims) {
  set.seed(iter * 1)
  pop <- SimSurvey::sim_abundance(ages = 1:10,
                                  years = 1:20,

                                  Z = sim_Z(log_mean = log(0.5),
                                            log_sd = 0.2,
                                            phi_age = 0.3,
                                            phi_year = 0.4),

                                  R = sim_R(log_mean = log(8e+09),
                                            log_sd = 0.2,
                                            random_walk = FALSE),

                                  N0 = sim_N0(N0 = "exp", plot = FALSE),
                                  growth = sim_vonB(Linf = 56, L0 = 0,
                                                    K = 0.13, log_sd = 0.12,
                                                    length_group = 2, digits = 0))|>

    SimSurvey::sim_distribution(grid = grid_with_mpa,
                                sim_ays_covar(sd = 2.8,
                                              range = 500,
                                              phi_age = 0.4,
                                              phi_year = 0.8,
                                              group_ages = 5:9),
                                depth_par = depth_mpa_par)
    set.seed(35)
    survey <- SimSurvey::sim_survey(pop, n_sims = n_sims,
                          q = sim_logistic(k = 1.6, x0 = 5.5, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}


