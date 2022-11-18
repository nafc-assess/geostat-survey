#' Function of cod-like species population simulation based on SimSurvey package
#' The parameters of population dynamics and spatial patterns are given in Supplementary Table 1
#'
#' @param iter: # of simulation
#'
population_cod <- function(iter) {
  set.seed(iter * 35)
  pop <- sim_abundance(ages = 1:20,
                       years = 1:20,
                       Z = sim_Z(log_mean = log(0.4),
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
                                                         sigma = 0.5, log_space = TRUE, plot=FALSE)) |>
    SimSurvey::sim_survey(n_sims = 1,
                          q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}
