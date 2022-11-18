#' Function of yellowtail-like species population simulation based on SimSurvey package
#' The parameters of population dynamics and spatial patterns are given in Supplementary Table 1
#'
#' @param iter: # of simulation
#'
population_yellowtail <- function(iter) {
  set.seed(iter * 35)
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
                                              #sigma_right = 0.20,
                                              log_space = TRUE)) |>
    SimSurvey::sim_survey(n_sims = 1,
                          q = sim_logistic(k = 1.6, x0 = 5.5, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}
