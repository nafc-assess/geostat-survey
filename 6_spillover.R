### Spillover effect

#####

#############  Packages
library(sp)
library(sf)
library(SimSurvey)
library(tidyverse)
library(sdmTMB)
library(tidyr)
library(future)
library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)
library(ggpubr)
library(mapview)


## Sim non linear function -----------------------------------------------------
## Idea is to increase the probability the fish will occur in the MPA

grid <- make_grid(x_range = c(-150, 150),
                  y_range = c(-150, 150),
                  res = c(10, 10),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1000),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 40),
                  strat_splits = 2,
                  method = "spline")

## Hack: add mpa column to grid object
grid_xy <- as.data.frame(grid)
grid_xy$mpa <- 0
grid_xy$mpa[grid_xy$x < 50 & grid_xy$x > -120 & grid_xy$y > -65 & grid_xy$y < 100] <- 0.8
grid_xy$mpa[grid_xy$x < 30 & grid_xy$x > -100 & grid_xy$y > -45 & grid_xy$y < 80] <- 1
grid$mpa <- grid_xy$mpa

plot(x ~ y, data = grid_xy, col = grid_xy$mpa)

sim <- list(ages = 1:20, years = 1:20)
grid_dat <- data.table::data.table(as.data.frame(grid))
grid_edat <- grid_dat
i <- rep(seq(nrow(grid_edat)), times = length(sim$ages))
a <- rep(sim$ages, each = nrow(grid_edat))
grid_edat <- grid_edat[i]
grid_edat$age <- a
i <- rep(seq(nrow(grid_edat)), times = length(sim$years))
y <- rep(sim$years, each = nrow(grid_edat))
grid_edat <- grid_edat[i]
grid_edat$year <- y
grid_edat <- grid_edat[order(grid_edat$cell, grid_edat$year, grid_edat$age), ] # sort to align with error array


depth_par <- sim_nlf(formula = ~ alpha + ifelse(year > 10, (beta * mpa * (year - 10) / 10), 0) - ((depth - mu) ^ 2) / (2 * sigma ^ 2),
                     coeff = list(alpha = 0, mu = 160, sigma = 60, beta = 5))


depth <- depth_par(grid_edat)

grid_edat$depth_effect <- depth

library(plotly)
grid_edat %>%
  plot_ly(x = ~depth, y = ~exp(depth_effect), split = ~mpa,
          frame = ~year) %>%
  add_lines()


sim <- sim_abundance(ages = 1:20, years = 1:20) %>%
  sim_distribution(grid = grid,
                   ays_covar = sim_ays_covar(phi_age = 0.8,
                                             phi_year = 0.1),
                   depth_par = depth_par) # |> sim_survey(n_sims = 1, q = sim_logistic(k = 2, x0 = 3))

plot_distribution(sim, ages = 5, years = 8:20, scale = "log", type = "heatmap")
# plot_survey(sim, which_year = 20, which_sim = 1)

# ################ 2 functions for spillover and no spillover effect

population_cod_spillover <- function(iter) {
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

    SimSurvey::sim_distribution(grid = grid,
                                ays_covar = sim_ays_covar(sd = 2.8,
                                                          range = 300,
                                                          phi_age = 0.4,
                                                          phi_year = 0.8,
                                                          group_ages = 5:20),
                                depth_par = depth_par)|>
    SimSurvey::sim_survey(n_sims = 1,
                          q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}


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
                                depth_par = sim_parabola(mu = 160,
                                                         sigma = 60, log_space = FALSE, plot=FALSE)) |>
    SimSurvey::sim_survey(n_sims = 1,
                          q = sim_logistic(k = 2, x0 = 3, plot = FALSE),
                          trawl_dim = c(3, 0.02),
                          resample_cells = FALSE,
                          binom_error = TRUE,
                          min_sets = 2,
                          set_den = 2/1000)
}

set.seed(1)
survey_cod_spillover <- furrr::future_map(seq_len(1), population_cod_spillover, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()

set.seed(1)
survey_cod <- furrr::future_map(seq_len(1), population_cod, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
gc()


#### Calculating the all age numbers per grid and year

# Spillover
pop_dist_age_cod_spill <- map(survey_cod_spillover, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod_spill <- NULL
for(i in seq_along(pop_dist_age_cod_spill)){
  pop_dist_cod_spill[[i]] <- pop_dist_age_cod_spill[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))|>
    mutate(status = "Spillover effect")}

# No spillover
pop_dist_age_cod <- map(survey_cod, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod <- NULL
for(i in seq_along(pop_dist_age_cod)){
  pop_dist_cod[[i]] <- pop_dist_age_cod[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I)) |>
    mutate(status = "No spillover effect")}


## plotting the map
combine_status <- rbind(pop_dist_cod[[1]], pop_dist_cod_spill[[1]])

combine_status|>
  filter(year > 10) |>
  ggplot(aes(x, y, fill = I_all_age)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_grid(status ~ year, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15",
                                                    "16" = "Year 16",
                                                    "17" = "Year 17",
                                                    "18" = "Year 18",
                                                    "19" = "Year 19",
                                                    "20" = "Year 20"))) +
  coord_fixed(expand = FALSE) +
  theme_bw()+
  theme(text = element_text(size = 14))+
  theme(legend.position = "bottom", legend.spacing.x = unit(0.5, 'cm'), legend.direction="horizontal")+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10))+
  theme(strip.background =element_rect(fill="grey97"))+
  theme(panel.spacing = unit(1, "lines"))+
  labs(title = "Available population distribution") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())


#######################

# in the mpa and buffer area

combine_status_list <- list(pop_dist_cod[[1]], pop_dist_cod_spill[[1]])

combine_status_list

for( i in seq_along(combine_status_list)){
  combine_status_list[[i]]$mpa <- 0
  combine_status_list[[i]]$mpa[combine_status_list[[i]]$x < 50 & combine_status_list[[i]]$x > -120 & combine_status_list[[i]]$y > -65 & combine_status_list[[i]]$y < 100] <- 0.8
  combine_status_list[[i]]$mpa[combine_status_list[[i]]$x < 30 & combine_status_list[[i]]$x > -100 & combine_status_list[[i]]$y > -45 & combine_status_list[[i]]$y < 80] <- 1
}

combine_status_list_db <- rbind(combine_status_list[[1]], combine_status_list[[2]])

a <- combine_status_list_db |>
  group_by(status, year) |>
  summarise(N_total = sum(N_all_age),
            I_total = sum(I_all_age))|>ggplot() +
  geom_line(aes(year, I_total, colour = status, group = status), linewidth = 1) +
  facet_grid(~ status)+
  labs(title = "Total available population") +
  theme_bw()

b <- combine_status_list_db |>
  group_by(mpa, status, year) |>
  summarise(N_total = sum(N_all_age),
            I_total = sum(I_all_age))|>
ggplot() +
  geom_line(aes(year, I_total, colour = status, group = status), linewidth = 1) +
  facet_grid(~ mpa, labeller = labeller(mpa =
                                         c("0" = "Outside MPA and buffer area",
                                           "0.8" = "Buffer area",
                                           "1" = "MPA"))) +
  labs(title = "Available population by MPA status") +
  theme_bw()


spillover_numbers <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
ggsave("spillover_numbers.png", plot = spillover_numbers, width = 10, height = 12, units = "in", dpi = 500, bg="white")

spillover_numbers
