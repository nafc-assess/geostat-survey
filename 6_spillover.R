### Spillover effect

#####

##  Survey scenarios

library(data.table)
library(sp)
library(sf)
library(SimSurvey)
library(tidyverse)
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
grid_xy$mpa[grid_xy$x < 50 & grid_xy$x > -120 & grid_xy$y > -65 & grid_xy$y < 100] <- 0.5
grid_xy$mpa[grid_xy$x < 60 & grid_xy$x > -110 & grid_xy$y > -55 & grid_xy$y < 90] <- 1
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
                     coeff = list(alpha = 0, mu = 200, sigma = 70, beta = 5))

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


