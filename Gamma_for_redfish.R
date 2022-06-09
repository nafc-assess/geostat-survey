# redfish

#############  Packages
library(SimSurvey)
library(sdmTMB)
library(tidyr)
library(future)
library(tictoc)
library(ggplot2)
library(data.table)
library(dplyr)
library(purrr)

plan(multisession, workers = floor(availableCores()/2))

set.seed(794)
population <- sim_abundance(ages = 1:50,
                       years = 1:20,
                       R = sim_R(log_mean = log(600000000),
                                 log_sd = 0.6,
                                 random_walk = F),
                       Z = sim_Z(log_mean = log(0.2),
                                 log_sd = 0.2,
                                 phi_age = 0.4,
                                 phi_year = 0.4),
                       N0 = sim_N0(N0 = "exp", plot = FALSE),
                       growth = sim_vonB(Linf = 30, L0 = 0,      #roughly based on Cadigan & Compana 2016
                                         K = 0.1, log_sd = 0.13,
                                         length_group = 1, digits = 0)) |>
    sim_distribution(grid = make_grid(x_range = c(-150, 150),
                                      y_range = c(-150, 150),
                                      res = c(10, 10),
                                      shelf_depth = 60,
                                      shelf_width = 170,
                                      depth_range = c(0, 1600),
                                      n_div = 2,
                                      strat_breaks = seq(0, 1600, by = 65),
                                      strat_splits = 4,
                                      method = "bezier"),
                     ays_covar = sim_ays_covar(sd = 2,
                                               range = 200,
                                               #lambda = 0.5,
                                               #model = "matern",
                                               phi_age = 0.5,
                                               phi_year = 0.9,
                                               #group_ages = c(1,20:24)
                     ),
                     depth_par = sim_parabola(mu = log(190),
                                              sigma = 0.3,
                                              log_space = TRUE))


survey_fn <- function(n_sims) {
    sim_survey(population,
               n_sims = n_sims,
               q = sim_logistic(k = 1, x0 = 6.5),
               trawl_dim = c(1.5, 0.02),
               resample_cells = FALSE,
               binom_error = TRUE,
               min_sets = 2,
               set_den = 1/1000,
               lengths_cap = 250,
               ages_cap = 20,
               age_sampling = "stratified",
               age_length_group = 1,
               age_space_group = "division") |>
    run_strat()
}

nsims=replicate(10, 100, FALSE) ## each run is for 100 sims within it - 10 runs in total. Number of simulation can be increased here.

tic()
survey <- furrr::future_map(nsims, survey_fn, .options = furrr::furrr_options(seed = TRUE, packages = "SimSurvey"))
toc()


### combining the design-based estimations from all runs

total_strats <- map(survey, function(x) {pluck(x, 'total_strat')})

for( i in seq_along(total_strats)){
  total_strats[[i]]$sim_number <- total_strats[[i]]$sim + ((as.numeric(i) - 1) * 100)
} # Updating the unique sim number (since each batch has 100 sims within it)

total_strats_dfr <- do.call(rbind.data.frame, total_strats)


### Filtering the year 20

total_strats_dfr_y20 <- total_strats_dfr |> filter(year == 20)

ggplot(total_strats_dfr_y20, aes(total)) +
  geom_histogram(bins = 20)


### Gamma prob. dist.

## Overall

# Gamma estimators
sigma = mean(total_strats_dfr_y20$sd * total_strats_dfr_y20$sampling_units)
scale = mean((sigma ^ 2) / total_strats_dfr_y20$total)
shape = mean(total_strats_dfr_y20$total / scale)

# plot 1
hist(as.numeric(total_strats_dfr_y20$total), probability = TRUE, breaks = 20,
     main = "Dist. of Abundance (1000 sims) at Year 20", xlab = "Abundance")
curve(dgamma(x, shape = shape, scale = scale), add = TRUE, col = "darkorange", lwd = 2)

# plot 2
qqplot(x = qgamma(ppoints(total_strats_dfr_y20$total),
                  shape = shape, scale = scale),
       y = total_strats_dfr_y20$total,
       #xlim = c(0, 1986736283), ylim = c(0, 1986736283),
       main = "QQ-Plot: Abundance (1000 sims), Gamma Distribution",
       xlab = "Theoretical Quantiles, Gamma Distribution",
       ylab = "Sample Quantiles, Abundance")
abline(a = 0, b = 1, col = "dodgerblue", lwd = 2)
grid()

## Individual sims: Strata level

total_strats_dfr_y20$sigma <- total_strats_dfr_y20$sampling_units * total_strats_dfr_y20$sd
total_strats_dfr_y20$scale <- (total_strats_dfr_y20$sigma ^ 2) / total_strats_dfr_y20$total
total_strats_dfr_y20$shape <- total_strats_dfr_y20$total / total_strats_dfr_y20$scale

add_curve <- function(plot, i) {
  return(plot + stat_function(aes(x=total_strats_dfr_y20$total[i], colour = "sim"), fun = dgamma, col=i, args = list(shape=total_strats_dfr_y20$shape[i], scale=total_strats_dfr_y20$scale[i])))
}

p1 <- ggplot(total_strats_dfr_y20, aes(x = total))+
  geom_histogram(aes(y=..density..),fill="grey",alpha=0.8, bins = 20)  +
  theme_minimal()


for (i in 1:1000){ # choose how many sims you want to use in the plot
  p1<- add_curve(p1, i)
}

p1


### Bootstapping

sumYst <- function(data, i = seq_len(nrow(data))) {
  data[i, ] |>
    ### stratum level
    group_by(year, strat, strat_area) |>
    summarise(meanYh = mean(n), tow_area = mean(tow_area), .groups = "drop_last") |>
    mutate(Nh = strat_area/(tow_area)) |>
    group_by(year) |>
    mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)|>
    ### year level
    summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") |>
    pull(sumYst)
}

boot_one_year <- function(data, reps) {
  b <- boot::boot(data, statistic = sumYst, strata = data$strat, R = reps)
  suppressWarnings(bci <- boot::boot.ci(b, type = "bca"))
  mean_boot <- mean(b$t)
  sd_boot <- sd(b$t)
  boot <- data.table(
    sim = mean(data$sim_number),
    mean_boot= mean_boot,
    sd_boot = sd_boot,
    lwr = bci$bca[[4]],
    upr = bci$bca[[5]],
    cv = sd_boot / mean_boot,
    type = "Bootstrapped"
  )
  return(boot)
}

setdets <- map(survey, function(x) {pluck(x, 'setdet')}) ### plucking only setdet element

for( i in seq_along(setdets)){
  setdets[[i]]$sim_number <- setdets[[i]]$sim + ((as.numeric(i) - 1) * 100)
} ### Updating the simulation numbers

for( i in seq_along(setdets)){
  setdets[[i]] <- setdets[[i]][setdets[[i]]$year == 20,]
} ### Subsetting only year 20

setdet_sim <- NULL
for( i in seq_along(setdets)){
  setdet_sim[[i]] <- setdets[[i]] |> group_split(sim_number)
} ### Splitting the sims for future map

setdet_sim <- unlist(setdet_sim,recursive=FALSE)


tic()
boot_index <- furrr::future_map_dfr(setdet_sim, boot_one_year, reps=1000, .options = furrr::furrr_options(seed = TRUE))
toc()

