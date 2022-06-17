
#############  Packages
library(SimSurvey)
library(tidyr)
library(future)
library(tictoc)
library(ggplot2)
library(ggridges)
library(dplyr)
library(purrr)
library(data.table)
library(NAFOdown)

plan(multisession, workers = floor(availableCores()/2))

n_sims <- 5
n_boot <- 1000

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


survey <- sim_survey(population,
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

## Density from the Gamma distribution

total_strat <- survey$total_strat |>
  mutate(sigma = sampling_units * sd,
         scale = sigma ^ 2 / total,
         shape = total / scale)

## Use gamma to generate density by sim and year
rng <- c(0.001, max(total_strat$total) * 2)
x <- seq(rng[1], rng[2], length.out = 100)
total_strat_den <- lapply(seq.int(nrow(total_strat)), function(i) {
  data.frame(sim = total_strat$sim[i],
             year = total_strat$year[i],
             total = x,
             den = dgamma(x, shape = total_strat$shape[i],
                          scale = total_strat$scale[i]))
}) |> dplyr::bind_rows()


### Density from bootstrapping

setdet <- survey$setdet

split_setdet <- split(setdet, paste0(setdet$year, "-", setdet$sim))

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
  boot <- data.table(b$t) |> dplyr::rename(total = V1) |>
    mutate(sim = mean(data$sim), year = mean(data$year))
  return(boot)
}

tic()
boot_index <- furrr::future_map_dfr(split_setdet, boot_one_year, reps = n_boot,
                                    .options = furrr::furrr_options(seed = TRUE))
toc()


den_plot <- ggplot() +
  geom_density_ridges(aes(x = total / 1e+08, y = as.numeric(year), group = factor(year)),
                      color = "grey90", fill = "steelblue", alpha = 0.7,
                      data = boot_index, scale = 1) +
  geom_density_ridges(aes(x = total / 1e+08, y = year, height = den, group = factor(year)),
                      stat = "identity", color = "grey90", fill = "red", alpha = 0.7,
                      data = total_strat_den, scale = -1) +
  coord_flip() + guides(fill = "none") +
  ylab("Year") + xlab("Abundance index") +
  xlim(0, 40) +
  facet_grid(rows = "sim") +
  theme_nafo()

saveRDS(total_strat, file = "Gamma_SCR/data/total_strat.rds")
saveRDS(total_strat_den, file = "Gamma_SCR/data/total_strat_den.rds")
saveRDS(boot_index, file = "Gamma_SCR/data/boot_index.rds")
saveRDS(den_plot, file = "Gamma_SCR/data/den_plot.rds")


