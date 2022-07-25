
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
n_boot <- 5000


## Simulation --------------------------------------------------------------------------------------

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


## Density from the Gamma distribution -------------------------------------------------------------

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


### Density from bootstrapping ---------------------------------------------------------------------

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
    mutate(samp = seq.int(reps), sim = mean(data$sim), year = mean(data$year))
  return(boot)
}

tic()
boot_index <- furrr::future_map_dfr(split_setdet, boot_one_year, reps = n_boot,
                                    .options = furrr::furrr_options(seed = TRUE))
toc()

quantile(boot_index$total, prob = c(0.001, 0.999))

den_plot <- ggplot() +
  geom_density_ridges(aes(x = total, y = as.numeric(year), group = factor(year)),
                      color = "grey90", fill = "steelblue", alpha = 0.7,
                      data = boot_index, scale = 1) +
  geom_density_ridges(aes(x = total, y = year, height = den, group = factor(year)),
                      stat = "identity", color = "grey90", fill = "red", alpha = 0.7,
                      data = total_strat_den, scale = -1) +
  coord_flip() + guides(fill = "none") +
  scale_x_continuous(labels = scales::label_number(suffix = "", scale = 1e-8),
                     limits = c(194587641, 5116017391)) +
  ylab("Year") + xlab("Abundance index") +
  facet_grid(rows = "sim") +
  theme_nafo()

saveRDS(total_strat, file = "Gamma_SCR/data/total_strat.rds")
saveRDS(total_strat_den, file = "Gamma_SCR/data/total_strat_den.rds")
saveRDS(boot_index, file = "Gamma_SCR/data/boot_index.rds")
saveRDS(den_plot, file = "Gamma_SCR/data/den_plot.rds")


## Relative status ---------------------------------------------------------------------------------

total_strat <- readRDS("Gamma_SCR/data/total_strat.rds")
total_strat_den <- readRDS("Gamma_SCR/data/total_strat_den.rds")
boot_index <- readRDS("Gamma_SCR/data/boot_index.rds")

sub_total_strat <- total_strat |>
  filter(sim == 1)

ref_est <- total_strat |>
  filter(sim == 1, year %in% 2:9) |>
  summarise(total = mean(total),
            sigma = sqrt(sum(sigma ^ 2) / (2 * n())),
            scale = sigma ^ 2 / total,
            shape = total / scale)

ref_boot <- boot_index |>
  filter(sim == 1, year %in% 2:9) |>
  group_by(samp) |>
  summarise(total = mean(total)) |>
  ungroup()

x <- seq(min(boot_index$total), max(boot_index$total), length.out = 100)
ref_den <- data.frame(total = x, den = dgamma(x, shape = ref_est$shape, scale = ref_est$scale))

t_est <- total_strat |>
  filter(sim == 1, year == 20)

t_den <- total_strat_den |>
  filter(sim == 1, year == 20)

t_boot <- boot_index |>
  filter(sim == 1, year == 20)

boot_prob <- mean((t_boot$total - ref_boot$total) < 0)
n_samp <- 100000
ref_samp <- rgamma(n_samp, shape = ref_est$shape, scale = ref_est$scale)
t_samp <- rgamma(n_samp, shape = t_est$shape, scale = t_est$scale)
gamma_prob <- mean((t_samp - ref_samp) < 0)

ref_plot <- ggplot() +
  geom_density(aes(x = total), data = ref_boot, fill = "steelblue", color = "steelblue", alpha = 0.5) +
  geom_area(aes(x = total, y = -den), data = ref_den, fill = "red", color = "red", alpha = 0.5) +
  geom_density(aes(x = total), data = t_boot, fill = NA, color = "steelblue", size = .nafo_lwd) +
  geom_area(aes(x = total, y = -den), data = t_den, fill = NA, color = "red", size = .nafo_lwd) +
  geom_text(aes(x = t_est$total, y = max(ref_den$den) * 1.2, label = "Terminal estimate"), hjust = 0, vjust = 0.5) +
  geom_text(aes(x = ref_est$total, y = max(ref_den$den) * 1.2, label = "Reference point"), hjust = 0, vjust = 1) +
  geom_text(aes(x = t_est$total, y = 0, label = round(boot_prob, 3)), hjust = -0.5, color = "steelblue") +
  geom_text(aes(x = t_est$total, y = 0, label = round(gamma_prob, 3)), hjust = 1.5, color = "red") +
  theme_nafo() +
  coord_flip() +
  scale_x_continuous(labels = scales::label_number(suffix = "", scale = 1e-8)) +
  ylab("") + xlab("Abundance index") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

save(ref_plot, t_est, ref_est, ref_den, boot_prob, gamma_prob, file = "Gamma_SCR/data/ref_plot.rda")






