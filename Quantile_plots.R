library(here)
library(ggplot2)
library(tidyr)
library(dplyr)

load(here("Data", "index_all_scenarios.Rdata"))

### Calculating MRE, RMSLE, and pulling AIC
index_all_scenarios_accuracy <-
  index_all_scenarios |>
  filter(year > 10) |>
  ungroup() |>
  filter(type != "Bootstrapped") |>
  group_by(pop, type, scenario, species) |>
  summarise(MRE = mean((log(N) - log(true)) / log(true)),
            RMSLE= sqrt(mean((log(N) - log(true))^2)),
            AIC = mean(AIC))

### the mean of values that are worse than the 90th percentile
metric_q90_worst <- index_all_scenarios_accuracy |>
  group_by(species, scenario) |>
  ### calculating the 90th percentile by species and scenario (averaging populations and estimation types)
  mutate(AIC_q90 = quantile(AIC, c(0.90), na.rm = TRUE),
         MRE_q90 = quantile(abs(MRE), c(0.90), na.rm = TRUE),
         RMSLE_q90 = quantile(abs(RMSLE), c(0.90), na.rm = TRUE)) |>
  ungroup() |>
  ### Only pooling populations, if the mean value is worse than the 90th percentile
  group_by(species, scenario, type) |>
  summarise(AIC_worse = mean(AIC > AIC_q90),
            MRE_worse = mean(abs(MRE) > MRE_q90),
            RMSLE_worse = mean(RMSLE > RMSLE_q90))

### long-data format to plot
metric_q90_worst_long <- metric_q90_worst |>
  pivot_longer(-c(scenario, species, type))

### Cod-like plot
metric_q90_worst_long |>
  filter(species == "Cod-like") |>
  ggplot(aes(x =  scenario, y =  value, group = type, color = type)) +
  geom_line(aes(), linewidth = 1) +
  geom_point(aes(), size = 2) +
  facet_grid(type ~ name) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position = NULL)+
  labs(x = "Scenario", y = "The ratio of the mean of values that are worse than the 90th percentile", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Cod-like species")

### Yellowtail-like plot
metric_q90_worst_long |>
  filter(species == "Yellowtail-like") |>
  ggplot(aes(x =  scenario, y =  value, group = type, color = type)) +
  geom_line(aes(), linewidth = 1) +
  geom_point(aes(), size = 2) +
  facet_grid(type ~ name) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position = NULL)+
  labs(x = "Scenario", y = "The ratio of the mean of values that are worse than the 90th percentile", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Yellowtail-like species")
