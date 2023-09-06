# ------------------------------------------------------------------------
# This R code demonstrates the methods for:
# Plotting the main figures in:
# Yalcin et al. (2023). "Exploring the limits of spatiotemporal and design-based index standardization under reduced survey coverage."
# ICES JMS. doi:
# ------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(here)
library(sf)
library(sp)
library(ggpubr)
library(tidyverse)

# Load data
load(here("data", "index_cod_all_scenarios_200L.Rdata"))
load(here("data", "index_yellowtail_all_scenarios_200L.Rdata"))

# Combining cod-like and yellowtail-like results
index_all_scenarios <- rbind(index_cod_all_scenarios, index_yellowtail_all_scenarios)

# Arranging the levels of the dat for plotting
index_all_scenarios$type <- factor(index_all_scenarios$type, levels = c("TW + Depth",
                                                                        "TW",
                                                                        "NB + Depth",
                                                                        "NB",
                                                                        "DG + Depth",
                                                                        "DG",
                                                                        "Design-based",
                                                                        "Bootstrapped"))


index_all_scenarios$scenario <- factor(index_all_scenarios$scenario, levels = c("Base",
                                                                                "Set reduction" ,
                                                                                "Strata removal",
                                                                                "Area blocked",
                                                                                "Recovery",
                                                                                "Recovery + Spillover"))

load(here("data", "cell.Rdata"))
load(here("data", "strat.Rdata"))

load(here("data", "setdet_cod_base.Rdata"))
load(here("data", "setdet_cod_SR.Rdata"))
load(here("data", "setdet_cod_r30.Rdata"))
load(here("data", "setdet_cod_b30.Rdata"))
load(here("data", "setdet_cod_rec.Rdata"))
load(here("data", "setdet_cod_so.Rdata"))

######################### Creating a polygon for the MPA

mpa <- st_polygon(list(cbind(c(-100, 50, 50, -100, -100), c(-70, -70, 110, 110, -70))))
buffer <- st_polygon(list(cbind(c(-110, 60, 60, -110, -110), c(-80, -80, 120, 120, -80))))
spill <- st_difference(buffer, mpa)

blocked_strat <- strat |> filter(strat %in% c(2,3,4,5,17,18,19,20))

# ------------------------------------------------------------------------

# Figure 2: Scenario representations

# ------------------------------------------------------------------------

setdet_cod_sf <- NULL # Base scenario
for(i in seq_along(setdet_cod)){
  coordinates(setdet_cod[[i]]) = cbind(setdet_cod[[i]]$x, setdet_cod[[i]]$y)
  setdet_cod_sf[[i]] <- st_as_sf(setdet_cod[[i]])
}

setdet_cod_sf_r30 <- NULL # 30% set reduction
for(i in seq_along(setdet_cod_r30)){
  coordinates(setdet_cod_r30[[i]]) = cbind(setdet_cod_r30[[i]]$x, setdet_cod_r30[[i]]$y)
  setdet_cod_sf_r30[[i]] <- st_as_sf(setdet_cod_r30[[i]])
}

setdet_cod_sf_b30 <- NULL # 30% area blocked
for(i in seq_along(setdet_cod_b30)){
  coordinates(setdet_cod_b30[[i]]) = cbind(setdet_cod_b30[[i]]$x, setdet_cod_b30[[i]]$y)
  setdet_cod_sf_b30[[i]] <- st_as_sf(setdet_cod_b30[[i]])
}

setdet_cod_sf_SR <- NULL # Strata removal scenario
for(i in seq_along(setdet_cod_SR)){
  coordinates(setdet_cod_SR[[i]]) = cbind(setdet_cod_SR[[i]]$x, setdet_cod_SR[[i]]$y)
  setdet_cod_sf_SR[[i]] <- st_as_sf(setdet_cod_SR[[i]])
}

setdet_cod_sf_rec <- NULL # Recovery
for(i in seq_along(setdet_cod_rec)){
  coordinates(setdet_cod_rec[[i]]) = cbind(setdet_cod_rec[[i]]$x, setdet_cod_rec[[i]]$y)
  setdet_cod_sf_rec[[i]] <- st_as_sf(setdet_cod_rec[[i]])
}

setdet_cod_sf_so <- NULL # Recovery + Spillover
for(i in seq_along(setdet_cod_so)){
  coordinates(setdet_cod_so[[i]]) = cbind(setdet_cod_so[[i]]$x, setdet_cod_so[[i]]$y)
  setdet_cod_sf_so[[i]] <- st_as_sf(setdet_cod_so[[i]])
}


base <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "a) Base scenario") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


r30 <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_r30[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_r30[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "b) Set reduction")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sr <- ggplot() +
  geom_sf(data = blocked_strat, fill = alpha("yellow",0.2)) +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_SR[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_SR[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "c) Strata removal")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


b30 <- ggplot() +
  geom_sf(data = mpa, fill = alpha("yellow",0.2)) +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_b30[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_b30[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange")) +
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right") +
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "d) Area blocked")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rec <- ggplot() +
  geom_sf(data = mpa, fill = alpha("yellow",0.2)) +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_rec[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_rec[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange")) +
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right") +
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "e) Recovery")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


so <- ggplot() +
  geom_sf(data = mpa, fill = alpha("yellow",0.2)) +
  geom_sf(data = spill, fill = alpha("yellow",0.5)) +
  geom_sf(data = strat, mapping = aes(), fill = "grey99", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_so[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_so[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange")) +
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right") +
  theme(text = element_text(size = 12),
        strip.background = element_rect(fill = "grey99"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 12)) +
  labs(title = "f) Recovery + Spillover") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


figure2 <- ggarrange(base, r30, sr, b30, rec, so, ncol = 6, nrow = 1, legend = "bottom", common.legend = TRUE)
figure2

ggsave("data/figure2_samples_66.pdf", plot = figure2, width = 15, height = 3, units = "in", dpi = 500, bg = "white")

# ------------------------------------------------------------------------

# Figure 3: Example time series of estimated abundance indices

# ------------------------------------------------------------------------

colour_pal <-  c("Bootstrapped" = "#FDBF6F",
                 "DG" = "#1F78B4",
                 "DG + Depth" = "#A6CEE3",
                 "NB" = "#33A02C",
                 "NB + Depth" =  "#B2DF8A",
                 "TW" = "#E31A1C",
                 "TW + Depth" = "#FB9A99",
                 "Bootstrapped" = "#FDBF6F")

cod_ts <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 20)|>
              filter(species == "Cod-like") |>
              filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), linewidth = 1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 1 | pop == 20)|>
                filter(species == "Cod-like") |>
                filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
              aes(x =  year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1) +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 20)|>
              filter(species == "Cod-like"), aes(year, true), linewidth = 1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator") +
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("1" = "Population 1",
                                                                       "20" = "Population 2"))) +
  scale_fill_manual(values = colour_pal[1:3]) +
  scale_colour_manual(values = colour_pal[1:3]) +
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Cod-like")

yt_ts <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 110 | pop == 2)|>
              filter(species == "Yellowtail-like") |>
              filter(type == "Bootstrapped" | type == "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), linewidth = 1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 110 | pop == 2)|>
                filter(species == "Yellowtail-like") |>
                filter(type == "Bootstrapped" | type == "DG" | type == "DG + Depth"),
              aes(x =  year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1)+
  geom_line(data = index_all_scenarios |>
              filter(pop == 110 | pop == 2)|>
              filter(species == "Yellowtail-like"), aes(year, true), linewidth = 1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator") +
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("110" = "Population 1",
                                                                       "2" = "Population 2"))) +
  scale_fill_manual(values = colour_pal[1:3]) +
  scale_colour_manual(values = colour_pal[1:3]) +
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Yellowtail-like")

figure3 <- ggarrange(cod_ts, yt_ts, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure3

ggsave("data/figure3_timeseries.pdf", plot = figure3, width = 12, height = 8, units = "in", dpi = 500, bg = "white")

# ------------------------------------------------------------------------

# Figure 4: Distributions of root mean squared log error (RMSLE) and mean relative error (MRE)

# ------------------------------------------------------------------------

index_all_scenarios_accuracy <-
  index_all_scenarios |>
  filter(year > 10) |>
  ungroup() |>
  filter(type !=  "Bootstrapped") |>
  filter(type !=  "Design-based" | scenario !=  "Strata removal") |>
  group_by(pop, type, scenario, species) |>
  summarise(MRE = mean((N - true) / true),
            RMSE =  sqrt(mean((log(N) - log(true))^2)))


index_all_scenarios_accuracy_long <- index_all_scenarios_accuracy |>
  pivot_longer(-c(type, pop, scenario, species))

metrics_mean <- index_all_scenarios_accuracy_long |>
  group_by(name, type, species, scenario) |>
  summarise(m = mean(value))

rmse_plot <- metrics_mean |>
  filter(name == "RMSE") |>
  ggplot(aes(x =  m, y =  type)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_violin(data = index_all_scenarios_accuracy_long |> filter(name == "RMSE"),
              aes(x =  value, y =  type, col = type), alpha = 0.5) +
  geom_point(aes(), colour = "black", size = 1) +
  #facet_grid(species ~ scenario) +
  facet_grid(species ~ scenario, scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position = NULL) +
  labs(x = "RMSE", y = "", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14)) +
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Root mean square error") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #xlim(NA, 2.5)

mre_plot <- metrics_mean |>
  filter(name == "MRE") |>
  ggplot(aes(x =  m, y =  type)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_violin(data = index_all_scenarios_accuracy_long |> filter(name == "MRE"),
              aes(x =  value, y =  type, col = type), alpha = 0.5) +
  geom_point(aes(), colour = "black", size = 1) +
  #facet_grid(species ~ scenario) +
  facet_grid(species ~ scenario, scales = "free") +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position =  NULL) +
  labs(x = "MRE", y = "", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14)) +
  theme(strip.background = element_rect(fill = "grey97")) +
  #xlim(-1, 2) +
  labs(title = "b) Mean relative error") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



figure4_performance <- ggarrange(rmse_plot, mre_plot, ncol = 1, nrow = 2, legend = "none")
figure4_performance

ggsave("data/figure4_performance.pdf", plot = figure4_performance, width = 12, height = 10, units = "in", dpi = 500, bg = "white")

# ------------------------------------------------------------------------

# Figure 5: Coverage of the 95% confidence intervals

# ------------------------------------------------------------------------

CI_coverage <- index_all_scenarios |>
  filter(type !=  "Design-based") |>
  filter(type !=  "Bootstrapped" | scenario !=  "Strata removal") |>
  mutate(covered = lwr < true & upr > true) |>
  group_by(type, scenario, species) |>
  summarise(mc = mean(covered)) |>
  group_by(scenario) |>
  arrange(mc)

CI_width <- index_all_scenarios |>
  filter(type !=  "Design-based") |>
  filter(type !=  "Bootstrapped" | scenario !=  "Strata removal") |>
  mutate(CI_width = upr-lwr) |>
  group_by(type, scenario, species) |>
  summarise(mCI_width = mean(CI_width))

CI_width_per_pop <- index_all_scenarios |>
  filter(type !=  "Design-based") |>
  filter(type !=  "Bootstrapped" | scenario !=  "Strata removal") |>
  mutate(CI_width = upr-lwr) |>
  group_by(pop, type, scenario, species) |>
  summarise(mCI_width = mean(CI_width))

CI_coverage_plot <- CI_coverage |>
  ggplot(aes(x = mc, y = type, fill = type, group = paste(scenario, type))) +
  geom_point(position = position_dodge(width = 0.6), mapping = aes(colour = type), size = 3) +
  geom_linerange(xmin = 0, mapping = aes(xmax = mc, colour = type), position = position_dodge(width = 0.6)) +
  facet_grid(species ~ scenario)+
  scale_x_continuous(breaks = seq(0, 2, 0.2))+
  scale_colour_manual(values = c("Bootstrapped" = "#FDBF6F",
                                 "DG" = "#1F78B4",
                                 "DG + Depth" = "#A6CEE3",
                                 "NB" = "#33A02C",
                                 "NB + Depth" =  "#B2DF8A",
                                 "TW" = "#E31A1C",
                                 "TW + Depth" = "#FB9A99"))+
  labs(y = "", x = "CI coverage", colour = "Estimator", fill = "Estimator")+
  labs(title = " a) Confidence interval coverage") +
  coord_cartesian(xlim = c(0.45, 1), expand = FALSE) +
  theme_bw() +
  geom_vline(xintercept = 0.95, linetype = 2)+
  theme(legend.position = NULL) +
  theme(text = element_text(size = 14)) +
  theme(strip.background = element_rect(fill = "grey97"))


CI_width_plot <- ggplot(CI_width) +
  geom_violin(data= CI_width_per_pop |> filter(type !=  "Bootstrapped" | scenario !=  "Strata removal"), aes(log(mCI_width), type, col = type),
              alpha = 0.5)+
  geom_point(aes(log(mCI_width), type), colour = "black", size=1)+
  facet_grid(species ~ scenario)+
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Bootstrapped" = "#FDBF6F"))+
  theme(legend.position =  NULL)+
  labs(x = "CI width", y = "", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Confidence interval width")


figure5_CI <- ggarrange(CI_coverage_plot, CI_width_plot, ncol = 1, nrow = 2, legend = "none")
figure5_CI


ggsave("data/figure5_CI.pdf", plot = figure5_CI, width = 12, height = 10, units = "in", dpi = 500, bg = "white")
