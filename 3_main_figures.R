library(ggplot2)
library(sp)
library(sf)
library(data.table)
library(raster)
library(dplyr)
library(ggpubr)
library(future)
library(purrr)
library(tidyr)
library(facetscales)
library(here)

######################### Loading datasets

here()

load(here("Data", "index_all_scenarios.Rdata"))
load(here("Data", "cell.Rdata"))
load(here("Data", "strat.Rdata"))

load(here("Data", "setdet_cod_base.Rdata"))
load(here("Data", "setdet_cod_SR.Rdata"))
load(here("Data", "setdet_cod_r30.Rdata"))
load(here("Data", "setdet_cod_b30.Rdata"))
load(here("Data", "block_poly_sf_30.Rdata"))

#### Figure 2: Scenario representations

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


base <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey90", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey90"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(title = "a) Base scenario")


r30 <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey90", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_r30[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_r30[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey90"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(title = "b) 30% set reduction")

b30 <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey90", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_b30[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_b30[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  geom_sf(data = block_poly_sf_30, fill = alpha("yellow",0.1)) +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange")) +
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right") +
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey90"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(title = "c) 30% area blocked")


blocked_strat <- strat |>
  filter(strat %in% c(2,3,4,5,17,18,19,20))

sr <- ggplot() +
  geom_sf(data = strat, mapping = aes(), fill = "grey90", alpha = 0.4, colour = "grey20", size = 1, inherit.aes = FALSE) +
  geom_sf(data = setdet_cod_sf_SR[[1]] |> filter(year == 11), mapping = aes(colour = "Positive catch"), show.legend = "point") +
  geom_sf(data = setdet_cod_sf_SR[[1]] |> filter(year == 11) |> filter(n == 0), mapping = aes(colour = "Zero catch"), show.legend = "point") +
  geom_sf(data = blocked_strat, fill = alpha("yellow",0.1)) +
  scale_colour_manual(values = c("Positive catch" = "blue", "Zero catch" = "orange"))+
  theme_bw()+
  coord_sf(expand = FALSE) +
  labs(color = 'Sample') +
  theme(legend.position = "right")+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey90"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(title = "d) Strata removal")


figure2 <- ggarrange(base, r30, b30, sr, ncol = 4, nrow = 1, legend = "bottom", common.legend = TRUE)
figure2

ggsave("figure2_samples.pdf", plot = figure2, width = 11, height = 3, units = "in", dpi = 500, bg = "white")


######################### Figure 3: Example time series of estimated abundance indices

colour_pal <-  c("Bootstrapped" = "#FDBF6F",
                 "DG" = "#1F78B4",
                 "DG + Depth" = "#A6CEE3",
                 "NB" = "#33A02C",
                 "NB + Depth" =  "#B2DF8A",
                 "TW" = "#E31A1C",
                 "TW + Depth" = "#FB9A99",
                 "Bootstrapped" = "#FDBF6F")

a <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Cod-like") |>
              filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), size = 1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 1 | pop == 2 | pop == 3)|>
                filter(species == "Cod-like") |>
                filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
              aes(x =  year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1) +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Cod-like"), aes(year, true), linewidth = 1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator") +
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("1" = "Population 1",
                                                                       "2" = "Population 2",
                                                                       "3" = "Population 3"))) +
  scale_fill_manual(values = colour_pal[1:3]) +
  scale_colour_manual(values = colour_pal[1:3]) +
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Cod-like")

b <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Yellowtail-like") |>
              filter(type == "Bootstrapped" | type == "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), size = 1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 1 | pop == 2 | pop == 3)|>
                filter(species == "Yellowtail-like") |>
                filter(type == "Bootstrapped" | type == "DG" | type == "DG + Depth"),
              aes(x =  year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1)+
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Yellowtail-like"), aes(year, true), linewidth = 1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator") +
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("1" = "Population 1",
                                                                       "2" = "Population 2",
                                                                       "3" = "Population 3"))) +
  scale_fill_manual(values = colour_pal[1:3]) +
  scale_colour_manual(values = colour_pal[1:3]) +
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 14),
        strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Yellowtail-like")

figure3 <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure3

ggsave("figure3.pdf", plot = figure3, width = 10, height = 12, units = "in", dpi = 500, bg = "white")


#### Figure 4: Distributions of root mean squared log error (RMSLE) and mean relative error (MRE)

index_all_scenarios_accuracy <-
  index_all_scenarios |>
  filter(year > 10) |>
  ungroup() |>
  filter(type !=  "Bootstrapped") |>
  filter(type !=  "Design-based" | scenario !=  "Strata removal") |>
  group_by(pop, type, scenario, species) |>
  summarise(MRE = mean((log(N) - log(true)) / log(true)),
            RMSLE =  sqrt(mean((log(N) - log(true))^2)))

index_all_scenarios_accuracy_long <- index_all_scenarios_accuracy |>
  pivot_longer(-c(type, pop, scenario, species))

result_table <- index_all_scenarios_accuracy_long |>
  group_by(name, type, species, scenario) |>
  summarise(m = mean(value))

result_table$type <- factor(result_table$type, levels = c("TW + Depth",
                                                        "TW",
                                                        "NB + Depth",
                                                        "NB",
                                                        "DG + Depth",
                                                        "DG",
                                                        "Design-based"))

index_all_scenarios_accuracy_long$type <- factor(index_all_scenarios_accuracy_long$type , levels = c("TW + Depth",
                                                                                                   "TW",
                                                                                                   "NB + Depth",
                                                                                                   "NB",
                                                                                                   "DG + Depth",
                                                                                                   "DG",
                                                                                                   "Design-based"))
a <- result_table |>
  filter(name == "RMSLE") |>
  ggplot(aes(x =  m, y =  type)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(data = index_all_scenarios_accuracy_long |> filter(name == "RMSLE"),
             aes(x =  value, y =  type), col = "grey", alpha = 0.5) +
  geom_point(aes(colour = type), size = 3) +
  facet_grid(species ~ scenario) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position = NULL)+
  labs(x = "RMSLE", y = "", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Root mean square log error")

b <- result_table |>
  filter(name == "MRE") |>
  ggplot(aes(x =  m, y =  type)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(data = index_all_scenarios_accuracy_long |> filter(name == "MRE"),
             aes(x =  value, y =  type), col = "grey", alpha = 0.5) +
  geom_point(aes(colour = type), size = 3) +
  facet_grid(species ~ scenario) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position =  NULL)+
  labs(x = "MRE", y = "", colour = "Estimator", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Mean relative error")



figure4_means <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "none")

figure4_means

ggsave("figure4_means.pdf", plot = figure4_means, width = 10, height = 12, units = "in", dpi = 500, bg = "white")

################### Figure 5: Coverage of the 95% confidence intervals

CI_coverage <- index_all_scenarios |>
  filter(type !=  "Design-based") |>
  filter(type !=  "Bootstrapped" | scenario !=  "Strata removal") |>
  mutate(covered = lwr < true & upr > true) |>
  group_by(type, scenario, species) |>
  summarise(mc = mean(covered)) |>
  group_by(scenario) |>
  arrange(mc)

CI_coverage |>
  group_by(species, scenario) |>
  top_n(1, mc)

a <- CI_coverage |>
  filter(species == "Cod-like") |>
  ggplot(aes(x = 0, y = mc, fill = type, group = paste(scenario, type))) +
  geom_point(position = position_dodge(width = 0.6), mapping = aes(colour = type), size = 5) +
  geom_linerange(ymin = 0, mapping = aes(ymax = mc, colour = type), position = position_dodge(width = 0.6)) +
  facet_grid(~ scenario)+
  scale_y_continuous(breaks = seq(0, 1, 0.2))+
  scale_colour_manual(values = c("Bootstrapped" = "#FDBF6F",
                               "DG" = "#1F78B4",
                               "DG + Depth" = "#A6CEE3",
                               "NB" = "#33A02C",
                               "NB + Depth" =  "#B2DF8A",
                               "TW" = "#E31A1C",
                               "TW + Depth" = "#FB9A99"))+
  labs(x = "", y = "Ratio of CI covered", colour = "Estimator", fill = "Estimator")+
  labs(y = "Ratio of CI covered") +
  labs(title = "a) Cod-like") +
  coord_cartesian(ylim = c(0,1), expand = FALSE) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = 2)+
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 14)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(fill = "grey97"))

b <- CI_coverage |>
  filter(species == "Yellowtail-like") |>
  ggplot(aes(x = 0, y = mc, fill = type, group = paste(scenario, type))) +
  geom_point(position = position_dodge(width = 0.6), mapping = aes(colour = type), size = 5) +
  geom_linerange(ymin = 0, mapping = aes(ymax = mc, colour = type), position = position_dodge(width = 0.6)) +
  facet_grid(~ scenario)+
  scale_y_continuous(breaks = seq(0, 1, 0.2))+
  scale_colour_manual(values = c("Bootstrapped" = "#FDBF6F",
                               "DG" = "#1F78B4",
                               "DG + Depth" = "#A6CEE3",
                               "NB" = "#33A02C",
                               "NB + Depth" =  "#B2DF8A",
                               "TW" = "#E31A1C",
                               "TW + Depth" = "#FB9A99"))+
  labs(x = "", y = "Ratio of CI covered", colour = "Estimator", fill = "Estimator")+
  labs(y = "Ratio of CI covered") +
  labs(title = "b) Yellowtail-like") +
  coord_cartesian(ylim = c(0,1), expand = FALSE) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = 2)+
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 14)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(fill = "grey97"))


figure5 <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure5

ggsave("figure5.pdf", plot = figure5, width = 10, height = 7, units = "in", dpi = 500, bg = "white")
