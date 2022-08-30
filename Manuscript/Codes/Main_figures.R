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



load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/cell.Rdata")
load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/strat.Rdata")

load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/setdet_cod_base.Rdata")
load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/setdet_yellowtail_base.Rdata")

load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/survey_cod_base.Rdata")
load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/survey_yellowtail_base.Rdata")

load("~/GeoStat - Rprojects/GeoStat-Survey/Manuscript/Data/index_all_scenarios.Rdata")

######################### Figure 2

### Cod-like

pop_dist_age_cod <- map(survey_cod, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod <- NULL

for(i in seq_along(pop_dist_age_cod)){
  pop_dist_cod[[i]] <- pop_dist_age_cod[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))}

pop_dist_sf_cod <- NULL
for(i in seq_along(pop_dist_cod)){
  pop_dist_sf_cod[[i]] <- merge(pop_dist_cod[[i]], cell, by=c("cell"))
  pop_dist_sf_cod[[i]] <- st_as_sf(pop_dist_sf_cod[[i]])}

st_crs(strat) = st_crs(pop_dist_sf_cod)

b_sf_cod <- NULL
for(i in seq_along(setdet_cod)){
  coordinates(setdet_cod[[i]]) = cbind(setdet_cod[[i]]$x, setdet_cod[[i]]$y)
  b_sf_cod[[i]] <- st_as_sf(setdet_cod[[i]])}

st_crs(strat) = st_crs(b_sf_cod)


cod_catch_map <- ggplot() +
  geom_raster(data = pop_dist_sf_cod[[15]] |> filter(year == 11 | year == 12 | year == 13), mapping=aes(x, y, fill = log(I_all_age))) +
  geom_sf(data = strat, mapping=aes(), fill = NA, alpha = 0.4, colour = "grey20", size = 1) +
  geom_sf(data = b_sf_cod[[15]] |> filter(year == 11 | year == 12 | year == 13), colour = "black", shape = 20, size=1) +
  geom_sf(data = b_sf_cod[[15]] |> filter(year == 11 | year == 12 | year == 13) |> filter(n == 0), colour = "white", shape = 20, size=1, inherit.aes = FALSE) +
  scale_fill_viridis_c(limits = c(0,15))+
  facet_wrap(~year, labeller = labeller(year =
                                          c("11" = "Year 11",
                                            "12" = "Year 12",
                                            "13" = "Year 13")))+
  theme_bw()+
  labs(x = "", y = "", colour = "", fill = "Actual abundance", size = "")+
  theme(legend.position = "right")+
  theme(text = element_text(size = 20),
        strip.background =element_rect(fill="grey97"))+
  coord_sf(expand = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing = unit(1, "lines"))+
  labs(title = "a) Cod-like")


### Yellowtail-like

pop_dist_age_yellowtail <- map(survey_yellowtail, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_yellowtail <- NULL

for(i in seq_along(pop_dist_age_yellowtail)){
  pop_dist_yellowtail[[i]] <- pop_dist_age_yellowtail[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))}

pop_dist_sf_yellowtail <- NULL
for(i in seq_along(pop_dist_yellowtail)){
  pop_dist_sf_yellowtail[[i]] <- merge(pop_dist_yellowtail[[i]], cell, by=c("cell"))
  pop_dist_sf_yellowtail[[i]] <- st_as_sf(pop_dist_sf_yellowtail[[i]])}

st_crs(strat) = st_crs(pop_dist_sf_yellowtail)

b_sf_yellowtail <- NULL
for(i in seq_along(setdet_yellowtail)){
  coordinates(setdet_yellowtail[[i]]) = cbind(setdet_yellowtail[[i]]$x, setdet_yellowtail[[i]]$y)
  b_sf_yellowtail[[i]] <- st_as_sf(setdet_yellowtail[[i]])}

st_crs(strat) = st_crs(b_sf_yellowtail)


yellowtail_catch_map <- ggplot() +
  geom_raster(data = pop_dist_sf_yellowtail[[11]] |> filter(year == 11 | year == 12 | year == 13), mapping=aes(x, y, fill = log(I_all_age))) +
  geom_sf(data = strat, mapping=aes(), fill = NA, alpha = 0.4, colour = "grey20", size = 1) +
  geom_sf(data = b_sf_yellowtail[[11]] |> filter(year == 11 | year == 12 | year == 13), colour = "black", shape = 20, size=1) +
  geom_sf(data = b_sf_yellowtail[[11]] |> filter(year == 11 | year == 12 | year == 13) |> filter(n == 0), colour = "white", shape = 20, size=1, inherit.aes = FALSE) +
  scale_fill_viridis_c(limits = c(0,15))+
  facet_wrap(~year, labeller = labeller(year =
                                          c("11" = "Year 11",
                                            "12" = "Year 12",
                                            "13" = "Year 13")))+
  labs(x = "", y = "", colour = "", fill = "Actual abundance", size = "")+
  theme_bw()+
  theme(legend.position = "right")+
  theme(text = element_text(size = 20),
        strip.background =element_rect(fill="grey97"))+
  coord_sf(expand = FALSE)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing = unit(1, "lines"))+
  labs(title = "b) Yellowtail-like")

figure2 <- ggarrange(cod_catch_map, yellowtail_catch_map, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure2

ggsave("figure2.png", plot = figure2, width = 10, height = 8, units = "in", dpi = 500)


######################### Figure 3

colour_pal <-  c("Bootstrapped" = "#FDBF6F",
                 "DG" = "#1F78B4",
                 "DG + Depth" = "#A6CEE3",
                 "NB" = "#33A02C",
                 "NB + Depth" =  "#B2DF8A",
                 "TW" = "#E31A1C",
                 "TW + Depth" = "#FB9A99",
                 "Bootstrapped" = "#FDBF6F")


## Abundance (raw) graphs

a <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Cod-like") |>
              filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), size=1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 1 | pop == 2 | pop == 3)|>
                filter(species == "Cod-like") |>
                filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
              aes(x= year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1)+
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Cod-like"), aes(year, true), size=1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator")+
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("1" = "Population 1",
                                                                       "2" = "Population 2",
                                                                       "3" = "Population 3")))+
  scale_fill_manual(values = colour_pal[1:3])+
  scale_colour_manual(values = colour_pal[1:3])+
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 20),
        strip.background =element_rect(fill="grey97"))+
  labs(title = "a) Cod-like")





b <- ggplot() +
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Yellowtail-like") |>
              filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
            aes(year, N, colour = type), size=1) +
  geom_ribbon(data = index_all_scenarios |>
                filter(pop == 1 | pop == 2 | pop == 3)|>
                filter(species == "Yellowtail-like") |>
                filter(type == "Bootstrapped" | type ==  "DG" | type == "DG + Depth"),
              aes(x= year, ymin = lwr, ymax = upr, fill = type), alpha = 0.1)+
  geom_line(data = index_all_scenarios |>
              filter(pop == 1 | pop == 2 | pop == 3)|>
              filter(species == "Yellowtail-like"), aes(year, true), size=1) +
  labs(x = "Year", y = "N", colour = "Estimator", fill = "Estimator", size = "Estimator")+
  facet_grid(pop ~scenario, scales = "free_y", labeller = labeller(pop =
                                                                     c("1" = "Population 1",
                                                                       "2" = "Population 2",
                                                                       "3" = "Population 3")))+
  scale_fill_manual(values = colour_pal[1:3])+
  scale_colour_manual(values = colour_pal[1:3])+
  theme_bw()+
  scale_y_log10()+
  theme(text = element_text(size = 20),
        strip.background =element_rect(fill="grey97"))+
  labs(title = "a) Yellowtail-like")

### Plot
figure3 <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure3

ggsave("figure3.png", plot = figure3, width = 12, height = 12, units = "in", dpi = 500)



############################### Figure 4


colour_pal <-  c("DG" = "#1F78B4",
                 "DG + Depth" = "#A6CEE3",
                 "NB" = "#33A02C",
                 "NB + Depth" =  "#B2DF8A",
                 "TW" = "#E31A1C",
                 "TW + Depth" = "#FB9A99",
                 "Design-based" = "#FDBF6F")


#### Figure 4

index_all_scenarios_accuracy <-
  index_all_scenarios |>
  filter(year > 10) |>
  ungroup() |>
  filter(type != "Bootstrapped") |>
  group_by(pop, type, scenario, species) |>
  summarise(E = mean(N - true), #bias
            MRE = mean((N - true)/ true), #relative bias ( = PE - 1) # when close to 0, good
            PE = mean(N/true), # relative bias, # when close to 1, good
            RMSE = sqrt(mean((N - true )^2)), #Sqrt mean squared error, # when close to 0, good
            MAPE = 100*mean(abs((N-true)/true)), # percent average absolute error (bias)
            MPE = 100*mean((N - true)/true),
            RMSLE= sqrt(mean((log(N) - log(true))^2))) # percent average error (bias)

index_all_scenarios_accuracy_long <- index_all_scenarios_accuracy |>
  pivot_longer(-c(type, pop, scenario, species))

library(facetscales)

scales_y <- list(
  `MRE` = scale_y_continuous(limits = c(-1, 2), breaks = seq(-1,2,1)),
  `RMSLE` = scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1))
)

a <- index_all_scenarios_accuracy_long |>
  filter(name == "RMSLE" | name == "MRE") |>
  filter(species == "Cod-like") |>
  ggplot(aes(as.factor(type), value, group = type)) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_violin(aes(colour = type), scale = "width")+
  geom_point(aes(colour = type),
             position = position_jitter(width = 0.1, height = 0), alpha=0.2)+
  stat_summary(fun=mean, geom="point", shape=20, size=3)+
  facet_grid_sc(name ~ scenario, scales = list(y = scales_y)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="none")+
  labs(x = "", y = "Performance metrics", colour = "Estimator", fill = "Estimator")+
  scale_colour_manual(values=colour_pal)+
  theme(text = element_text(size = 20))+
  labs(title = "a) Cod-like")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.background =element_rect(fill="grey97"))

b <- index_all_scenarios_accuracy_long |>
  filter(name == "RMSLE" | name == "MRE") |>
  filter(species == "Yellowtail-like") |>
  ggplot(aes(as.factor(type), value, group = type)) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_violin(aes(colour = type), scale = "width")+
  geom_point(aes(colour = type),
             position = position_jitter(width = 0.1, height = 0), alpha=0.2)+
  stat_summary(fun=mean, geom="point", shape=20, size=3)+
  #stat_summary(aes(label=round(..y.., 2)), fun=mean, geom="text", size=2) +
  facet_grid_sc(name ~ scenario, scales = list(y = scales_y)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position="bottom")+
  labs(x = "", y = "Performance metrics", colour = "Estimator", fill = "Estimator")+
  scale_colour_manual(values=colour_pal)+
  theme(text = element_text(size = 20))+
  labs(title = "b) Yellowtail-like")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.background =element_rect(fill="grey97"))

figure4 <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)

figure4

ggsave("figure4.png", plot = figure4, width = 12, height = 12, units = "in", dpi = 500)


################### Figure 5

colour_pal <-  c(#"Design-based" = "#FDBF6F",
  "Bootstrapped" = "#FDBF6F",
  "DG" = "#1F78B4",
  "DG + Depth" = "#A6CEE3",
  "NB" = "#33A02C",
  "NB + Depth" =  "#B2DF8A",
  "TW" = "#E31A1C",
  "TW + Depth" = "#FB9A99")

### CI coverage

CI_coverage <- index_all_scenarios |>
  filter(type != "Design-based") |>
  mutate(covered = lwr < true & upr > true) |>
  group_by(type, scenario, species) |>
  summarise(mc = mean(covered)) |>
  group_by(scenario) |>
  arrange(mc)

a <- CI_coverage |>
  filter(species == "Cod-like") |>
  ggplot(aes(x = 0, y = mc, fill = type, group=paste(scenario, type))) +
  geom_point(position = position_dodge(width = 0.6), mapping = aes(colour=type), size=5) +
  #geom_col(position = "dodge") +
  geom_linerange(ymin = 0, mapping = aes(ymax = mc, colour = type), position = position_dodge(width = 0.6)) +
  labs(y = "Ratio of CI covered", x = "Scenario") +
  theme_bw()+
  theme(legend.position="bottom")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(text = element_text(size = 20))+
  facet_grid(~ scenario)+
  labs(x = "", y = "Ratio of CI covered", colour = "Estimator", fill = "Estimator")+
  scale_y_continuous(breaks=seq(0,1, 0.1))+
  scale_colour_manual(values=colour_pal)+
  coord_cartesian(ylim=c(0,1), expand = FALSE)+
  labs(title = "a) Cod-like")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.background =element_rect(fill="grey97"))


b <- CI_coverage |>
  filter(species == "Yellowtail-like") |>
  ggplot(aes(x = 0, y = mc, fill = type, group=paste(scenario, type))) +
  geom_point(position = position_dodge(width = 0.6), mapping = aes(colour=type), size=5) +
  #geom_col(position = "dodge") +
  geom_linerange(ymin = 0, mapping = aes(ymax = mc, colour = type), position = position_dodge(width = 0.6)) +
  labs(y = "Ratio of CI covered", x = "Scenario") +
  theme_bw()+
  theme(legend.position="bottom")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(text = element_text(size = 20))+
  facet_grid(~ scenario)+
  labs(x = "", y = "Ratio of CI covered", colour = "Estimator", fill = "Estimator")+
  scale_y_continuous(breaks=seq(0,1, 0.1))+
  scale_colour_manual(values=colour_pal)+
  coord_cartesian(ylim=c(0,1), expand = FALSE)+
  labs(title = "b) Yellowtail-like")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        strip.background =element_rect(fill="grey97"))


figure5 <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)
figure5

ggsave("figure5.png", plot = figure5, width = 12, height = 8, units = "in", dpi = 500)
