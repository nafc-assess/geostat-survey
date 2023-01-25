##Appendix 1 – Additional figures

## Loading libraries
library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggpubr)
library(sp)
library(sf)
library(sdmTMB)
library(purrr)
library(SimSurvey)
library(raster)
library(here)

## Loading the data
load(here("Data", "cell.Rdata"))
load(here("Data", "strat.Rdata"))
load(here("Data", "survey_cod_base.Rdata"))
load(here("Data", "sdm_data_cod_SR.Rdata"))
load(here("Data", "sdm_newdata_cod.Rdata"))
load(here("Data", "mesh_sdm_cod_SR.Rdata"))
load(here("Data", "setdet_yellowtail_SR.Rdata"))
load(here("Data", "survey_yellowtail_base.Rdata"))
load(here("Data", "sdm_data_yellowtail_SR.Rdata"))
load(here("Data", "sdm_newdata_yellowtail.Rdata"))
load(here("Data", "mesh_sdm_yellowtail_SR.Rdata"))
load(here("Data", "setdet_cod_SR.Rdata"))
load(here("Data", "index_all_scenarios.Rdata"))

######################## Figure S1: Depth profile + species depth preferences

grid <- make_grid(x_range = c(-150, 150),
                  y_range = c(-150, 150),
                  res = c(1, 1), # increased resolution to make smoother lines
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1000),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 40),
                  strat_splits = 2,
                  method = "spline")
cod_depth <- sim_parabola(mu = log(160),
                          sigma = 0.5, log_space = TRUE, plot=FALSE)
yt_depth <- sim_parabola(mu = log(90),
                         sigma = 0.3,
                         #sigma_right = 0.20,
                         log_space = TRUE)

xyz <- as.data.frame(grid) |> cbind(coordinates(grid))

(profile <- ggplot(data = xyz) +
    geom_ribbon(aes(x = x, ymin = max(depth), ymax = depth), fill = "grey60", color = "grey50") +
    scale_y_reverse(expand = c(0, 0), limits = c(max(xyz$depth), 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    ylab("Depth (m)") +
    theme_bw()) +
  theme(text = element_text(size = 14))

depth <- sort(unique(xyz$depth))
yt_effect <- yt_depth(data.frame(age = 1, depth = depth))
cod_effect <- cod_depth(data.frame(age = 1, depth = depth))
depth_effect <- data.frame(depth = c(depth, depth),
                           species = c(rep("Cod-like", length(depth)),
                                       rep("Yellowtail-like", length(depth))),
                           effect = c(cod_effect, yt_effect))

(parabola <- ggplot(data = depth_effect) +
    geom_path(aes(x = exp(effect), y = depth, color = species), linewidth = 1) +
    scale_y_reverse(expand = c(0, 0), limits = c(max(xyz$depth), 0)) +
    scale_color_brewer(palette = "Set1", name = "Species") +
    xlab("Depth effect") +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size = 14)))

library(patchwork)

FigureS_depth <- profile + parabola +
  plot_layout(widths = c(2, 1))

FigureS_depth

ggsave("FigureS_depth.pdf", plot = FigureS_depth, width = 8, height = 6, units = "in", dpi = 500, bg = "white")

######################## Figure S2 Cod-like species maps

## Delta-gamma model

sdm_DG_IID_dg_cod <- sdmTMB(count ~ 0 + as.factor(year),
                                 data = sdm_data_cod_SR[[1]],
                                 mesh = mesh_sdm_cod_SR[[1]],
                                 offset = sdm_data_cod_SR[[1]]$offset,
                                 time = "year",
                                 family = delta_gamma(),
                                 spatial = TRUE,
                                 spatiotemporal = list("off", "IID"),
                                 priors = sdmTMBpriors(
                                   matern_s = pc_matern(range_gt = 75, sigma_lt = 7.5),
                                   matern_st = pc_matern(range_gt = 75, sigma_lt = 7.5)),
                                 share_range = FALSE,
                                 control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

sdm_DG_IID_dg_pred_cod  <- predict(sdm_DG_IID_dg_cod,
                                        newdata = sdm_newdata_cod,
                                        return_tmb_object = TRUE)

dg_table_cod <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(pop == 1) |>
  filter(type == "DG") |>
  filter(scenario == "Strata removal") |>
  filter(species == "Cod-like") |>
  dplyr::select(year, N)

## Delta-gamma + depth model

formula1 = count ~ 0 + as.factor(year)
formula2 = count ~ 0 + as.factor(year) + poly(log(depth), 2)
sdm_DG_IID_dg_depth_cod <- sdmTMB(list(formula1, formula2),
                                       data = sdm_data_cod_SR[[1]],
                                       mesh = mesh_sdm_cod_SR[[1]],
                                       offset = sdm_data_cod_SR[[1]]$offset,
                                       time = "year",
                                       family = delta_gamma(),
                                       spatial = TRUE,
                                       spatiotemporal = list("off", "IID"),
                                       priors = sdmTMBpriors(
                                         matern_s = pc_matern(range_gt = 75, sigma_lt = 12.5),
                                         matern_st = pc_matern(range_gt = 75, sigma_lt = 12.5)),
                                       share_range = FALSE,
                                       control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

sdm_DG_IID_dg_pred_depth_cod <-  predict(sdm_DG_IID_dg_depth_cod,
                                              newdata = sdm_newdata_cod,
                                              return_tmb_object = TRUE)

dg_table_depth_cod <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(pop == 1) |>
  filter(type == "DG + Depth") |>
  filter(scenario == "Strata removal") |>
  filter(species == "Cod-like") |>
  dplyr::select(year, N)

## Actual abundance distribution

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

true_table_cod <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(species == "Cod-like" & pop == 1) |>
  dplyr::select(year, true)

## Design-based

sumYst <- function(data, i = seq_len(nrow(data))) {
  data[i,] |>
    ### stratum level
    group_by(year, strat, strat_area) |>
    summarise(meanYh = mean(n), tow_area = mean(tow_area), .groups = "drop_last") |>
    mutate(Nh = strat_area/(tow_area)) |>
    group_by(year, strat) |>
    mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh, sumYst = N * WhmeanYh)
}

Years <- data.frame(year = c(1:20))
strat <- merge(strat, Years, all=TRUE)

blocked_strat <- strat |>
  filter(strat %in% c(2,3,4,5,17,18,19,20))

cod_sr <- sumYst(setdet_cod_SR[[1]])
strat_cod_sr <- full_join(strat, cod_sr)

design_table_cod <- strat_cod_sr |>
  filter(year > 10 & year < 16) |>
  group_by(year) |>
  summarise(Nsum = sum(sumYst, na.rm = TRUE))

a <- pop_dist_sf_cod[[1]] |>
  filter(year > 10 & year < 16) |>
  ggplot(aes(x, y, fill = I_all_age)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=true_table_cod,aes(x=0, y=-120, label=paste0(round(true/1e+06, 1), " * 1e+06")), col= "black", size = 4, inherit.aes = FALSE)+

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
  labs(title = "a) Actual distribution") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

a2 <- strat_cod_sr |>
  filter(year > 10 & year < 16) |>
  ggplot() +
  geom_sf(aes(fill = sumYst)) +
  #geom_sf_label(aes(label = strat)) +
  geom_sf(data=blocked_strat|>
            filter(year > 10 & year < 16), colour="red", fill=NA) +
  coord_sf(expand=FALSE) +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=design_table_cod,aes(x=0, y=-120, label=paste0(round(Nsum/1e+06, 1), " * 1e+06")), col= "black", size = 4, inherit.aes = FALSE)+
  #coord_fixed(expand = FALSE) +
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
  labs(title = "b) Distribution from design-based")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

b <- sdm_DG_IID_dg_pred_cod[[1]] |>
  filter(year > 10 & year < 16) |>
  ggplot(aes(x, y, fill = plogis(est1) * exp(est2))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10")+
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=dg_table_cod,aes(x=0, y=-120, label=paste0(round(N/1e+06, 1), " * 1e+06")), col= "black", size = 4, inherit.aes = FALSE) +
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
  labs(title = "c) Distribution from DG model")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

c <- sdm_DG_IID_dg_pred_depth_cod[[1]] |>
  filter(year > 10 & year < 16) |>
  ggplot(aes(x, y, fill = plogis(est1) * exp(est2))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=dg_table_depth_cod,aes(x=0, y=-120, label=paste0(round(N/1e+06, 1), " * 1e+06")), col= "black", size = 4, inherit.aes = FALSE) +
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
  labs(title = "d) Distribution from DG + depth effect model")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

FigureS_cod_maps <- ggarrange(a, a2, b, c, ncol = 1, nrow = 4, legend = "bottom", common.legend = TRUE)

FigureS_cod_maps

ggsave("FigureS_cod_maps.pdf", plot = FigureS_cod_maps, width = 10, height = 12, units = "in", dpi = 500, bg="white")

######################## Figure S3 Yellowtail-like species maps

## Delta-gamma model

sdm_DG_IID_dg_yt <- sdmTMB(count ~ 0 + as.factor(year),
                           data = sdm_data_yellowtail_SR[[1]],
                           mesh = mesh_sdm_yellowtail_SR[[1]],
                           offset = sdm_data_yellowtail_SR[[1]]$offset,
                           time = "year",
                           family = delta_gamma(),
                           spatial = TRUE,
                           spatiotemporal = list("off", "IID"),
                           priors = sdmTMBpriors(
                             matern_s = pc_matern(range_gt = 125, sigma_lt = 12.5),
                             matern_st = pc_matern(range_gt = 125, sigma_lt = 12.5)),
                           share_range = FALSE,
                           control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

sdm_DG_IID_dg_pred_yt <- predict(sdm_DG_IID_dg_yt,
                                 newdata = sdm_newdata_yellowtail,
                                 return_tmb_object = TRUE)

dg_table_yt <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(pop == 1) |>
  filter(type == "DG") |>
  filter(scenario == "Strata removal") |>
  filter(species == "Yellowtail-like") |>
  dplyr::select(year, N)

## Delta-gamma + depth model

formula1 = count ~ 0 + as.factor(year)
formula2 = count ~ 0 + as.factor(year) + poly(log(depth), 2)

sdm_DG_IID_dg_depth_yt <- sdmTMB(list(formula1, formula2),
                                 data = sdm_data_yellowtail_SR[[1]],
                                 mesh = mesh_sdm_yellowtail_SR[[1]],
                                 offset = sdm_data_yellowtail_SR[[1]]$offset,
                                 time = "year",
                                 family = delta_gamma(),
                                 spatial = TRUE,
                                 spatiotemporal = list("off", "IID"),
                                 priors = sdmTMBpriors(
                                   matern_s = pc_matern(range_gt = 125, sigma_lt = 12.5),
                                   matern_st = pc_matern(range_gt = 125, sigma_lt = 12.5)),
                                 share_range = FALSE,
                                 control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

sdm_DG_IID_dg_pred_depth_yt <-  predict(sdm_DG_IID_dg_depth_yt,
                                        newdata = sdm_newdata_yellowtail,
                                        return_tmb_object = TRUE)

dg_table_depth_yt <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(pop == 1) |>
  filter(type == "DG + Depth") |>
  filter(scenario == "Strata removal") |>
  filter(species == "Yellowtail-like") |>
  dplyr::select(year, N)

## Actual abundance distribution

pop_dist_age_yt <- map(survey_yellowtail, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})
gc()

pop_dist_yt <- NULL
for(i in seq_along(pop_dist_age_yt)){
  pop_dist_yt[[i]] <- pop_dist_age_yt[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))}

pop_dist_sf_yt <- NULL
for(i in seq_along(pop_dist_yt)){
  pop_dist_sf_yt[[i]] <- merge(pop_dist_yt[[i]], cell, by=c("cell"))
  pop_dist_sf_yt[[i]] <- st_as_sf(pop_dist_sf_yt[[i]])}

true_table_yt <- index_all_scenarios |>
  filter(year > 10 & year < 16) |>
  filter(species == "Yellowtail-like" & pop == 1) |>
  dplyr::select(year, true)

## Design-based

sumYst <- function(data, i = seq_len(nrow(data))) {
  data[i,] |>
    ### stratum level
    group_by(year, strat, strat_area) |>
    summarise(meanYh = mean(n), tow_area = mean(tow_area), .groups = "drop_last") |>
    mutate(Nh = strat_area/(tow_area)) |>
    group_by(year, strat) |>
    mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh, sumYst = N * WhmeanYh)
}

Years <- data.frame(year = c(1:20))
strat <- merge(strat, Years, all=TRUE)

blocked_strat <- strat |>
  filter(strat %in% c(2,3,4,5,17,18,19,20))

yt_sr <- sumYst(setdet_yellowtail_SR[[1]])
strat_yt_sr <- full_join(strat, yt_sr)

design_table_yt <- strat_yt_sr |>
  filter(year > 10 & year < 16) |>
  group_by(year) |>
  summarise(Nsum = sum(sumYst, na.rm = TRUE))

a <- pop_dist_sf_yt[[1]] |>
  filter(year > 10 & year < 16) |>
ggplot(aes(x, y, fill = I_all_age)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=true_table_yt,aes(x=0, y=-120, label=paste0(round(true/1e+08, 1), " * 1e+08")), col= "black", size = 4, inherit.aes = FALSE)+
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
  labs(title = "a) Actual distribution") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

a2 <- strat_yt_sr |>
  filter(year > 10 & year < 16) |>
  ggplot() +
  geom_sf(aes(fill = sumYst)) +
  #geom_sf_label(aes(label = strat)) +
  geom_sf(data=blocked_strat|>
            filter(year > 10 & year < 16), colour="red", fill=NA) +
  coord_sf(expand=FALSE) +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=design_table_yt,aes(x=0, y=-120, label=paste0(round(Nsum/1e+08, 1), " * 1e+08")), col= "black", size = 4, inherit.aes = FALSE)+
  #coord_fixed(expand = FALSE) +
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
  labs(title = "b) Distribution from design-based")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

b <- sdm_DG_IID_dg_pred_yt[[1]] |>
  filter(year > 10 & year < 16) |>
ggplot(aes(x, y, fill = plogis(est1) * exp(est2))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10")+
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
  geom_text(data=dg_table_yt,aes(x=0, y=-120, label=paste0(round(N/1e+08, 1), " * 1e+08")), col= "black", size = 4, inherit.aes = FALSE) +
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
  labs(title = "c) Distribution from DG model")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

c <- sdm_DG_IID_dg_pred_depth_yt[[1]] |>
  filter(year > 10 & year < 16) |>
ggplot(aes(x, y, fill = plogis(est1) * exp(est2))) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  labs(fill="Abundance")+
  facet_wrap(~year, nrow=1, labeller = labeller(year =
                                                  c("11" = "Year 11",
                                                    "12" = "Year 12",
                                                    "13" = "Year 13",
                                                    "14" = "Year 14",
                                                    "15" = "Year 15"))) +
geom_text(data=dg_table_depth_yt,aes(x=0, y=-120, label=paste0(round(N/1e+08, 1), " * 1e+08")), col= "black", size = 4, inherit.aes = FALSE) +
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
  labs(title = "d) Distribution from DG + depth effect model")+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

FigureS_yt_maps <- ggarrange(a, a2, b, c, ncol = 1, nrow = 4, legend = "bottom", common.legend = TRUE)

FigureS_yt_maps

ggsave("FigureS_yt_maps.pdf", plot = FigureS_yt_maps, width = 10, height = 12, units = "in", dpi = 500, bg="white")


######################## Figure S4 Yellowtail-like species maps

data_reg <- index_all_scenarios |>
  arrange(species, pop, type, scenario, year) |>
  group_by(species, pop, type, scenario) |>
  mutate(re = (log(N) - log(true))/ log(true))

all_regress <-  data_reg |>
  group_by(species, pop, sim, type, scenario) %>%
  do(mod = lm(re ~ year, .)) %>% ungroup()

summary_lm <-all_regress %>% mutate(tidy = map(mod, conf.int = TRUE, broom::tidy),
                                    glance = map(mod, broom::glance),
                                    augment = map(mod, broom::augment),
                                    rsq = glance |> map_dbl('r.squared'),
                                    slope = tidy |> map_dbl(function(x) x$estimate[2]),
                                    c_low = tidy |>  map_dbl(function(x) x$conf.low[2]),
                                    c_high = tidy |>  map_dbl(function(x) x$conf.high[2]))

summary_lm$type <- factor(summary_lm$type, levels=c("TW + Depth",
                                                    "TW",
                                                    "NB + Depth",
                                                    "NB",
                                                    "DG + Depth",
                                                    "DG",
                                                    "Design-based"))
result_table <- summary_lm |>
  group_by(type, species, scenario) |>
  summarise(m=mean(slope))

result_table$type <- factor(result_table$type, levels=c("TW + Depth",
                                                        "TW",
                                                        "NB + Depth",
                                                        "NB",
                                                        "DG + Depth",
                                                        "DG",
                                                        "Design-based"))

trend_bias_cod <- result_table |>
  filter(!type=="Bootstrapped") |>
  filter(type != "Design-based" | scenario != "Strata removal") |>
  filter(species == "Cod-like")|>
  ggplot(aes(x= m, y= type)) +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_point(data=summary_lm |>
               filter(!type=="Bootstrapped") |>
               filter(type != "Design-based" | scenario != "Strata removal") |>
              filter(species == "Cod-like"),
             aes(x= slope, y= type), col = "grey", alpha=0.5)+
  geom_point(aes(colour = type), size = 3)+
  facet_grid( ~ scenario) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position="bottom")+
  labs(x = "Slope of ln(RE) ~ year", y = "", colour = "Estimator", fill = "Estimator")+
  theme(text = element_text(size = 14))+
  theme(strip.background =element_rect(fill="grey97"))+
  labs(title = "a) Cod-like species")

trend_bias_yt <- result_table |>
  filter(!type=="Bootstrapped") |>
  filter(type != "Design-based" | scenario != "Strata removal") |>
  filter(species == "Yellowtail-like") |>
  ggplot(aes(x= m, y= type)) +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_point(data=summary_lm |>
               filter(!type=="Bootstrapped") |>
               filter(type != "Design-based" | scenario != "Strata removal") |>
               filter(species == "Yellowtail-like"),
             aes(x= slope, y= type), col = "grey", alpha=0.5)+
  geom_point(aes(colour = type), size = 3)+
  facet_grid( ~ scenario) +
  theme_bw() +
  scale_color_manual(values = c("DG" = "#1F78B4",
                                "DG + Depth" = "#A6CEE3",
                                "NB" = "#33A02C",
                                "NB + Depth" =  "#B2DF8A",
                                "TW" = "#E31A1C",
                                "TW + Depth" = "#FB9A99",
                                "Design-based" = "#FDBF6F")) +
  theme(legend.position="bottom")+
  labs(x = "Slope of ln(RE) ~ year", y = "", colour = "Estimator", fill = "Estimator")+
  theme(text = element_text(size = 14))+
  theme(strip.background =element_rect(fill="grey97"))+
  labs(title = "b) Yellowtail-like species")

FigureS_trend_bias <- ggarrange(trend_bias_cod, trend_bias_yt, ncol = 1, nrow = 2, legend = "none")

FigureS_trend_bias

ggsave("FigureS_trend_bias.pdf", plot = FigureS_trend_bias, width = 10, height = 6, units = "in", dpi = 500)

######################## Figure S5 delta-AIC

all_AIC <-
  index_all_scenarios |>
  #filter(year > 10) |>
  ungroup() |>
  filter(type != "Bootstrapped") |>
  filter(type != "Design-based") |>
  group_by(pop, type, scenario, species) |>
  summarise(AIC = mean(AIC))

best_models <- map_dfr(unique(all_AIC$pop),function(i) {
  all_AIC_dAIC <- all_AIC |>
    filter(pop == i) |>
    group_by(species, scenario) |>
    slice_min(AIC) |>
    mutate(mAIC=AIC)
  all_AIC_dAIC$AIC <- NULL
  return(all_AIC_dAIC)
})

AIC_w_best <- merge(all_AIC, best_models, by=c("pop", "scenario", "species"))

AIC_w_best <- AIC_w_best |>
  mutate(dAIC = AIC - mAIC)

a <- AIC_w_best |>
  filter(type.x != "NB + Depth") |>
  filter(species == "Cod-like") |>
  ggplot() +
  geom_violin(aes(y = reorder(type.x, desc(dAIC)), x = dAIC, fill = type.x), colour= NA) +
  #geom_point(data = AIC_w_best |> filter(species == "Cod-like"),
  #aes(x = reorder(type.x, desc(dAIC)), y = dAIC))+
  facet_grid( ~ scenario, scales = "free_y")+
  theme_bw()+
  scale_fill_manual(values = c("DG" = "#1F78B4",
                               "DG + Depth" = "#A6CEE3",
                               "NB" = "#33A02C",
                               "NB + Depth" =  "#B2DF8A",
                               "TW" = "#E31A1C",
                               "TW + Depth" = "#FB9A99"))+
  labs(y = "") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Cod-like species")

b <- AIC_w_best |>
  filter(type.x != "NB + Depth") |>
  filter(species == "Yellowtail-like") |>
  ggplot() +
  geom_violin(aes(y = reorder(type.x, desc(dAIC)), x = dAIC, fill = type.x), colour= NA) +
  #geom_point(data = AIC_w_best |> filter(species == "Cod-like"),
  #aes(x = reorder(type.x, desc(dAIC)), y = dAIC))+
  facet_grid( ~ scenario, scales = "free_y")+
  theme_bw()+
  scale_fill_manual(values = c("DG" = "#1F78B4",
                               "DG + Depth" = "#A6CEE3",
                               "NB" = "#33A02C",
                               "NB + Depth" =  "#B2DF8A",
                               "TW" = "#E31A1C",
                               "TW + Depth" = "#FB9A99"))+
  labs(y = "") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97"))+
  labs(title = "b) Yellowtail-like species")

dAIC <- ggarrange(a, b, ncol = 1, nrow = 2, legend = "none")

dAIC

ggsave("FigureS_dAIC.pdf", plot = dAIC, width = 10, height = 8, units = "in", dpi = 500, bg = "white")

######################## Figure S6 and S7 Yellowtail-like a model results

yt_nb_sr <- sdmTMB(count ~ 0 + as.factor(year),
                   data = sdm_data_yellowtail_SR[[60]],
                   mesh = mesh_sdm_yellowtail_SR[[60]],
                   offset = sdm_data_yellowtail_SR[[60]]$offset,
                   time = "year",
                   family = nbinom2(),
                   spatial = TRUE,
                   spatiotemporal = list("IID"),
                   priors = sdmTMBpriors(
                     matern_s = pc_matern(range_gt = 125, sigma_lt = 12.5),
                     matern_st = pc_matern(range_gt = 125, sigma_lt = 12.5)),
                   share_range = FALSE,
                   control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

yt_nb_sr_pred <-  predict(yt_nb_sr,
                          newdata = sdm_newdata_yellowtail,
                          return_tmb_object = TRUE)

yt_nb_sr_depth <- sdmTMB(count ~ 0 + as.factor(year) + poly(log(depth), 2),
                         data = sdm_data_yellowtail_SR[[60]],
                         mesh = mesh_sdm_yellowtail_SR[[60]],
                         offset = sdm_data_yellowtail_SR[[60]]$offset,
                         time = "year",
                         family = nbinom2(),
                         spatial = TRUE,
                         spatiotemporal = list("IID"),
                         priors = sdmTMBpriors(
                           matern_s = pc_matern(range_gt = 125, sigma_lt = 12.5),
                           matern_st = pc_matern(range_gt = 125, sigma_lt = 12.5)),
                         share_range = FALSE,
                         control = sdmTMBcontrol(newton_loops = 1L, nlminb_loops = 1L))

yt_nb_sr_depth_pred <-  predict(yt_nb_sr_depth,
                                newdata = sdm_newdata_yellowtail,
                                return_tmb_object = TRUE)

#index_nb <- get_index(yt_nb_sr_pred, area = yt_nb_sr_pred$data$area, bias_correct = TRUE)
#index_nb_dpeth <- get_index(yt_nb_sr_depth_pred, area = yt_nb_sr_depth_pred$data$area, bias_correct = TRUE)
#ggplot() +
#geom_line(data=index_nb, aes(x=year, y=est), colour="blue")+
#geom_line(data=index_nb_dpeth, aes(x=year, y=est), colour="red")+
#theme_bw()

yt_nb_pred_long <- yt_nb_sr_pred$data |>
  pivot_longer(-c(x, y, depth, year, offset, area))

yt_nb_depth_pred_long <- yt_nb_sr_depth_pred$data |>
  pivot_longer(-c(x, y, depth, year, offset, area))

######################## Figure S6 Estimates with fix and random effects

yt_models_year12 <- rbind(yt_nb_pred_long |>
                            mutate(type = "NB") |>
                            filter(year == 12) |>
                            filter(name != "est_rf"),
                          yt_nb_depth_pred_long |>
                            mutate(type = "NB + Depth") |>
                            filter(year == 12) |>
                            filter(name != "est_rf"))

model_res <- c(
  `est` = "a) Estimate",
  `est_non_rf` = "b) Year",
  `omega_s` = "c) Spatial RF",
  `epsilon_st` = "d) Spatiotemporal RF"
)

nb <- ggplot(yt_models_year12 |> filter(type == "NB")) +
  geom_raster(aes(x, y, fill = exp(value))) +
  scale_fill_viridis_c(trans = "log10") +
  geom_sf(data=blocked_strat|>
            filter(year == 12), colour="#ffffff50", fill=NA) +
  #facet_grid(rows= vars(year), cols=vars(type))+
  facet_wrap(~ factor(name, levels=c('est','est_non_rf','omega_s', 'epsilon_st')), ncol = 4, labeller = as_labeller(model_res)) +
  theme_bw() +
  coord_sf(expand=FALSE) +
  labs(fill = "") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "Negative binomial model") +
  theme(legend.position = "bottom", legend.spacing.x = unit(0.5, 'cm'), legend.direction="horizontal")+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

model_res2 <- c(
  `est` = "e) Estimate",
  `est_non_rf` = "f) Year + Depth",
  `omega_s` = "g) Spatial RF",
  `epsilon_st` = "h) Spatiotemporal RF"
)

nb_depth <- ggplot(yt_models_year12 |> filter(type == "NB + Depth")) +
  geom_raster(aes(x, y, fill = exp(value))) +
  scale_fill_viridis_c(trans = "log10") +
  geom_sf(data=blocked_strat|>
            filter(year == 12 ), colour="#ffffff50", fill=NA) +
  coord_sf(expand=FALSE) +
  #facet_grid(rows= vars(year), cols=vars(type))+
  facet_wrap(~ factor(name, levels=c('est','est_non_rf','omega_s', 'epsilon_st')), ncol = 4, labeller = as_labeller(model_res2)) +
  theme_bw() +
  labs(fill = "") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "Negative binomial model including depth covariate") +
  theme(legend.position = "bottom", legend.spacing.x = unit(0.5, 'cm'), legend.direction="horizontal")+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

FigureS_RF <- ggarrange(nb, nb_depth, ncol = 1, nrow = 2, common.legend = TRUE, legend = "bottom")
FigureS_RF

ggsave("FigureS_RF.pdf", plot = FigureS_RF, width = 10, height = 8, units = "in", dpi = 500, bg = "white")

######################## Figure S7 Estimates with fix and random effects

a <- yt_nb_pred_long |>
  filter(year < 6 | year > 15) |>
  filter(name == "epsilon_st") |>
  ggplot() +
  geom_raster(aes(x, y, fill = exp(value))) +
  scale_fill_viridis_c(trans = "log10") +
  geom_sf(data=blocked_strat|>
            filter(year < 6 | year > 15), colour="#ffffff50", fill=NA) +
  coord_sf(expand=FALSE) +
  facet_wrap(~year, ncol=5, labeller = labeller(year =
                                                  c("1" = "Year 1",
                                                    "2" = "Year 2",
                                                    "3" = "Year 3",
                                                    "4" = "Year 4",
                                                    "5" = "Year 5",
                                                    "16" = "Year 16",
                                                    "17" = "Year 17",
                                                    "18" = "Year 18",
                                                    "19" = "Year 19",
                                                    "20" = "Year 20")))+
  theme(text = element_text(size = 14))+
  labs(fill="Spatiotemporal random effect")+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "a) Negative binomial model") +
  theme(legend.position = "bottom", legend.spacing.x = unit(0.5, 'cm'), legend.direction="horizontal")+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

b <- yt_nb_depth_pred_long |>
  filter(year < 6 | year > 15) |>
  filter(name == "epsilon_st") |>
  ggplot() +
  geom_raster(aes(x, y, fill = exp(value))) +
  scale_fill_viridis_c(trans = "log10") +
  geom_sf(data=blocked_strat|>
            filter(year < 6 | year > 15), colour="#ffffff50", fill=NA, size = 0.4, alpha = 0.3) +
  coord_sf(expand=FALSE) +
  facet_wrap(~year, ncol=5, labeller = labeller(year =
                                                  c("1" = "Year 1",
                                                    "2" = "Year 2",
                                                    "3" = "Year 3",
                                                    "4" = "Year 4",
                                                    "5" = "Year 5",
                                                    "16" = "Year 16",
                                                    "17" = "Year 17",
                                                    "18" = "Year 18",
                                                    "19" = "Year 19",
                                                    "20" = "Year 20"))) +
  theme(text = element_text(size = 14))+
  labs(fill="Spatiotemporal random effect")+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs(title = "b) Negative binomial model including depth covariate") +
  theme(legend.position = "bottom", legend.spacing.x = unit(0.5, 'cm'), legend.direction="horizontal")+
  theme(legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

FigureS_st_maps <- ggarrange(a + rremove("ylab") + rremove("xlab"),
                             b + rremove("ylab") + rremove("xlab"),
                             labels = NULL,
                             ncol = 1, nrow = 2,
                             common.legend = TRUE, legend = "bottom",
                             align = "hv",
                             font.label = list(size = 14, color = "black", face = "bold", family = NULL, position = "top"))

FigureS_st_maps

ggsave("FigureS_st_maps.pdf", plot = FigureS_st_maps, width = 10, height = 12, units = "in", dpi = 500, bg = "white")
