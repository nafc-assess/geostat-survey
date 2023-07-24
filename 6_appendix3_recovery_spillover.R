#Appendix 3 --- MPA Recovery and Spillover Effect maps

library(here)
load((here("data", "survey_yellowtail_base.Rdata")))
load((here("data", "survey_yellowtail_rec.Rdata")))
load((here("data", "survey_yellowtail_so.Rdata")))

load((here("data", "survey_cod_base.Rdata")))
load((here("data", "survey_cod_rec.Rdata")))
load((here("data", "survey_cod_so.Rdata")))


library(purrr)
library(dplyr)
library(ggplot2)
library(sf)
library(ggpubr)

##### Yellowtail maps
# Base
pop_dist_age_base_yt <- map(survey_yellowtail, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_yellowtail_base <- NULL
for(i in seq_along(pop_dist_age_base_yt)){
  pop_dist_yellowtail_base[[i]] <- pop_dist_age_base_yt[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))|>
    mutate(scenario = "Base", pop = as.numeric(i)) |>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}

# Recovery
pop_dist_age_rec_yt <- map(survey_yellowtail_rec, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_yellowtail_rec <- NULL
for(i in seq_along(pop_dist_age_rec_yt)){
  pop_dist_yellowtail_rec[[i]] <- pop_dist_age_rec_yt[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))|>
    mutate(scenario = "Recovery", pop = as.numeric(i))|>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}

# Spillover
pop_dist_age_so_yt <- map(survey_yellowtail_so, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_yellowtail_so <- NULL
for(i in seq_along(pop_dist_age_so_yt)){
  pop_dist_yellowtail_so[[i]] <- pop_dist_age_so_yt[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age = sum(I)) |>
    mutate(scenario = "Spillover", pop = as.numeric(i))|>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}


## plotting the map
combine_status2 <- rbind(pop_dist_yellowtail_base[[30]],
                         pop_dist_yellowtail_rec[[30]],
                         pop_dist_yellowtail_so[[30]])

maps <- ggplot(combine_status2|>
                 filter(year > 10), aes(x, y, fill = N_all_age)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_grid(scenario ~ year, labeller = labeller(year =
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
  labs(title = "b) Yellowtail-like") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

maps

#################

##### Cod maps

# Base
pop_dist_age_base_cod <- map(survey_cod, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod_base <- NULL
for(i in seq_along(pop_dist_age_base_cod)){
  pop_dist_cod_base[[i]] <- pop_dist_age_base_cod[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))|>
    mutate(scenario = "Base", pop = as.numeric(i)) |>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}

# Recovery
pop_dist_age_rec_cod <- map(survey_cod_rec, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod_rec <- NULL
for(i in seq_along(pop_dist_age_rec_cod)){
  pop_dist_cod_rec[[i]] <- pop_dist_age_rec_cod[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age= sum(I))|>
    mutate(scenario = "Recovery", pop = as.numeric(i))|>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}

# Spillover
pop_dist_age_so_cod <- map(survey_cod_so, function(x) {merge(x$grid_xy, x$sp_N, by=c("cell"))})

pop_dist_cod_so <- NULL
for(i in seq_along(pop_dist_age_so_cod)){
  pop_dist_cod_so[[i]] <- pop_dist_age_so_cod[[i]] |>
    group_by(cell, x, y, division, strat, year) |>
    summarise(N_all_age = sum(N), I_all_age = sum(I)) |>
    mutate(scenario = "Spillover", pop = as.numeric(i))|>
    mutate(status = "Unprotected") |>
    mutate(status = ifelse(x < 60 & x > -110 & y > -80 & y < 120, "Buffer", status)) |>
    mutate(status = ifelse(x < 50 & x > -100 & y > -70 & y < 110, "MPA", status))}

## plotting the map
combine_status <- rbind(pop_dist_cod_base[[30]],
                        pop_dist_cod_rec[[30]],
                        pop_dist_cod_so[[30]])

maps2 <- ggplot(combine_status|>
                  filter(year > 10), aes(x, y, fill = N_all_age)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10", na.value="white") +
  labs(fill="Abundance")+
  facet_grid(scenario ~ year, labeller = labeller(year =
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
  labs(title = "a) Cod-like") +
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks = element_blank())

maps2

figure_rec_spill <- ggarrange(maps2, maps, ncol = 1, nrow = 2, legend = "bottom", common.legend = TRUE)

figure_rec_spill

ggsave("data/figure_rec_spill.pdf", plot = figure_rec_spill, width = 12, height = 8, units = "in", dpi = 500, bg = "white")
