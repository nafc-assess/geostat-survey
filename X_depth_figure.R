
library(ggplot2)
library(patchwork)
library(SimSurvey)

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


xyz <- data.frame(grid) # works with new stars based raster object (latest GitHub version)
(profile <- ggplot(data = xyz) +
  geom_ribbon(aes(x = x, ymin = max(depth), ymax = depth), fill = "grey60", color = "grey50") +
  scale_y_reverse(expand = c(0, 0), limits = c(max(xyz$depth), 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  ylab("Depth (m)") +
  theme_bw())

depth <- sort(unique(xyz$depth))
yt_effect <- yt_depth(data.frame(age = 1, depth = depth))
cod_effect <- cod_depth(data.frame(age = 1, depth = depth))
depth_effect <- data.frame(depth = c(depth, depth),
                           species = c(rep("cod-like", length(depth)),
                                       rep("yellowtail-like", length(depth))),
                           effect = c(cod_effect, yt_effect))

(parabola <- ggplot(data = depth_effect) +
  geom_path(aes(x = exp(effect), y = depth, color = species), size = 1) +
  scale_y_reverse(expand = c(0, 0), limits = c(max(xyz$depth), 0)) +
  scale_color_brewer(palette = "Set1", name = "Species") +
  xlab("Depth effect") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()))

profile + parabola +
  plot_layout(widths = c(2, 1))

