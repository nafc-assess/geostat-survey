
library(plotly)
library(Rstrap)
library(dplyr)
library(NAFOdown)

n_sim <- 10000

RF_spring <- strat.fun(setdet = con.setdet, lf = con.lf, program = "strat2 & strat1",
                       data.series = c("Engel","Campelen"), species = 794,
                       survey.year = c(1991:2005, 2007:2021),
                       season = "spring", NAFOdiv = c("3O"),
                       strat = c(329,332:337,339,354:356,717:722),
                       sex = c("male", "female", "unsexed"),
                       length.group = 1, group.by = "length", export = NULL, plot.results = F)


RF_fall <-strat.fun(setdet = con.setdet, lf = con.lf, program = "strat2 & strat1",
                    data.series = c("Engel","Campelen"), species = 794,
                    survey.year = c(1991:2013,2015:2021),
                    season = "fall", NAFOdiv = c("3O"),
                    strat = c(329,332:337,339,354:356,717:722),
                    sex = c("male", "female", "unsexed"),
                    length.group = 1, group.by = "length", export = NULL, plot.results = F)

spring_totals <- RF_spring$strat2$biomass$summary[, c("survey.year", "total", "var", "sample.units", "df")]
names(spring_totals) <- c("year", "spring_total", "spring_var", "spring_units", "spring_df")
fall_totals <- RF_fall$strat2$biomass$summary[, c("survey.year", "total", "var", "sample.units", "df")]
names(fall_totals) <- c("year", "fall_total", "fall_var", "fall_units", "fall_df")
totals <- merge(spring_totals, fall_totals, by = "year", all = TRUE)
totals$mean_total <- rowMeans(totals[, c("spring_total", "fall_total")], na.rm = TRUE)
totals$mean_var <- rowSums(totals[, c("spring_var", "fall_var")], na.rm = TRUE) / (2 ^ 2)
totals$mean_units <- rowMeans(totals[, c("spring_units", "fall_units")], na.rm = TRUE)
totals$mean_df <- rowMeans(totals[, c("spring_df", "fall_df")], na.rm = TRUE)

totals <- totals |>
  mutate(mean_lcl = (mean_total - abs(qt(0.025, mean_df)) * mean_units * sqrt(mean_var)),
         mean_ucl = (mean_total + abs(qt(0.025, mean_df)) * mean_units * sqrt(mean_var)),
         mean_sigma = mean_units * sqrt(mean_var),
         mean_scale = (mean_sigma ^ 2) / mean_total,
         mean_shape = mean_total / mean_scale)

bmsy_proxy <- mean(totals$mean_total)
bmsy_proxy_units <- mean(totals$mean_units)
bmsy_proxy_df <- mean(totals$mean_df)
bmsy_proxy_var <- sum(totals$mean_var) / (length(totals$mean_var) ^ 2)
bmsy_proxy_sigma <- (bmsy_proxy_units * sqrt(bmsy_proxy_var))
bmsy_proxy_scale <- (bmsy_proxy_sigma ^ 2) / bmsy_proxy # gamma dist parameter theta
bmsy_proxy_shape <- bmsy_proxy / bmsy_proxy_scale       # gamma dist parameter k
bmsy_proxy_lcl <- (bmsy_proxy - abs(qt(0.025, bmsy_proxy_df)) * bmsy_proxy_units * sqrt(bmsy_proxy_var))
bmsy_proxy_ucl <- (bmsy_proxy + abs(qt(0.025, bmsy_proxy_df)) * bmsy_proxy_units * sqrt(bmsy_proxy_var))

## Using the t-distribution does not result in realistic intervals
plot_ly(data = totals) |>
  add_markers(x = ~year, y = ~mean_total, name = "Bhat", color = "Bhat") |>
  add_segments(x = ~year, xend = ~year, y = ~mean_lcl, yend = ~mean_ucl, name = "Bhat_ci",
               color = "Bhat") |>
  add_segments(x = ~min(year), xend = ~max(year),
               y = bmsy_proxy, yend = bmsy_proxy, color = "Bmsy", name = "Bmsy") |>
  add_segments(x = ~min(year), xend = ~max(year),
               y = bmsy_proxy_lcl, yend = bmsy_proxy_lcl, name = "Bmsy_lcl",
               color = "Bmsy", linetype = I(3)) |>
  add_segments(x = ~min(year), xend = ~max(year),
               y = bmsy_proxy_ucl, yend = bmsy_proxy_ucl, name = "Bmsy_ucl",
               color = "Bmsy", linetype = I(3)) |>
  layout(xaxis = list(title = "Year"),
         yaxis = list(title = "Biomass Index"))


## Simulate distribution using gamma distribution
gamma_res <- rgamma(n_sim,
                    shape = bmsy_proxy_shape,
                    scale = bmsy_proxy_scale)

trend <- plot_ly(data = totals) |>
  add_lines(x = ~year, y = ~mean_total, color = I("steelblue")) |>
  layout(yaxis = list(title = "Biomass Index", range = c(0, max(totals$mean_total, na.rm = TRUE))))

hist <- plot_ly() |>
  add_histogram(y = gamma_res, color = I("steelblue")) |>
  add_histogram(y = gamma_res * 0.3, color = I("red"))

subplot(trend, hist, nrows = 1, widths = c(0.8, 0.2), shareY = TRUE, titleY = TRUE) |>
  hide_guides()

## Simulate biomass index for each year | mean and variance estimates to approximate CI
sim_totals <- lapply(seq.int(nrow(totals)), function(i) {
  data.frame(year = totals$year[i], sim_total = rgamma(n_sim, shape = totals$mean_shape[i],
                                                       scale = totals$mean_scale[i]))
}) |> data.table::rbindlist()

plot_ly(sim_totals, x = ~year, y = ~sim_total) |>
  add_histogram2d()

totals <- sim_totals |>
  group_by(year) |>
  summarise(gamma_median = median(sim_total),
            gamma_lcl = quantile(sim_total, probs = 0.1), # 80% intervals
            gamma_ucl = quantile(sim_total, probs = 0.9)) |>
  left_join(totals, by = "year")


totals$p_bmsy <- sapply(split(sim_totals, sim_totals$year), function(x) {
  mean((gamma_res - x$sim_total) > 0)
})
totals$p_lrp <- sapply(split(sim_totals, sim_totals$year), function(x) {
  mean((gamma_res * 0.3 - x$sim_total) > 0)
})
totals$p_fixed_lrp <- sapply(split(sim_totals, sim_totals$year), function(x) {
  mean(x$sim_total < (bmsy_proxy) * 0.3)
})

## Calculate ecdf to generate status background colours that are a graident, not discrete
bmsy_ecdf <- ecdf(gamma_res)
blim_ecdf <- ecdf(gamma_res * 0.3)

y <- seq(0, max(totals$gamma_ucl, na.rm = TRUE), length.out = 1000)
p_bmsy <- bmsy_ecdf(y)
p_blim <- blim_ecdf(y)
p_vec <- p_bmsy + p_blim

p_dat <- data.frame(year = rep(totals$year, each = length(p_vec)),
                    mean_total = y,
                    prob = rep(p_vec, length(unique(totals$year))))

p1 <- plot_ly() |>
  add_heatmap(x = ~year, y = ~mean_total, z = ~prob, data = p_dat,
              colors = rev(c("#00c853", "#ffd600", "#d50000"))) |>
  add_lines(x = ~year, y = ~mean_total, color = I("black"), data = totals) |>
  add_lines(x = ~year, y = ~gamma_lcl, linetype = I(3), color = I("black"),
            legendgroup = "gamma sim", showlegend = FALSE, data = totals) |>
  add_lines(x = ~year, y = ~gamma_ucl, linetype = I(3), color = I("black"),
            legendgroup = "gamma sim", showlegend = FALSE, data = totals) |>
  layout(yaxis = list(title = "Biomass Index", range = c(0, max(totals$gamma_ucl, na.rm = TRUE))),
         xaxis = list(title = "Year")) |>
  hide_guides()

p2 <- plot_ly(totals, x = ~year) |>
  add_lines(y = ~p_lrp, color = I("black")) |>
  layout(yaxis = list(title = "P(B < Blim)"),
         xaxis = list(title = "Year")) |>
  hide_guides()

subplot(p1, p2, nrows = 2, shareX = TRUE, titleY = TRUE)

## Stats for terminal year
delta <- gamma_res * 0.3 - sim_totals$sim_total[sim_totals$year == max(sim_totals$year)]
mean(delta > 0)

plot_ly() |>
  add_histogram(x = gamma_res * 0.3, name = "Blim") |>
  add_histogram(x = sim_totals$sim_total[sim_totals$year == 2018], name = "B") |>
  add_histogram(x = delta, name = "delta")


p1 <- ggplot() +
  geom_ribbon(aes(x = year, ymin = gamma_lcl, ymax = gamma_ucl), data = totals,
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(x = year, y = mean_total), data = totals, size = .nafo_lwd, color = "steelblue") +
  geom_hline(aes(yintercept = bmsy_proxy * 0.3), size = .nafo_lwd, color = "red") +
  geom_hline(aes(yintercept = quantile(gamma_res * 0.3, probs = 0.1)), linetype = 3, size = .nafo_lwd, color = "red") +
  geom_hline(aes(yintercept = quantile(gamma_res * 0.3, probs = 0.9)), linetype = 3, size = .nafo_lwd, color = "red") +
  theme_nafo() +
  ylab("Combined Biomass Index") +
  scale_y_continuous(labels = scales::label_number(suffix = "", scale = 1e-7)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

p2 <- ggplot() +
  geom_bar(aes(x = year, y = p_lrp, fill = p_lrp), data = totals, stat = "identity") +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  guides(fill = "none") +
  xlab("Year") + ylab("P(B < Blim)") +
  theme_nafo()

red_plot <- cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(0.6, 0.3), align = "v")

save(red_plot, totals, bmsy_proxy, gamma_res, file = "Gamma_SCR/data/red_plot.rda")

write.csv(totals, file = "Gamma_SCR/data/red_vals.csv", row.names = FALSE)
