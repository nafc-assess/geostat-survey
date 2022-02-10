#####################################
### Creating Survey Data ###


library(SimSurvey)
library(tidyverse)

set.seed(17)
pop <- sim_abundance(ages = 1:20, #fish age vector
                     years = 1:10, #year vector
                     Z = sim_Z(log_mean = log(0.6), #Total mortality matrix -- rnorm
                               log_sd = 0.1,
                               phi_age = 0.9,
                               phi_year = 0.6,
                               plot = FALSE),
                     R = sim_R(log_mean = log(3e+07), #Recruitment vector -- rnorm
                               log_sd = 0.5,
                               random_walk = TRUE,
                               plot = FALSE),
                     N0 = sim_N0(N0 = "exp", #Starting abundance vector -- exponential Z; N0 at age 1 is R0
                                 plot = FALSE),
                     growth = sim_vonB(Linf = 100, #abundance-at-age matrix
                                       L0 = 5,
                                       K = 0.2,
                                       log_sd = 0.05,
                                       plot = FALSE,
                                       length_group = 3))  %>%
  sim_distribution(grid = make_grid(x_range = c(-150, 150), # km
                                    y_range = c(-150, 150), # km
                                    res = c(10, 10), # km
                                    shelf_depth = 200, #Approximate depth of the shelf in m
                                    shelf_width = 50, #Approximate width of the shelf in km
                                    depth_range = c(10, 1000), #m
                                    n_div = 1, #number of division
                                    strat_breaks = seq(0, 1000, by = 40), #strata depth breaks, meter
                                    strat_splits = 2, # horizontally split strat
                                    method = "spline"), #spline, loess, linear interpolation options
                   ays_covar = sim_ays_covar(range = 300, #adding noise correlated across space, year, and age dimensions -- covariance matrix, higher value higher spatial correlation
                                             phi_age = 0.8, #strong correlation across ages
                                             phi_year = 0.6),# correlation across years
                   depth_par = sim_parabola(mu = 200, #defining relationship between abundance and depth, i.e., adding a parabolic depth 'preference', fish prefers occurring around 200 meters
                                            sigma = 70, plot=FALSE))

survey <- sim_survey(
  pop, #sim_distribution
  n_sims = 1, #number of survey
  q = sim_logistic(k = 2, x0 = 3, plot = FALSE), #function of catchability at age (k= The steepness of the curve, x0= x-value of the sigmoid's midpoint, curve max value=1)
  #trawl_dim = c(1.5, 0.02), #Trawl width and distance (default=same units as grid)
  #resample_cells = FALSE, #Allow resampling of sampling units?
  binom_error = TRUE, #Impose binomial error?
  #min_sets = 2, #Minimum number of sets per strat
  set_den = 2/1000, #set per strata ~ strata area with a set density rule, sample size = sum(set numbers per strat) * years
  #lengths_cap = 500, #Maximum number of lengths measured per set
  #ages_cap = 10,
  #age_sampling = "stratified",
  #age_length_group = 1,
  #age_space_group = "division",
  #light = TRUE
) %>%
  SimSurvey::run_strat()%>%
  SimSurvey::strat_error()


#### Index calculation to compare later

strat_abund <- tibble::as_tibble(survey$total_strat) %>%
  mutate(N = total, upr =total_ucl, lwl = total_lcl, type = "Design-based") %>%
  select(year, N, type, lwl, upr)

#####################################
### Stratified bootstrap analysis ###


#survey$setdet
survey_year1 <- subset(survey$setdet, year==1)

# Samples are in stratum-level n (number of fish) meanYh = mean(n)
# To calculate to total number of fish in the area, we need sumYst := N * meanYst
# meanYst = sum(Wh * meanYh)
# Wh = strat_area/(1.5 * 0.02)
# N = sum(Nh)

#### The calculation
survey_year1%>%
  ### stratum level
  group_by(year, strat, strat_area) %>%
  summarise(meanYh = mean(n), sumYh = sum(n), nh=n(), .groups = "drop_last") %>%
  mutate(Nh = strat_area/(1.5 * 0.02)) %>%
  group_by(year) %>% #each year grouping for sum(Nh)
  mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)%>%
  ### year level
  summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") %>% ## N is actually same for each strata (total grid area), there should be another way to call N...
  pull(sumYst)

#### Statistic for boot
sumYst <- function(dat, i = seq_len(nrow(dat))) {
  dat[i, ] %>%
    ### stratum level
    group_by(year, strat, strat_area) %>%
    summarise(meanYh = mean(n), sumYh = sum(n), nh=n(), .groups = "drop_last") %>%
    mutate(Nh = strat_area/(1.5 * 0.02)) %>%
    group_by(year) %>%
    mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)%>%
    ### year level
    summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") %>%
    pull(sumYst)
}

#### Testing

sumYst(survey$setdet) ### It's working

library(boot)

b <- boot(survey$setdet, statistic = sumYst, strata = survey$setdet$strat, R = 2000) #t0(originals) are matching with Simsurvey index.
#b$t, index is year level, index=1 means first year
#R= 1000 didn't work, "with a small number of replications, empinf sometimes fails and returns a vector of NA values" issue

bci_1 <- boot.ci(b, index = 1, conf = 0.95, type = "bca") # index=1 means first year

#### Boot function

boot_one_year <- function(x, reps) {
  b <- boot::boot(x, statistic = sumYst, strata = x$strat, R = reps)
  suppressWarnings(bci <- boot::boot.ci(b, type = "bca"))
  tibble::tibble(
    index = sumYst(x),
    mean_boot = mean(b$t),
    median_boot = median(b$t),
    lwr = bci$bca[[4]],
    upr = bci$bca[[5]],
    cv = sd(b$t) / mean(b$t),
  )
}

#### Boot wrapper for each year

boot_wrapper <- function(dat, reps) {
  out <- dat %>%
    split(dat$year) %>%
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  out
}


#### Stratified bootstrapping


set.seed(15)
boot_abund <- boot_wrapper(survey$setdet, 2000)

boot_abund
strat_abund

### Let's make a graph
boot_abund_i <- boot_abund %>%
  mutate(N = mean_boot, upr =upr, lwl = lwr, type = "Bootstrap") %>%
  select(year, N, type, lwl, upr)

index <- rbind(boot_abund_i, strat_abund)

library(ggplot2)
library(plotly)

index %>%
  plot_ly(x = ~year, color = ~type, legendgroup = ~type) %>%
  add_ribbons(ymin = ~lwl, ymax = ~upr, line = list(width = 0), showlegend = FALSE) %>%
  add_lines(y = ~N)






