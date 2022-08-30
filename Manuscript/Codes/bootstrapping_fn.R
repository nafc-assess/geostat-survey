sumYst <- function(data, i = seq_len(nrow(data))) {
  data[i,] |>
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
  b <- boot::boot(data, statistic = sumYst, strata = data$strat, R = reps, parallel = "multicore", ncpus=6)
  suppressWarnings(bci <- boot::boot.ci(b, type = "bca"))
  mean_boot <- mean(b$t)
  boot <- data.table(
    sim = mean(data$sim),
    N = sumYst(data),
    pop = mean(data$pop),
    mean_boot= mean_boot,
    lwr = bci$bca[[4]],
    upr = bci$bca[[5]],
    cv = sd(b$t) / mean_boot,
    type = "Bootstrapped")
  return(boot)
  gc()
}

boot_wrapper <- function(data, reps) {
  out <- data |>
    split(data$year) |>
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  return(out)
}



