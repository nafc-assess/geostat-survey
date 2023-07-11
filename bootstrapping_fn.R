#' Bootstrapping loop functions to calculate CI of design-based estimates

# sumYst <- function(data, i = seq_len(nrow(data))) {
#   data[i,] |>
#     ### stratum level
#     group_by(year, strat, strat_area) |>
#     summarise(meanYh = mean(n), tow_area = mean(tow_area), .groups = "drop_last") |>
#     mutate(Nh = strat_area/(tow_area)) |>
#     group_by(year) |>
#     mutate(N = sum(Nh), Wh = Nh/N, WhmeanYh = Wh * meanYh)|>
#     ### year level
#     summarise(sumYst= mean(N) * sum(WhmeanYh), .groups = "drop_last") |>
#     pull(sumYst)
# }

sumYst <- function(data, i = seq_len(nrow(data))) {
  dt <- data.table(data)
  dt <- dt[i,]

  # stratum level
  dt <- dt[, .(meanYh = mean(n), tow_area = mean(tow_area)), by = .(year, strat, strat_area)][
    , Nh := strat_area / tow_area][
      , .(N = sum(Nh), Wh = Nh / sum(Nh), WhmeanYh = Nh / sum(Nh) * meanYh), by = year]

  # year level
  sumYst <- dt[, .(sumYst = mean(N) * sum(WhmeanYh)), by = year][
    , sumYst]

  return(sumYst)
}

boot_one_year <- function(data, reps) {
  b <- boot::boot(data, statistic = sumYst, strata = data$strat, R = reps)
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
}

boot_wrapper <- function(data, reps) {
  out <- data |>
    split(data$year) |>
    purrr::map_dfr(boot_one_year, reps = reps, .id = "year")
  out$year <- as.numeric(out$year)
  return(out)
}




