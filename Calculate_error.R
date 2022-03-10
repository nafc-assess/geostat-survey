hessian <- Sdreport$pdHess
gradient <- rich$opt$max_gradient <= 0.001
hessian == TRUE & gradient == TRUE


error_stats <- function(error) {
  c(ME = mean(error),
    MAE = mean(abs(error)),
    MSE = mean(error ^ 2),
    RMSE = sqrt(mean(error ^ 2)))
}

############# Calculating errors

result_base$error <- result_base$N - result_base$true




strat_error <- function(sim) {

  total <- age <- NULL

  ## total_strat
  I_hat <- sim$total_strat[, list(sim, year, total)]
  names(I_hat) <- c("sim", "year", "I_hat")
  I <- data.frame(year = sim$years, I = colSums(sim$I))
  comp <- merge(I_hat, I, by = "year")
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$total_strat_error <- comp
  sim$total_strat_error_stats <- means

  ## length_strat
  I_hat <- sim$length_strat[, list(sim, year, length, total)]
  names(I_hat) <- c("sim", "year", "length", "I_hat")
  sly <- expand.grid(sim = seq(max(sim$total_strat$sim)),
                     year = sim$years, length = sim$lengths)
  I_hat <- merge(sly, I_hat, by = c("sim", "year", "length"), all = TRUE) # expand to all lengths
  I_hat$I_hat[is.na(I_hat$I_hat)] <- 0                                    # fill missing lengths with zero
  I <- as.data.frame.table(sim$I_at_length, responseName = "I")
  I$year <- as.numeric(as.character(I$year))
  I$length <- as.numeric(as.character(I$length))
  comp <- merge(data.table(I_hat), data.table(I), by = c("year", "length"))
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$length_strat_error <- comp
  sim$length_strat_error_stats <- means

  ## age_strat
  I_hat <- sim$age_strat[, list(sim, year, age, total)]
  names(I_hat) <- c("sim", "year", "age", "I_hat")
  say <- expand.grid(sim = seq(max(sim$total_strat$sim)),
                     year = sim$years, age = sim$ages)
  I_hat <- merge(say, I_hat, by = c("sim", "year", "age"), all = TRUE) # expand to all ages
  I_hat$I_hat[is.na(I_hat$I_hat)] <- 0                                 # fill missing ages with zero
  I <- as.data.frame.table(sim$I, responseName = "I")
  I$year <- as.numeric(as.character(I$year))
  I$age <- as.numeric(as.character(I$age))
  comp <- merge(data.table(I_hat), data.table(I), by = c("year", "age"))
  comp$error <- comp$I_hat - comp$I
  means <- error_stats(comp$error)
  sim$age_strat_error <- comp
  sim$age_strat_error_stats <- means

  sim

}
