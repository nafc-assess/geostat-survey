library(data.table)

# strat_means <- function(data = NULL, metric = NULL, strat_groups = NULL,
#                         survey_groups = NULL, confidence = 95) {
#
#   Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <-
#     meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <-
#     sumYst_ucl <- NULL
#
#   ## set-up

  lc <- (100 - 95) / 200
  uc <- (100 - 95) / 200 + (confidence / 100)

  d <- data.table::copy(design_index_sc2_20[[5]]$setdet) # make a copy of the provided data.table

  metric="n"
  strat_groups= c("year", "strat", "strat_area")
  survey_groups= "year"


  d <- d[,  c(strat_groups,  metric),  with = FALSE]

  setnames(d, names(d), c(strat_groups, "metric"))

  setkeyv(d,  strat_groups)

  ## strat.tab includes strat-level means,  variances,  etc.

  strat_tab <- d[, list(sumYh = sum(metric), meanYh = mean(metric),
                        varYh = stats::var(metric), nh = .N), by = strat_groups]

  strat_tab[, Nh := strat_area / 0.03] # number of sample units in each strat

  strat_tab[, Wh := Nh / sum(Nh), by = survey_groups]
  strat_tab[, total := Nh * sumYh / nh]
  strat_tab[, gh := Nh * (Nh - nh) / nh]

  ## survey.tab includes large-scale means,  such as mean weight per tow,  and totals,  such as total biomass + associated confidence intervals

  survey_tab <- strat_tab[, list(n = sum(nh), N = sum(Nh), meanYst = sum(Wh * meanYh),
                                 varYst = (1 / ((sum(Nh)) ^ 2)) * sum(gh * varYh),
                                 df = ((sum(gh * varYh)) ^ 2) / (sum((gh ^ 2 * varYh ^ 2) / (nh - 1)))),
                          by = survey_groups]
  survey_tab[, meanYst_lcl := (meanYst - (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, meanYst_ucl := (meanYst + (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
  survey_tab[, sumYst := N * meanYst]
  survey_tab[, sumYst_lcl := (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[, sumYst_ucl := (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst))]
  survey_tab[sapply(survey_tab,  is.nan)] <- NA

  ## Rename cols
  survey_tab <- survey_tab[, c(survey_groups,  "n", "N", "df", "varYst",
                               "meanYst",  "meanYst_lcl", "meanYst_ucl",
                               "sumYst", "sumYst_lcl", "sumYst_ucl"),  with = FALSE]
  survey_tab$varYst <- sqrt(survey_tab$varYst) # convert to sd
  setnames(survey_tab, names(survey_tab), c(survey_groups, "sets", "sampling_units", "df", "sd",
                                            "mean", "mean_lcl", "mean_ucl",
                                            "total", "total_lcl", "total_ucl"))

  # ## return results (here I only return the survey_tab to minimize details and object size)
  # survey_tab
  #}

strat_tab[is.na(strat_tab$varYh),]

design_index_sc2_20[[5]]$setdet[design_index_sc2_20[[5]]$setdet$year==11 & design_index_sc2_20[[5]]$setdet$strat==9,]

design_index_sc2_20[[5]]$setdet[design_index_sc2_20[[5]]$setdet$year==13 & design_index_sc2_20[[5]]$setdet$strat==12,]

survey_tab

a <- strat_tab[strat_tab$year==11,]

b<- na.omit(a)

#varYst ("sd")

(1 / ((sum(a$Nh)) ^ 2)) * sum(a$gh * a$varYh)

(1 / ((sum(b$Nh)) ^ 2)) * sum(b$gh * b$varYh)

survey_tab_b <- b[, list(n = sum(nh), N = sum(Nh), meanYst = sum(Wh * meanYh),
                               varYst = (1 / ((sum(Nh)) ^ 2)) * sum(gh * varYh),
                               df = ((sum(gh * varYh)) ^ 2) / (sum((gh ^ 2 * varYh ^ 2) / (nh - 1)))),
                        by = survey_groups]
survey_tab_b[, meanYst_lcl := (meanYst - (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
survey_tab_b[, meanYst_ucl := (meanYst + (sqrt(varYst)) * abs(stats::qt(lc,  df)))]
survey_tab_b[, sumYst := N * meanYst]
survey_tab_b[, sumYst_lcl := (sumYst - abs(stats::qt(lc, df)) * N * sqrt(varYst))]
survey_tab_b[, sumYst_ucl := (sumYst + abs(stats::qt(lc, df)) * N * sqrt(varYst))]
setnames(survey_tab_b, names(survey_tab_b), c(survey_groups, "sets", "sampling_units", "df", "sd",
                                          "mean", "mean_lcl", "mean_ucl",
                                          "total", "total_lcl", "total_ucl"))
