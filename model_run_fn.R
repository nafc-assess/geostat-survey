#' Model functions for parallel runs
#' Model formulas:
formula1 = count ~ 0 + as.factor(year)
formula2 = count ~ 0 + as.factor(year) + poly(log(depth), 2)

#' Negative-binomial spatiotemporal model
#' @param data: sdmTMB dataset created
#' @param mesh: sdmTMB mesh created
#' @param formula: model formula defined above
#' @param range_gt: prior (km)
#' @param sigma_lt: prior (km)
#' @param type: modelling approach classification based on formula (NB or NB + depth)
#' @param scenario: survey change scenario
#' @param species: species name
#' @param newdata: prediction dataset created
#
#' @return: population index value (N) in each year, convergence report, AIC value for each population and scenarios
model_run_NB <- function(data, mesh, formula, range_gt, sigma_lt, type, scenario, species, newdata){

  fit <- try({sdmTMB(formula,
                     data = data,
                     mesh = mesh,
                     offset = data$offset,
                     time = "year",
                     family = nbinom2(link = "log"),
                     spatial = TRUE,
                     spatiotemporal = "IID",
                     priors = sdmTMBpriors(
                       matern_s = pc_matern(range_gt = range_gt, sigma_lt = sigma_lt),
                       matern_st = pc_matern(range_gt = range_gt, sigma_lt = sigma_lt)),
                     share_range = FALSE,
                     control = sdmTMBcontrol(newton_loops = 1L)
  )})

  if(isFALSE(inherits(fit, "sdmTMB"))){
    fit_fail <- data.frame(type = type, scenario = scenario, species = species,
                           sim = unique(data$sim), pop = unique(data$pop), status = "fail")
    return(fit_fail)}

  else{

    if (max(fit$gradients) > 0.001) {
      fit <- try({
        run_extra_optimization(fit, newton_loops = 1L, nlminb_loops = 1L)
      })
    }
    convergence_report <- tryCatch(expr = {sanity(fit)},
                                   error = function(e){
                                     return(convergence_report = list("hessian_ok"= NA,
                                                                      "eigen_values_ok" = NA,
                                                                      "nlminb_ok"= NA,
                                                                      "range_ok"= NA,
                                                                      "gradients_ok"= NA,
                                                                      "se_magnitude_ok"= NA,
                                                                      "se_na_ok"= NA,
                                                                      "sigmas_ok"= NA,
                                                                      "all_ok"= NA))
                                   }
    )

    pred <- predict(fit,
                    newdata = newdata,
                    return_tmb_object = TRUE)

    index <- get_index(pred, area = pred$data$area, bias_correct = TRUE) |>
      cbind(convergence_report) |>
      cbind(AIC = AIC(fit)) |>
      mutate(type = type,
             N = est,
             scenario = scenario,
             species = species,
             sim = unique(data$sim),
             pop = unique(data$pop),
             status = "pass")
    return(index)
  }
}

#' Tweedie spatiotemporal model
#' @param data: sdmTMB dataset created
#' @param mesh: sdmTMB mesh created
#' @param formula: model formula defined above
#' @param range_gt: prior (km)
#' @param sigma_lt: prior (km)
#' @param type: modelling approach classification based on formula (TW or TW + depth)
#' @param scenario: survey change scenario
#' @param species: species name
#' @param newdata: prediction dataset created
#
#' @return: population index value (N) in each year, convergence report, AIC value for each population and scenarios
model_run_TW <- function(data, mesh, formula, range_gt, sigma_lt, type, scenario, species, newdata){
  matern = pc_matern(range_gt = range_gt, sigma_lt = sigma_lt)

  fit <- try({sdmTMB(formula,
                     data = data,
                     mesh = mesh,
                     offset = data$offset,
                     time = "year",
                     family = tweedie(),
                     spatial = TRUE,
                     spatiotemporal = "IID",
                     priors = sdmTMBpriors(
                       matern_s = matern,
                       matern_st = matern),
                     share_range = FALSE,
                     control = sdmTMBcontrol(newton_loops = 1L)
  )})

  if(isFALSE(inherits(fit, "sdmTMB"))){
    fit_fail <- data.frame(type = type, scenario = scenario, species = species,
                           sim = unique(data$sim), pop = unique(data$pop), status = "fail")
    return(fit_fail)}

  else{

    if (max(fit$gradients) > 0.001) {
      fit <- try({
        run_extra_optimization(fit, newton_loops = 1L, nlminb_loops = 1L)
      })
    }
    convergence_report <- tryCatch(expr = {sanity(fit)},
                                   error = function(e){
                                     return(convergence_report = list("hessian_ok"= NA,
                                                                      "eigen_values_ok" = NA,
                                                                      "nlminb_ok"= NA,
                                                                      "range_ok"= NA,
                                                                      "gradients_ok"= NA,
                                                                      "se_magnitude_ok"= NA,
                                                                      "se_na_ok"= NA,
                                                                      "sigmas_ok"= NA,
                                                                      "all_ok"= NA))
                                   }
    )

    pred <- predict(fit,
                    newdata = newdata,
                    return_tmb_object = TRUE)

    index <- get_index(pred, area = pred$data$area, bias_correct = TRUE) |>
      cbind(convergence_report) |>
      cbind(AIC = AIC(fit)) |>
      mutate(type = type,
             N = est,
             scenario = scenario,
             species = species,
             sim = unique(data$sim),
             pop = unique(data$pop),
             status = "pass")
    return(index)
  }
}

#' Delta-Gamma spatiotemporal model
#' @param data: sdmTMB dataset created
#' @param mesh: sdmTMB mesh created
#' @param formula: model formula defined above
#' @param range_gt: prior (km)
#' @param sigma_lt: prior (km)
#' @param type: modelling approach classification based on formula (DG or DG + depth)
#' @param scenario: survey change scenario
#' @param species: species name
#' @param newdata: prediction dataset created
#
#' @return: population index value (N) in each year, convergence report, AIC value for each population and scenarios
model_run_DG <- function(data, mesh, formula, range_gt, sigma_lt, type, scenario, species, newdata){

  fit <- try({sdmTMB(formula,
                     data = data,
                     mesh = mesh,
                     offset = data$offset,
                     time = "year",
                     family = delta_gamma(),
                     spatial = TRUE,
                     spatiotemporal = list("off", "IID"),
                     priors = sdmTMBpriors(
                       matern_s = pc_matern(range_gt = range_gt, sigma_lt = sigma_lt),
                       matern_st = pc_matern(range_gt = range_gt, sigma_lt = sigma_lt)),
                     share_range = FALSE,
                     control = sdmTMBcontrol(newton_loops = 1L)
  )})

  if(isFALSE(inherits(fit, "sdmTMB"))){
    fit_fail <- data.frame(type = type, scenario = scenario, species = species,
                           sim = unique(data$sim), pop = unique(data$pop), status = "fail")
    return(fit_fail)}

  else{

    if (max(fit$gradients) > 0.001) {
      fit <- try({
        run_extra_optimization(fit, newton_loops = 1L, nlminb_loops = 1L)
      })
    }
    convergence_report <- tryCatch(expr = {sanity(fit)},
                                   error = function(e){
                                     return(convergence_report = list("hessian_ok"= NA,
                                                                      "eigen_values_ok" = NA,
                                                                      "nlminb_ok"= NA,
                                                                      "range_ok"= NA,
                                                                      "gradients_ok"= NA,
                                                                      "se_magnitude_ok"= NA,
                                                                      "se_na_ok"= NA,
                                                                      "sigmas_ok"= NA,
                                                                      "all_ok"= NA))
                                   }
    )

    pred <- predict(fit,
                    newdata = newdata,
                    return_tmb_object = TRUE)

    index <- get_index(pred, area = pred$data$area, bias_correct = TRUE) |>
      cbind(convergence_report) |>
      cbind(AIC = AIC(fit)) |>
      mutate(type = type,
             N = est,
             scenario = scenario,
             species = species,
             sim = unique(data$sim),
             pop = unique(data$pop),
             status = "pass")
    return(index)
  }
}
