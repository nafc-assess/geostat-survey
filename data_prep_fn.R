sdm_data_fn <- function(x) {
  dat <- as_tibble(x) |>
    dplyr::select(x, y, set, sim, pop, year, count = n, tow_area, area = cell_area, depth) |>
    mutate(offset = log(tow_area), density = count / tow_area)
}

mesh_sdm_fn <- function(sdm_data){
  mesh <- sdmTMB::make_mesh(sdm_data,
                            xy_cols = c("x", "y"),
                            cutoff = 45)
}

sdm_newdata_fn <- function(survey, sdm_data){
  newdata <- as_tibble(dplyr::select(survey$grid_xy, x, y, depth)) |> distinct()
  newdata <- purrr::map_dfr(sort(unique(sdm_data$year)), ~ bind_cols(newdata, year = .)) |>
    mutate(offset = 0, area = sdm_data$area[1])
}
