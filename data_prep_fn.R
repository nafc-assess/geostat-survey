sdm_data_fn <- function(x) {
  dat <- as_tibble(x) |>
    dplyr::select(x, y, set, sim, pop, year, count = n, tow_area, area = cell_area, depth) |>
    mutate(offset = log(tow_area), density = count / tow_area)
}

mesh_sdm_fn <- function(sdm_data, existing_mesh = NULL){
  if (is.null(existing_mesh)) {
    mesh <- sdmTMB::make_mesh(sdm_data,
      xy_cols = c("x", "y"),
      cutoff = 45)
  } else {
    mesh <- sdmTMB::make_mesh(sdm_data,
      xy_cols = c("x", "y"),
      mesh = existing_mesh)
  }
  mesh
}

sdm_newdata_fn <- function(survey, sdm_data){
  newdata <- as_tibble(dplyr::select(survey$grid_xy, x, y, depth)) |> distinct()
  newdata <- purrr::map_dfr(sort(unique(sdm_data$year)), ~ bind_cols(newdata, year = .)) |>
    mutate(offset = 0, area = sdm_data$area[1])
}

make_standard_mesh <- function() {
  loc.bnd <- matrix(c(-150, -150, 150, -150, 150, 150, -150, 150), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd)
  inla_mesh <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(50, 80),
    offset = c(0, 25),
    cutoff = 10
  )
  # plot(inla_mesh)
  # axis(1);axis(2)
  # points(temp$x, temp$y, col = "red")
  # inla_mesh$n
  inla_mesh
}
