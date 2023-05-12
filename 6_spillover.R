### Spillover effect

#####

##  Survey scenarios

library(data.table)
library(raster)
library(sp)
library(sf)
library(SimSurvey)
library(tidyverse)
set.seed(1)

# Making a grid with SimSurvey
grid <- make_grid(x_range = c(-150, 150),
                  y_range = c(-150, 150),
                  res = c(10, 10),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1000),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 40),
                  strat_splits = 2,
                  method = "spline")

grid_sf <- st_as_sf(grid)

# Creating a buffer area

grid_sf_buffered <- function(x, fraction) {
  width_half <- sqrt(st_area(st_union(x)) * fraction)/2
  buffer <- st_buffer(st_geometry(st_union(x)), -width_half)
  st_geometry(x)[!st_is_empty(buffer)] = buffer[!st_is_empty(buffer)]
  return(st_union(x))
}
grid_sf_buffered_perc30 <- grid_sf_buffered(grid_sf, 0.3)

### selecting a random point and creating the blocked area

set.seed(30)

block_poly_sf_fn <- function(x, fraction, buffer) {
  grid_centroid <- st_centroid(x)
  pointpool <- st_intersection(grid_centroid, buffer)

  selected <-  pointpool |>
    #group_by()|>
    slice_sample(n = 1) |>
    st_coordinates(selected)|>
    data.table()

  width_half <- sqrt(as.numeric(st_area(st_union(x))) * fraction)/2

  yplus <- selected$Y + width_half
  xplus <- selected$X + width_half
  yminus <- selected$Y - width_half
  xminus <- selected$X - width_half

  x1 <- c(xminus,xplus,xplus,xminus,xminus) # a(SW), b(SE), c(NE), d(NW), a(SW) - closing the vertices
  y1 <- c(yplus,yplus,yminus,yminus,yplus)

  block <- sp::Polygon(cbind(x1,y1))
  block <- sp::Polygons(list(block), ID = "1")
  block_poly <- sp::SpatialPolygons(list(block))
  block_poly_sf <- st_as_sf(block_poly)

  return(block_poly_sf)
}

block_poly_sf_30 <- block_poly_sf_fn(grid_sf, 0.3, grid_sf_buffered_perc30)

st_crs(block_poly_sf_30) = st_crs(grid_sf)

mpa <- block_poly_sf_30


# ###### multrings
# st_multiringbuffer <- function(x,n,d){
#   buf<-function(dist){st_buffer(x,dist) %>% mutate(dist=dist)}
#   out<-purrr::map_df(seq(1:n)*d, buf) |> mutate(sp = 1/row_number())
# }
#
# layer <- st_multiringbuffer(x=mpa, n=10, d=5)
# layer
#
# mapview(layer) + mapview(mpa)
#
# layer_buff_only <- st_difference(layer, mpa)
# mapview(layer_buff_only)
# mpa$sp <- 1L

# mpa_sp <- st_join(mpa, layer_buff_only, by_feature = FALSE)
# mpa_sp
# mapview(mpa_sp)
##############

##########

mpa$mpa <- 1L
mpa

cell_mpa <- st_intersection(mpa, grid_sf)
mapview(cell_mpa)

grid_sf$mpa <- 0
cell_wo_mpa <- st_difference(grid_sf, mpa)
mapview(cell_wo_mpa)

output <- cell_mpa %>%
  bind_rows(cell_wo_mpa)

output$mpa.1 <- NULL
mapview(output)

#save(output, file = "./Data/output.Rdata")

####

buff1 <- st_buffer(mpa, 10)
mapview(buff1)
buff1_mpa <- st_difference(buff1, mpa)
mapview(buff1_mpa)
buff1_mpa_cell <- st_intersection(buff1_mpa, grid_sf)
mapview(buff1_mpa_cell)
buff1_mpa_cell$mpa <- 0.5
buff1_mpa_cell$mpa.1 <- NULL
buff1_mpa_cell$mpa.2 <- NULL
mapview(buff1_mpa_cell)
buff1_mpa_cell

buff_mpa_cell <- cell_mpa %>%
  bind_rows(buff1_mpa_cell)
mapview(buff_mpa_cell)

cell_wo_mpa_buff <- st_difference(grid_sf, buff1)

output2 <- buff_mpa_cell %>%
  bind_rows(cell_wo_mpa_buff)

mapview(output2)

## Sim non linear function -----------------------------------------------------
## Idea is to increase the probability the fish will occur in the MPA

grid <- make_grid(x_range = c(-150, 150),
                  y_range = c(-150, 150),
                  res = c(10, 10),
                  shelf_depth = 200,
                  shelf_width = 100,
                  depth_range = c(0, 1000),
                  n_div = 1,
                  strat_breaks = seq(0, 1000, by = 40),
                  strat_splits = 2,
                  method = "spline")

## Hack: add mpa column to grid object
grid_xy <- as.data.frame(grid)
grid_xy$mpa <- 0
grid_xy$mpa[grid_xy$y < 50 & grid_xy$y > -120 & grid_xy$x > -65 & grid_xy$x < 100] <- 1
grid$mpa <- grid_xy$mpa

plot(x ~ y, data = grid_xy, col = grid_xy$mpa)

sim <- list(ages = 1:20, years = 1:20)
grid_dat <- data.table::data.table(as.data.frame(grid))
grid_edat <- grid_dat
i <- rep(seq(nrow(grid_edat)), times = length(sim$ages))
a <- rep(sim$ages, each = nrow(grid_edat))
grid_edat <- grid_edat[i]
grid_edat$age <- a
i <- rep(seq(nrow(grid_edat)), times = length(sim$years))
y <- rep(sim$years, each = nrow(grid_edat))
grid_edat <- grid_edat[i]
grid_edat$year <- y
grid_edat <- grid_edat[order(grid_edat$cell, grid_edat$year, grid_edat$age), ] # sort to align with error array


depth_par <- sim_nlf(formula = ~ alpha + ifelse(year > 10, (beta * mpa * (year - 10) / 10), 0) - ((depth - mu) ^ 2) / (2 * sigma ^ 2),
                     coeff = list(alpha = 0, mu = 200, sigma = 70, beta = 0.2))

depth <- depth_par(grid_edat)

grid_edat$depth_effect <- depth

library(plotly)
grid_edat %>%
  plot_ly(x = ~depth, y = ~exp(depth_effect), split = ~mpa,
          frame = ~year) %>%
  add_lines()



