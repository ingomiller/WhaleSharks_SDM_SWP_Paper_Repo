

#_____________ use blockCV package 




library(blockCV)
library(tmap)
library(tidyverse)
library(terra)
library(sf)
library(dynamicSDM)
library(patchwork)
source("R/00_Helper_Functions.R") 





# tracks <- readRDS( "MODEL_DATA.rds") 


tracks_m <- tracks |>
  sf::st_make_valid() |>
  sf::st_zm(drop = TRUE, what = "ZM") |>
  sf::st_transform("EPSG:3577")


input_occ <- tracks_m |> 
  dplyr::transmute(id, 
                   date, 
                   PA, 
                   depth,
                   slope,
                   roughness, 
                   dist2000, 
                   thetao, 
                   uv, 
                   mltost, 
                   chl, 
                   wz
  )




#load monthly mean raster stack of continous variables 
input_raster_stack <- terra::rast("PREDICTORS.tif")


plot(input_raster_stack)


# crop to extent of Occurrence data 
v <- terra::vect(tracks_m)
v
bb <- terra::ext(v)
buff <- 10000 

bb <- sf::st_bbox(tracks)


bb_exp <- c(
  xmin = base::max(-180, base::as.numeric(bb["xmin"]) - 0.5),
  ymin = base::max( -90, base::as.numeric(bb["ymin"]) - 0.5),
  xmax = base::min( 180, base::as.numeric(bb["xmax"]) + 0.5),
  ymax = base::min(  90, base::as.numeric(bb["ymax"]) + 0.5)
)

bbox_poly_wgs <- sf::st_as_sfc(sf::st_bbox(bb_exp, crs = sf::st_crs(4326)))
bbox_vect_wgs <- terra::vect(bbox_poly_wgs)
bbox_vect_r   <- bbox_vect_wgs |> terra::project(terra::crs(input_raster_stack))

raster_crop <- input_raster_stack |> terra::crop(bbox_vect_r, snap = "out")
terra::crs(raster_crop) <- "EPSG:4326"

# reporject to metric 
raster_m <- raster_crop |>
  terra::project("EPSG:3577", method = "bilinear")


terra::plot(
  raster_m[[1]]
)
terra::points(v[v$PA == 0, ], pch = 21, cex = 0.5, col = "black", bg = "yellow")
terra::points(v[v$PA == 1, ], pch = 21, cex = 0.5, col = "black", bg = "red")


input_raster <- raster_m


tmap::tm_shape(input_raster[[names(input_raster) != "id"]]) +
  tm_raster(
    col.scale = tm_scale_continuous(values = gray.colors(10)),
    col.legend = tm_legend_hide()
  ) +
  tmap::tm_shape(input_occ) +
  tmap::tm_dots(
    fill = "PA",
    fill.scale = tm_scale_categorical(),
    size = 0.5,
    fill_alpha = 0.5
  )



# test spatial autocreelation distance
set.seed(44)


sac <- blockCV::cv_spatial_autocor(r = input_raster,
                                    x = input_occ, 
                                    column =  "PA", 
                                    num_sample = 10000,
                                    plot = TRUE)
sac2$range





## Spatial Blocks

range = 500000

sb <- blockCV::cv_spatial(x = input_occ,
                           column = "PA", 
                           r = input_raster,
                           k = 5, 
                           size = range, 
                           hexagon = TRUE,
                           selection = "random", 
                           iteration = 100, 
                           progress = TRUE,
                           seed = 666,
                           biomod2 = FALSE) 


sb$records

blockCV::cv_plot(cv = sb, 
                 x = input_occ,
                 r = input_raster,
                 num_plots = 1:5,
                 nrow = 2, 
                 points_alpha = 0.5)



sb_sim <- blockCV::cv_similarity(cv = sb,
                                  x = input_occ,
                                  r = input_raster,
                                  num_plot = 1:5,
                                  method = "MESS",
                                  num_sample = 10000,
                                  jitter_width = 0.2,
                                  points_size = 1,
                                  points_alpha = 0.5,
                                  progress = TRUE)



