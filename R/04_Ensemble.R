


#_____________________________________________________________________________
#                        Models: Ensemble
#_____________________________________________________________________________


# remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
# citation("SDMtune")
library(tidyverse)
library(raster)
library(terra)
library(ncdf4)
library(rerddap)
library(rerddapXtracto)
library(stringr)
library(lubridate)
library(SDMtune)
library(patchwork)
source("R/00_Helper_Functions.R") 
library(rasterVis)
library(tidyterra)
library(ggspatial)
library(sf)
library(basemaps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(grid)
library(patchwork)




# Import Predictions ------------------------------------------------------


season_2_max <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_Maxent_seasons2_means_rev_mp_crwPA.tif")


season_2_gam <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_GAMM_seasons2_means_rev_mp_crwPA.tif")


season_2_brt <- terra::rast("/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/SDM_whalesharks_Tracks_BRT_seasons2_means_rev_mp_crwPA.tif")




seasonal_mess_stack_2 <- terra::rast( "/Volumes/Ingo_PhD/PhD_Data_Analysis/PhD_WhaleSharks_SDMs_Enviro_Layers/Chapter2/SDM_Outputs_Rev/Mess_map_gamm_track_mp_seasons_2_crwPA.tif")



names(seasonal_mess_stack_2) <- c("Monsoon Season (Nov - Apr)", 
                                  "Dry Season (May - Oct)")
seasonal_mess_stack_2_extrap <- terra::ifel(seasonal_mess_stack_2 < 0, 1, NA_real_)

seasonal_mess_stack_2_extrap_df <- seasonal_mess_stack_2_extrap|>
  terra::as.data.frame(xy = TRUE) |>
  tidyr::pivot_longer(
    cols = -c(x, y),
    names_to  = "lyr",
    values_to = "flag"
  ) |>
  dplyr::filter(!is.na(flag))




# Ensemble  ---------------------------------------------------------------


## using TSS test weights 
weights <- c(0.33, 0.33, 0.28)



season2 <- c(season_2_gam, season_2_max, season_2_brt)


monsoon <- season2[[c(1, 3, 5)]]
dry <- season2[[c(2, 4, 6)]]


monsoon_ensemble <- terra::weighted.mean(monsoon, w = weights)
dry_ensemble <- terra::weighted.mean(dry, w = weights)


season_2_ensemble <- c(monsoon_ensemble, dry_ensemble)
names(season_2_ensemble) <- c("Monsoon_Season_Ensemble_Tracking", "Dry_Season_Ensemble_Tracking")



