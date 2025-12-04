#_____________________________________________________________________________
#                        Models: m
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
library(blockCV)
library(tmap)
source("R/00_Helper_Functions.R")






model_dt <- readRDS( "MODEL_DATA.rds")


swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = model_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = model_dt |>
    dplyr::select(thetao,
                  mltost,
                  chl,
                  uv,
                  wz,
                  slope,
                  depth, 
                  roughness,
                  dist2000,
                  month) |> 
    as.data.frame(),
  pa = model_dt$PA)

swd




# Using Cross Validation Method


set.seed(28)
cv_m0 <- SDMtune::train(method = "BRT", # Repeat for Maxent
                          data = swd,
                          folds = sp_blocks)


cv_m0

model <- cv_m0

cat("Training AUC: ", SDMtune::auc(model))
cat("Testing AUC: ", SDMtune::auc(model, test = TRUE))
cat("Training TSS: ", SDMtune::tss(model))
cat("Testing TSS: ", SDMtune::tss(model, test = TRUE))



# Variable Correlations ---------------------------------------------------




predictors_mean <- terra::rast("PREDICTORS.tif")


bg_dt <- model_dt |> dplyr::filter(PA == 0)
bg4cor <- SDMtune::SWD(
  species = "Bgs",
  coords  = bg_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data    = bg_dt |> dplyr::select(thetao, mltost, chl, uv, wz, 
                                   depth, 
                                   slope,
                                   roughness,
                                   dist2000,
                                   month) |> as.data.frame(),
  pa      = bg_dt$PA)



SDMtune::plotCor(bg4cor, 
                 method = "spearman", 
                 cor_th = NULL)

SDMtune::corVar(bg4cor, 
                method = "spearman", 
                cor_th = 0.6)


cv_m1 <- SDMtune::varSel(cv_m0, 
                           metric = "tss", 
                           test = TRUE, 
                           bg4cor = bg4cor,
                           method = "spearman", 
                           cor_th = 0.6,
                           permut = 10)

cv_m1


cat("Training AUC: ", SDMtune::auc(cv_m1))
cat("Testing AUC: ", SDMtune::auc(cv_m1, test = TRUE))
cat("Training TSS: ", SDMtune::tss(cv_m1))
cat("Testing TSS: ", SDMtune::tss(cv_m1, test = TRUE))

cat("Training AUC: ", SDMtune::auc(cv_m0))
cat("Testing AUC: ", SDMtune::auc(cv_m0, test = TRUE))
cat("Training TSS: ", SDMtune::tss(cv_m0))
cat("Testing TSS: ", SDMtune::tss(cv_m0, test = TRUE))



# Fine tune model  --------------------------------------------------------


SDMtune::getTunableArgs(cv_m1)


# for BRT
h_list <- list(
  n.trees           = c(2500L, 5000L, 10000L, 20000L),
  interaction.depth = c(1L, 2L, 3L, 4L, 5L),
  shrinkage         = c(0.01, 0.005, 0.001),
  bag.fraction      = c(0.5, 0.75),
  distribution      = "bernoulli"
)

# For Maxent
h_list <- list(reg = seq(0.5, 5, by = 0.5), 
          fc = c("l", "lq", "lh", "lqp", "lqph"),
          iter = c(2500, 5000, 7500, 10000)
)



expected_fits(pop = 20, 
              gen = 5, 
              keep_best = 0.4,
              keep_random = 0.2,
              rounding= "round")


cv_m2 <- SDMtune::optimizeModel(cv_m1,
                                    hypers = h_list,
                                    metric = "tss",
                                    test = TRUE,
                                    pop = 20,
                                    gen = 5,
                                    keep_best = 0.4,
                                    keep_random = 0.2,
                                    mutation_chance = 0.4,
                                    interactive = TRUE,
                                    progress = TRUE,
                                    seed = 694
)

cv_m2


cv_m2@results

index <- terra::which.max(cv_m2@results$test_TSS)
# index <- 2

best_cv_m2 <- cv_m2@models[[index]]


m <- best_cv_m2
m


cat("Training AUC: ", SDMtune::auc(m))
cat("Testing AUC: ", SDMtune::auc(m, test = TRUE))
cat("Training TSS: ", SDMtune::tss(m))
cat("Testing TSS: ", SDMtune::tss(m, test = TRUE))


ROCplots <- lapply(seq_len(ncol(m@folds$test)), function(i){
  idx <- m@folds$test[, i]
  test_i <- SDMtune::SWD(
    species = m@data@species,
    coords  = m@data@coords[idx, , drop = FALSE],
    data    = m@data@data  [idx, , drop = FALSE],
    pa      = m@data@pa    [idx]
  )
  SDMtune::plotROC(m@models[[i]], test = test_i) + ggplot2::ggtitle(paste("Fold", i))
})



ROCplots_combined <- wrap_plots(ROCplots) +
  plot_layout(ncol = 3, nrow = 2, axes = "collect") &
  ggplot2::theme(
    legend.position = c(0.95, 0.03),  # x, y coordinates (right–bottom)
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6), colour = NA),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ROCplots_combined





# Variable Importance

m_vi <- SDMtune::varImp(m, permut = 10)

# Make df
m_vi_df <- m_vi |>
  tibble::as_tibble(rownames = "variable") |>
  dplyr::rename(importance = 2) |>
  dplyr::mutate(algorithm = "m")


m_vi_df




# helper to build an SWD from row indices of the original SWD
.build_swd <- function(obj, idx) {
  SDMtune::SWD(
    species = obj@data@species,
    coords  = obj@data@coords[idx, , drop = FALSE],
    data    = obj@data@data  [idx, , drop = FALSE],
    pa      = obj@data@pa    [idx]
  )
}

m <- cv_m0

k <- ncol(m@folds$test)

metrics_per_fold <- data.frame(
  fold      = seq_len(k),
  auc_train = vapply(seq_len(k), function(i) {
    idx_tr <- m@folds$train[, i]
    swd_tr <- .build_swd(m, idx_tr)
    SDMtune::auc(m@models[[i]], test = swd_tr)
  }, numeric(1)),
  auc_test  = vapply(seq_len(k), function(i) {
    idx_te <- m@folds$test[, i]
    swd_te <- .build_swd(m, idx_te)
    SDMtune::auc(m@models[[i]], test = swd_te)
  }, numeric(1)),
  tss_train = vapply(seq_len(k), function(i) {
    idx_tr <- m@folds$train[, i]
    swd_tr <- .build_swd(m, idx_tr)
    SDMtune::tss(m@models[[i]], test = swd_tr)
  }, numeric(1)),
  tss_test  = vapply(seq_len(k), function(i) {
    idx_te <- m@folds$test[, i]
    swd_te <- .build_swd(m, idx_te)
    SDMtune::tss(m@models[[i]], test = swd_te)
  }, numeric(1))
)

# rounded view and fold means (optional)
metrics_per_fold_rounded <- metrics_per_fold |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ base::round(.x, 2)))


cv_summary <- base::list(
  metrics_per_fold |>
    dplyr::summarise(
      dplyr::across(c(auc_train, auc_test, tss_train, tss_test), ~ base::mean(.x))
    ) |>
    dplyr::mutate(stat = "mean"),
  metrics_per_fold |>
    dplyr::summarise(
      dplyr::across(c(auc_train, auc_test, tss_train, tss_test), ~ stats::sd(.x))
    ) |>
    dplyr::mutate(stat = "sd")
) |>
  dplyr::bind_rows() |>
  dplyr::mutate(dplyr::across(-stat, ~ base::round(.x, 3))) |>
  dplyr::relocate(stat)

metrics_per_fold_rounded
cv_summary




# get final model 
set.seed(25)
final_m <- SDMtune::combineCV(m)



# m <- best_cv_m2
SDMtune::plotResponse(m, 
                               var = "thetao", 
                               type = "cloglog", 
                               only_presence = TRUE, 
                               marginal = TRUE, 
                               rug = TRUE,
                               fun = mean,
                               color = "steelblue") +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression("sst ("*degree*"C)")) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))





# create extrnal validation data from sightings 
val_dt <- readRDS("SIGHTINGS.rds")


val.swd <- SDMtune::SWD(
  species = "Rhincodon_typus",
  coords = val_dt |> dplyr::select(lon, lat) |> as.data.frame(),
  data = val_dt |>
    dplyr::select(thetao, mltost, chl, uv, wz, depth, 
                  slope,
                  dist2000, month) |> as.data.frame(),
  pa = val_dt$PA)

val.swd


SDMtune::auc(final_m)
SDMtune::tss(final_m)

SDMtune::auc(final_m, 
             test = val.swd)
SDMtune::tss(final_m, 
             test = val.swd)




# Monthly predictions -----------------------------------------------------


relevant_vars <- names(final_m@data@data)
relevant_vars


unique(final_m@data@data$month)


# Import Monthly Predicotr Stack
input_list <- monthly_stacks_lst_0.1_trans

input_list <- lapply(input_list, function(stack) {
  stack[[relevant_vars]] 
})

plot(input_list[[1]])


# Get start/end from model data
rng <- range(as.Date(model_dt$Date), na.rm = TRUE)
start_mon <- as.Date(format(rng[1], "%Y-%m-01"))
end_mon   <- as.Date(format(rng[2], "%Y-%m-01"))
want <- seq(start_mon, end_mon, by = "1 month")

# Align monthly rasters to that range
avail <- as.Date(names(input_list))  
idx <- match(want, avail)

monthly_predictors   <- input_list[idx[!is.na(idx)]]
monthly_dates  <- avail[idx[!is.na(idx)]]



e <- terra::ext(c(140, 170, -40, 0))
e
tictoc::tic("Monthly dynamic predictions took: " )
monthly_predictions_stack <- sdmtune_predict_monthly(
  model = final_m,
  monthly_list = monthly_predictors,
  dates = monthly_dates, 
  type = "cloglog",
  extent = e,
  verbose = TRUE
)
tictoc::toc()




### Monthly Means over tracking period

# Define list for monthly averages:
monthly_means <- list()

# Loop through each month:
for (i in 1:12){
  
  # Extract the layers corresponding to the current month
  month_layers <- monthly_predictions_stack[[grep(sprintf("%02d", i), names(monthly_predictions_stack))]]
  
  # Check if layers exist for the current month
  if (length(month_layers) > 0) {
    
    # Calculate the mean across all layers for the current month
    #monthly_mean <- calc(month_layers, mean, na.rm = TRUE)
    monthly_mean <- terra::app(month_layers, fun = mean, na.rm = TRUE)
    
    # Assign meaningful names:
    names(monthly_mean) <- paste(month.abb[i], "Mean")
    
    # Print current layer being processed
    print(paste("Calculating mean for", month.abb[i]))
    
    # Store in list
    monthly_means[[i]] <- monthly_mean
  }
}



monthly_means_stack <- terra::rast(monthly_means)
names(monthly_means_stack) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


# plot(monthly_means_stack, range = c(0, 1))




## Overall mean 



mean_climate <- terra::app(monthly_means_stack, fun = mean, na.rm = TRUE)






### Monsoon vs Trade Wind::



seasons_2 <- c("Monsoon", "Dry")

# Define the months corresponding to each season
season_months_2 <- list(
  Monsoon = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr"),
  Dry = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
)


# Initialize an empty stack to store seasonal means
seasonal_means_stack_2 <- terra::rast()

# Loop through each season
for (season in seasons_2) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months_2[[season]]
  season_stack <- monthly_means_stack[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"Mean", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_means_stack_2 <- c(seasonal_means_stack_2, seasonal_mean)
}


seasonal_means_stack_2


plot(seasonal_means_stack_2, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))






seasons <- c("Q1", "Q2", "Q3", "Q4")
season_months <- list(
  Q1 = c("Jan", "Feb", "Mar"),
  Q2 = c("Apr", "May", "Jun"),
  Q3 = c("Jul", "Aug", "Sep"),
  Q4 = c("Oct", "Nov", "Dec"))


# Initialize an empty stack to store seasonal means
seasonal_means_stack_4 <- terra::rast()

# Loop through each season
for (season in seasons) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months[[season]]
  season_stack <- monthly_means_stack[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"Mean", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_means_stack_4 <- c(seasonal_means_stack_4, seasonal_mean)
}




plot(seasonal_means_stack_4, range = c(0, 1), col = colorRampPalette(c("blue4", "dodgerblue2", "cyan2", "green4", "yellow", "orange", "firebrick1"))(100))
