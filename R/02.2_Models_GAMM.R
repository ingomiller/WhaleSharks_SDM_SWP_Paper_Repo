#_____________________________________________________________________________
#                        Models: Tracking - GAMM
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
library(mgcv)
library(patchwork)
library(MuMIn)
source("R/00_Helper_Functions.R") 





# Build Model -------------------------------------------------------------


model_dt <- readRDS( "MODEL_DATA.rds")



set.seed(6666)


weights_PA <- ifelse(model_data$PA == 1, 1, 0.1)
weights_PA

k_general <- 5



# 
fml <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(mltost,   month, bs = c("cr","cr"), k = c(4,3))



sp_blocks <- sb
gamm_cv_out <- gam_cv(
  dat   = model_dt,             
  fid   = sb$folds_ids,                   
  fml   = fml,          
  family = binomial(link = "cloglog"),
  base_preds  = c("thetao", "mltost", "chl", "uv", "wz", "slope", "depth", "dist2000", "month"),
  re_term    = "s(id)"
)


gamm_cv_out
gamm_cv_out$summary

gamm_vi_df <- gamm_cv_out$varimp |>
  tibble::as_tibble() |>
  dplyr::mutate(algorithm = "GAMM")

gamm_vi_df




# Build full gamm


form <- PA ~
  s(thetao,  bs = "tp", k = 6) +
  s(uv,      bs = "cr", k = 5) +
  s(wz,      bs = "cr", k = 5) +
  s(chl,     bs = "cr", k = 5) +
  s(depth,   bs = "cr", k = 5) +
  s(mltost, bs = "cr", k = 5) +
  s(slope,   bs = "cr", k = 5) +
  s(dist2000, bs = "cr", k = 5) +
  # Month interactions for *dynamic* vars (small k, shrinkable)
  ti(thetao,   month, bs = c("tp","cr"), k = c(4,3)) +
  ti(chl,      month, bs = c("cr","cr"), k = c(4,3)) +
  ti(uv,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(wz,       month, bs = c("cr","cr"), k = c(4,3)) +
  ti(mltost,   month, bs = c("cr","cr"), k = c(4,3)) +
  s(id, bs = "re")


gam <- mgcv::bam(form,
                      family = binomial(link="cloglog"), 
                      method="fREML", 
                      control =  mgcv::gam.control(maxit = 500, epsilon = 1e-5, nthreads = 6),
                      discrete = TRUE,
                      weights = weights_PA,
                      na.action = na.fail,
                      data = model_dt,
                      select = TRUE,
                      gamma = 1
)


summary(gam)




# Response curves ---------------------------------------------------------
model <- gam
rug = 4

p.depth <- sjPlot::plot_model(model, type = "pred", terms = "depth", colors = c("turquoise1"),
                              pred.type = "fe",
                              axis.lim =  c(0, 1),
                              show.data = FALSE,
                              show.intercept = FALSE) +
  labs(x = expression("depth (m)"), y = "Probability of Presence", title = NULL) +
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = depth),
    inherit.aes = FALSE,
    sides = "t",                    
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.6
  ) +
  ggplot2::geom_rug(
    data = model_dt |>  dplyr::filter(PA == 1),
    mapping = ggplot2::aes(x = depth),
    inherit.aes = FALSE,
    sides = "b",                   
    outside = FALSE,
    length = grid::unit(rug, "pt"),
    alpha = 0.3
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  theme_bw() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

p.depth




# create extrnal validation data from sightings 

sight <- readRDS("SIGHTINGS.rds")


val_dt <- sight |>
  plyr::mutate(
    chl = log(chl),
    dist2000 = dist2000/1000) |> 
  sf::st_drop_geometry() |> 
  dplyr::filter(!month %in% c(7, 8, 9, 10)) |>
  as.data.frame()






pred_prob <- predict(
  gam,
  newdata = val_dt,
  type    = "response",
  exclude = "s(id)"
)

p <- pred_prob[val_dt$PA == 1]
a <- pred_prob[val_dt$PA == 0]

ev <- dismo::evaluate(p = p, a = a)

# AUC
auc_ext <- ev@auc
auc_ext

roc_obj <- pROC::roc(response = val_dt$PA, predictor = pred_prob, quiet = TRUE, plot = TRUE, col = "steelblue", lwd = 2, legacy.axes = TRUE)
roc_obj
auc_ext <- as.numeric(pROC::auc(roc_obj))
auc_ext

# TSS at the max-sum threshold (sens + spec)
thr     <- dismo::threshold(ev, stat = "spec_sens")  # same operating point SDMtune uses by default
tss_ext <- ev@TPR[which.max(ev@TPR + ev@TNR)] + ev@TNR[which.max(ev@TPR + ev@TNR)] - 1
tss_ext

roc_df <- data.frame(FPR = 1 - ev@TNR, TPR = ev@TPR)
gg_ext <- ggplot2::ggplot(roc_df, ggplot2::aes(FPR, TPR)) +
  ggplot2::geom_line() +
  ggplot2::geom_abline(linetype = 2) +
  ggplot2::coord_equal() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = sprintf("External validation (AUC = %.2f)", auc_ext),
                x = "False positive rate", y = "True positive rate")

print(gg_ext)


cat(sprintf("External AUC = %.3f | TSS = %.3f | Threshold = %.3f\n", auc_ext, tss_ext, thr))





# Monthlky predictions ----------------------------------------------------


model <- gam
model
relevant_vars <- all.vars(model$formula)[-1] # Excluding the intercept
relevant_vars <- relevant_vars[ !(relevant_vars %in% c("k_general")) ]
relevant_vars

# Extract relevant layers from each raster stack in the list

min(tracks$date) 
max(tracks$date) 


# import monhtly predictor stack raster 
input_list <- monthly_stacks_lst_0.1_trans
names(input_list)

input_list <- lapply(input_list, function(stack) {
  stack[[relevant_vars]] # Extract layers that match the model variables
})


dates <- as.Date(names(input_list))

monthly_predictions <- list()
N <- length(input_list)

tictoc::tic("Loop run took: " )
for (i in seq_along(input_list)) {
  # Extract the current monthly raster stack
  monthly_stack <- input_list[[i]]
  
  # Extract the current raster's date"
  Month_Date <- names(input_list[i])
  
  # Print information about the current prediction
  cat("Prediction", i, "of", N, "\n")
  
  # Generate predictions for the current monthly raster stack
  predictions <- terra::predict(monthly_stack,
                                model,
                                # const = (data.frame(Year = "0")),
                                #cores = 6,
                                type = "response"
  )
  
  # Assign Date as layer name:
  names(predictions) <- Month_Date
  
  # Store the predictions in the list
  monthly_predictions[[i]] <- predictions
}
tictoc::toc()

monthly_predictions

monthly_predictions_stack <- terra::rast(monthly_predictions)


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



