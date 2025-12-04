


#_____________________________________________________________________________
#                        Models: MESS
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
library(dismo)
source("R/00_Helper_Functions.R") 




# Import DATA -------------------------------------------------------------


model_input <- readRDS( "MODEL_DATA.rds")

model <- readRDS("FINAL_MODEL_BRT.rds")

# monthly predictor stack
input_list <- monthly_stacks_lst_0.1_trans



# prepare dta -------------------------------------------------------------

# variables used in models 
relevant_vars <- all.vars(model$formula)[-1]
relevant_vars
relevant_vars <- setdiff(relevant_vars, c("month", "id"))
relevant_vars

# reference environment in training data 
env_train <- model_dt[, relevant_vars, drop = FALSE]
env_train
names(env_train)



# MESS computation --------------------------------------------------------

vars <- colnames(env_train)
stopifnot(length(vars) > 0)

compute_mess <- function(spatr, env_tbl) {
  preds_r <- raster::stack(spatr)         
  dismo::mess(preds_r, env_tbl) |>
    terra::rast()
}

# only month present in training data range 

input_list <- lapply(input_list, function(stack) {
  stack[[relevant_vars]] # Extract layers that match the model variables
})

input_list <- input_list[119:186]


names(input_list[[1]])

N <- length(input_list)
mess_list <- vector("list", N)
names(input_list[[1]])
names(mess_list) <- names(input_list)
mess_list




for (i in seq_along(input_list)) {
  nm <- names(input_list)[i]
  if (is.na(nm) || nm == "") nm <- paste0("m", i)
  
  s <- input_list[[i]]
  
  # 0) Sanity: non-empty stack
  if (terra::nlyr(s) == 0L) {
    stop("Monthly stack has 0 layers at index ", i, " (name: ", nm, ").")
  }
  
  # 1) Name alignment
  rn <- names(s)
  missing <- setdiff(vars, rn)
  extra   <- setdiff(rn, vars)
  
  if (length(missing)) {
    stop("Missing layers in month ", nm, ": ", paste(missing, collapse = ", "),
         ". Layers present: ", paste(rn, collapse = ", "))
  }
  # Reorder to match env_train columns (ignore extras if any snuck in)
  s <- s[[vars]]
  
  # 2) Column/layer count must match exactly
  if (NCOL(env_train) != terra::nlyr(s)) {
    stop("Mismatch in counts for month ", nm, ": ",
         "NCOL(env_train)=", NCOL(env_train), " vs nlyr(stack)=", terra::nlyr(s))
  }
  
  cat("MESS", i, "of", N, ":", nm, "\n")
  r <- compute_mess(s, env_train)
  names(r) <- nm
  mess_list[[i]] <- r
}



mess_monthly_stack <- terra::rast(mess_list)

mess_monthly_stack[is.infinite(mess_monthly_stack)] <- NA

mess_monthly_stack

# Monthly means 
dates <- as.Date(names(mess_monthly_stack))
monthly_means <- vector("list", 12)

for (i in 1:12) {
  # pick layers for this calendar month
  idx <- which(format(dates, "%m") == sprintf("%02d", i))
  if (length(idx) > 0) {
    month_layers <- mess_monthly_stack[[idx]]
    monthly_mean <- terra::app(month_layers, fun = mean, na.rm = TRUE)
    names(monthly_mean) <- paste(month.abb[i])
    print(paste("Calculating mean for", month.abb[i]))
    monthly_means[[i]] <- monthly_mean
  }
}

# Combine into a 12-layer stack (skips NULLs if any months missing)
mess_mean_calendar <- terra::rast(Filter(Negate(is.null), monthly_means))



plot(mess_mean_calendar, range = c(-300, 0))



mess_mean_climate <- terra::app(mess_mean_calendar, fun = mean, na.rm = TRUE)

seasons_2 <- c("Monsoon", "Dry")

# Define the months corresponding to each season
season_months_2 <- list(
  Monsoon = c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr"),
  Dry = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
)


# Initialize an empty stack to store seasonal means
seasonal_mess_stack_2 <- terra::rast()

# Loop through each season
for (season in seasons_2) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months_2[[season]]
  season_stack <- mess_mean_calendar[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"MESS", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_mess_stack_2 <- c(seasonal_mess_stack_2, seasonal_mean)
}

seasonal_mess_stack_2

plot(seasonal_mess_stack_2, range = c(-300, 0))


seasons_4 <- c("Q1", "Q2", "Q3", "Q4")
season_months_4 <- list(
  Q1 = c("Jan", "Feb", "Mar"),
  Q2 = c("Apr", "May", "Jun"),
  Q3 = c("Jul", "Aug", "Sep"),
  Q4 = c("Oct", "Nov", "Dec"))


# Initialize an empty stack to store seasonal means
seasonal_mess_stack_4 <- terra::rast()


# Loop through each season
for (season in seasons_4) {
  # Subset the monthly means stack to include only the layers corresponding to the current season
  season_layers <- season_months_4[[season]]
  season_stack <- mess_mean_calendar[[season_layers]]
  
  # Calculate the mean across the layers for the current season
  seasonal_mean <- terra::app(season_stack, fun = mean, na.rm = TRUE)
  
  # Assign a meaningful name to the seasonal mean
  names(seasonal_mean) <- paste(season,"MESS", sep = "_")
  
  # Add the seasonal mean raster to the stack
  seasonal_mess_stack_4 <- c(seasonal_mess_stack_4, seasonal_mean)
}

seasonal_mess_stack_4

plot(seasonal_mess_stack_4, range = c(-300, 0))



