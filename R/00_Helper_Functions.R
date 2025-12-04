__________________________________________________________________
#                          HELPER FUNCTIONS 
#_____________________________________________________________________________







# CMEMS extractions  ------------------------------------------------------

## Function to extract surface layer values from CMEMS files

CMEMS_extract_by_date_surf <- function(rstack, sf_points, date_col = "Date", varname = "value", buffer_m = 10000) {
  stopifnot(inherits(rstack, "SpatRaster"))
  stopifnot(inherits(sf_points, "sf"))
  stopifnot(date_col %in% names(sf_points))
  
  # 1) normalise dates to match layer names ("YYYY-MM-DD")
  loc_dates_chr <- format(as.POSIXct(sf_points[[date_col]], tz = "UTC"), "%Y-%m-%d")
  lyr_names     <- names(rstack)
  
  # 2) find layer index for each row
  lyr_idx <- match(loc_dates_chr, lyr_names)
  
  if (all(is.na(lyr_idx))) {
    stop("No dates in locations matched any raster layer names", call. = FALSE)
  }
  if (any(is.na(lyr_idx))) {
    warning("Some locations had no matching date in the raster stack; returning NA for those rows.", call. = FALSE)
  }
  
  # 3) ensure geometries align
  v <- terra::vect(sf_points)
  if (!terra::same.crs(v, rstack)) {
    v <- terra::project(v, rstack)
  }
  
  # 4) extract per-layer
  vals <- rep(NA_real_, nrow(sf_points))
  idx_split <- split(seq_len(nrow(sf_points)), lyr_idx, drop = TRUE)
  
  buffer_m <- as.numeric(buffer_m)
  
  for (k in names(idx_split)) {
    li  <- as.integer(k)
    ids <- idx_split[[k]]
    vi  <- v[ids]
    val <- terra::extract(rstack[[li]], vi, method = "bilinear", search_radius = buffer_m)[, 2]
    vals[ids] <- val
  }
  
  # 5) append with dynamic column name
  sf_points[[varname]] <- vals
  sf_points
}




# Bluelink BRAN2020 -------------------------------------------------------

windDir <-function(u, v){
  (270 - atan2(u, v) * 180 / pi) %% 360 
}

windSpd <- function(u, v){
  sqrt(u^2 + v^2)
}

dataDownload <- function(type, year, month = NULL, dir, varname, quiet = TRUE) {
  if (type %in% c("day", "month") == FALSE) {
    stop("Please provide type as 'day' or 'month'.")
  }
  if (is.null(year)) {
    stop("Please provide a year of interest to download data.")
  }  
  if(!varname %in% c('ocean_temp', 'ocean_salt', 'ocean_u', 'ocean_v', 'ocean_w', 'ocean_eta_t', 'ocean_mld', 'atm_flux_diag')){
    stop("Environmental variable not recognised, options include:\n'ocean_temp', 'ocean_salt', 'ocean_u', 'ocean_v', 'ocean_w', 'ocean_eta_t', 'ocean_mld', 'atm_flux_diag'")}
  # Change variable name to download correct netCDF file from BRAN
  if (varname %in% c('air_wind')) {
    varname <- "atm_flux_diag"
  }
  # Increase timeout for large file downloads and slow internet connexions
  options(timeout = 1000000000) 
  # Daily resolution
  if (type == "day") {
    if (year < 2023) {
      if (is.null(month)) {
        month <- 1:12
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      } else {
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      }      
    } else {
      if (is.null(month)) {
        month <- 1:12 # 2023 data now goes to december!
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      } else {
        for (i in 1:length(month)) {
          if (month[i] < 10) {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_0", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          } else {
            download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/daily", varname, "_", year, "_", month[i], ".nc"), 
                          destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
            gc()
          }
        }
      }
    }
  }
  # Monthly resolution
  if (type == "month") {
    if (is.null(month)) {
      month <- 1:12
      for (i in 1:length(month)) {
        if (month[i] < 10) {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_0", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        } else {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        }
      }
    } else {
      for (i in 1:length(month)) {
        if (month[i] < 10) {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_0", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_0", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        } else {
          download.file(paste0("https://thredds.nci.org.au/thredds/dodsC/gb6/BRAN/BRAN2020/month/", varname, "_mth_", year, "_", month[i], ".nc"), 
                        destfile = paste0(dir, "/", varname, "_", year, "_", month[i], ".nc"), mode = 'wb', quiet = quiet)
          gc()
        }
      }
    }
  }   
}   







extractWz <- function(df,
                      X,
                      Y,
                      datetime,
                      folder_name,                 # directory with files like ocean_w_YYYY_MM.nc
                      export_path = NULL,
                      max_depth = -200,
                      fill_gaps = TRUE,
                      buffer = 10000,             # metres (geodesic with s2)
                      verbose = TRUE,
                      export_step = TRUE) {
  
  # ---- checks ---------------------------------------------------------------
  if (!base::dir.exists(folder_name)) {
    base::stop(base::sprintf("NetCDF directory not found: %s", folder_name), call. = FALSE)
  }
  if (!("Depth" %in% base::names(df))) {
    base::warning("Depth column not found; returning input df unchanged.")
    return(df)
  }
  
  
  # helpers for reading NC time units and converting offsets -> POSIXct
  .read_nc_time_info <- function(nc_path) {
    nc <- ncdf4::nc_open(nc_path)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    time_name <- base::intersect(base::names(nc$dim), c("Time", "time", "TIME"))
    if (base::length(time_name) == 0L) time_name <- "time"
    units_att <- try(ncdf4::ncatt_get(nc, time_name, "units")$value, silent = TRUE)
    if (base::inherits(units_att, "try-error") || base::is.null(units_att)) units_att <- ""
    list(units = units_att)
  }
  
  .convert_offsets_to_datetime <- function(offset_numeric, units_att) {
    # examples: "days since 1990-01-01 00:00:00", "hours since 1979-01-01"
    if (!base::grepl("since", units_att, ignore.case = TRUE)) {
      origin <- base::as.POSIXct("1979-01-01 00:00:00", tz = "UTC")
      return(origin + base::as.difftime(offset_numeric, units = "days"))
    }
    unit_part <- base::tolower(base::sub(" since .*", "", units_att))
    origin_str <- base::sub(".*since\\s+", "", units_att)
    origin     <- base::as.POSIXct(origin_str, tz = "UTC")
    if (base::is.na(origin)) origin <- base::as.POSIXct(base::paste0(origin_str, " 00:00:00"), tz = "UTC")
    
    if (unit_part %in% c("day", "days")) {
      origin + base::as.difftime(offset_numeric, units = "days")
    } else if (unit_part %in% c("hour", "hours")) {
      origin + base::as.difftime(offset_numeric, units = "hours")
    } else if (unit_part %in% c("sec", "second", "seconds")) {
      origin + base::as.difftime(offset_numeric, units = "secs")
    } else {
      origin + base::as.difftime(offset_numeric, units = "days")
    }
  }
  
  
  
  # clamp depths and prepare aux.date (YYYY-MM)
  df <- df |>
    dplyr::mutate(
      Depth_Wz = base::ifelse(Depth > -5, -6, base::ifelse(Depth < max_depth, max_depth, Depth)),
      Depth_Wz = base::round(Depth_Wz, 1),
      aux.date = base::substr(base::as.character(.data[[datetime]]), 1, 7)
    )
  
  if ("Wz" %in% base::names(df)) {
    if (verbose) base::message("Previous Wz values found. Continuing only for rows with NA Wz...")
    dates <- base::unique(df$aux.date[base::is.na(df$Wz)])
  } else {
    df$Wz <- base::as.numeric(NA)
    dates <- base::unique(df$aux.date)
  }
  if (base::length(dates) == 0L) {
    if (verbose) base::message("Nothing to do (no dates needing computation).")
    return(df |> dplyr::select(-.data$aux.date))
  }
  
  # progress bar
  if (verbose) {
    base::message("Processing Wz from existing NetCDF files...")
    pb <- utils::txtProgressBar(min = 0, max = base::length(dates), initial = 0, style = 3, width = 60)
  }
  
  # main loop over months present in df
  for (i in base::seq_along(dates)) {
    aux.date <- dates[i]
    yyyy <- base::substr(aux.date, 1, 4)
    mm   <- base::substr(aux.date, 6, 7)
    
    nc_file <- base::file.path(folder_name, base::sprintf("ocean_w_%s_%s.nc", yyyy, mm))
    if (!base::file.exists(nc_file)) {
      base::stop(base::sprintf("Missing NetCDF for %s: %s", aux.date, nc_file), call. = FALSE)
    }
    
    # read NC as SpatRaster
    nc.bran <- terra::rast(nc_file)
    
    nm <- terra::names(nc.bran)
    
    nm_mat <- stringr::str_split_fixed(
      nm |>
        stringr::str_remove("^w_sw_ocean=") |>
        stringr::str_replace("_Time=", "_"),
      pattern = "_",
      n = 2
    )
    
    aux.depth <- base::as.numeric(nm_mat[, 1])   
    aux.time  <- base::as.numeric(nm_mat[, 2])   
    
    
    time_info <- .read_nc_time_info(nc_file)
    meta_time <- .convert_offsets_to_datetime(aux.time, time_info$units)
    
    # meta table: depth, POSIXct time, day, and exact layer index
    meta <- dplyr::tibble(
      Depth = -1 * aux.depth,
      Time  = meta_time,
      Date  = base::as.Date(meta_time),
      lyr   = base::seq_along(aux.depth)  # exact layer ids in nc.bran
    )
    # rows for this month
    idx_month <- which(df$aux.date == aux.date)
    
    for (ii in base::seq_along(idx_month)) {
      r <- idx_month[ii]
      dt_i <- base::as.Date(df[[datetime]][r])
      
      # match by calendar day 
      sel <- meta[meta$Date == dt_i, , drop = FALSE]
      if (!base::nrow(sel)) { df$Wz[r] <- 0; next }
      
      # keep layers down to local bottom (Depth_Wz is negative)
      sel <- sel[sel$Depth >= df$Depth_Wz[r], , drop = FALSE]
      if (!base::nrow(sel)) { df$Wz[r] <- 0; next }
      
      # sort shallow (near 0) to deeper 
      sel <- sel[base::order(sel$Depth, decreasing = TRUE), , drop = FALSE]
      
      # extract W at each retained layer (unchanged behaviour; per-layer, with optional buffer)
      Wvals <- base::rep(NA_real_, base::nrow(sel))
      for (jj in base::seq_len(base::nrow(sel))) {
        aux.layer <- nc.bran[[sel$lyr[jj]]]
        val <- terra::extract(x = aux.layer, y = df[r, c(X, Y)])[1, 2]
        
        if (base::is.na(val) && fill_gaps) {
          pos_sf <- df[r, c(X, Y)] |>
            sf::st_as_sf(coords = base::c(1, 2), crs = 4326, remove = FALSE)
          pos_sf <- suppressMessages(sf::st_buffer(pos_sf, buffer))
          buf_vals <- terra::extract(x = aux.layer, y = pos_sf)
          if (base::is.data.frame(buf_vals) && base::ncol(buf_vals) >= 2) {
            base::names(buf_vals)[2] <- "Buffer"
            val <- buf_vals |>
              dplyr::group_by(.data$ID) |>
              dplyr::summarise(Buffer_mean = base::mean(.data$Buffer, na.rm = TRUE)) |>
              dplyr::pull(.data$Buffer_mean)
            if (base::is.nan(val)) val <- NA_real_
          }
        }
        Wvals[jj] <- val
      }
      
      # layer thickness 
      Height <- base::c(0 - sel$Depth[1], sel$Depth[-base::nrow(sel)] - sel$Depth[-1])
      
      # weighted mean upward velocity 
      if (base::all(base::is.na(Wvals)) || base::sum(Height, na.rm = TRUE) == 0) {
        aux.wz <- 0
      } else {
        Each   <- Wvals * Height
        aux.wz <- base::sum(Each, na.rm = TRUE) / base::sum(Height, na.rm = TRUE)
        if (base::is.na(aux.wz)) aux.wz <- 0
      }
      df$Wz[r] <- aux.wz
      
      
    }
    
    if (base::isTRUE(export_step) && !base::is.null(export_path)) {
      utils::write.csv(df |> dplyr::select(-.data$aux.date),
                       base::paste0(export_path, ".csv"),
                       row.names = FALSE)
    }
    
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  
  if (verbose) base::close(pb)
  df |> dplyr::select(-.data$aux.date)
}




### This one for Copernicus Marine vertical current data to complement 2024 data to Bluelink BRAN 2020 data which will not be updated anymore until new product released

## Note: The NetCDF files need to be downloaded from Copernicus Mairne sepearately in Python using the Copernicus Marine Toolbox

# Everything else in this function is in accordance with the Bluelink function 

extractWz_CM <- function(df, X, Y, datetime, bathy_path, folder_name, export_path, max_depth = -200, 
                         fill_gaps = TRUE, buffer = 10000, verbose = TRUE, export_step = TRUE, keep_nc_files = TRUE) {
  
  if (!dir.exists(folder_name)) {
    base::stop(base::sprintf("NetCDF directory not found: %s", folder_name), call. = FALSE)
  }
  if (!("Depth" %in% base::names(df))) {
    base::warning("Depth column not found; returning input df unchanged.")
    return(df)
  }
  
  df <- df |>
    dplyr::mutate(
      Depth_Wz = base::ifelse(Depth > -5, -6, base::ifelse(Depth < max_depth, max_depth, Depth)),
      Depth_Wz = base::round(Depth_Wz, 1),
      aux.date = base::substr(base::as.character(.data[[datetime]]), 1, 7)
    )
  
  if (dir.exists(folder_name) == FALSE) 
    dir.create(folder_name)
  df$aux.date <- substr(df[,which(names(df) == datetime)], 1, 7)
  
  if ("Wz" %in% names(df)) {
    message("Previous Wz calculations found. Continuing...")
    dates <- unique(df$aux.date[is.na(df$Wz)])  
  } else {
    dates <- unique(df$aux.date)
    df$Wz <- NA
  }
  
  nc_files <- list.files(folder_name, pattern = "\\.nc$", full.names = TRUE)
  
  # cache all rasters once
  nc_idx <- lapply(nc_files, function(nc_file) {
    r <- terra::rast(nc_file, subds = "ocean_w")
    if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- "OGC:CRS84"
    list(
      file = nc_file,
      rast = r,
      tvec = as.Date(terra::time(r), tz = "UTC"),  # per-layer timestamps (length = nlyr)
      dvec = as.numeric(terra::depth(r))           # matching depth per layer (length = nlyr)
    )
  })
  
  if (verbose) {
    message("Processing Wz data")
    pb <- txtProgressBar(min = 0, max = nrow(df), initial = 0, style = 3, width = 60)
  }
  
  for (ri in seq_len(nrow(df))) {
    
    # skip if already computed
    if (!is.na(df$Wz[ri])) { if (verbose) setTxtProgressBar(pb, ri); next }
    
    dt_i     <- as.Date(df[[datetime]][ri])
    bottom_m <- abs(df$Depth_Wz[ri])
    found    <- FALSE
    
    # try each file until we find the date
    for (k in seq_along(nc_idx)) {
      idx <- nc_idx[[k]]
      
      lyr_day <- which(idx$tvec == dt_i)
      if (!length(lyr_day)) next
      
      # order depths shallow->deep for that day and keep to bottom_m
      d_day <- idx$dvec[lyr_day]
      o     <- order(d_day, na.last = NA)
      lyr_day <- lyr_day[o]
      d_day   <- d_day[o]
      
      keep <- which(d_day <= bottom_m)
      if (!length(keep)) next
      
      L      <- lyr_day[keep]
      d_keep <- d_day[keep]
      
      # extract W at each kept layer (with optional buffer fill)
      Wvals <- rep(NA_real_, length(L))
      for (j in seq_along(L)) {
        rlyr <- idx$rast[[L[j]]]
        val  <- terra::extract(rlyr, df[ri, c(X, Y)])[1, 2]
        
        if (is.na(val) && fill_gaps) {
          pt  <- sf::st_as_sf(df[ri, c(X, Y)], coords = c(1, 2), crs = 4326, remove = FALSE)
          buf <- suppressMessages(suppressWarnings(sf::st_buffer(pt, buffer)))  # degrees by default
          ex  <- terra::extract(rlyr, buf)
          if (is.data.frame(ex) && ncol(ex) >= 2) {
            val <- mean(ex[[2]], na.rm = TRUE)
            if (is.nan(val)) val <- NA_real_
          }
        }
        Wvals[j] <- val
      }
      
      # positive layer thicknesses (m), assuming d_keep is ascending
      Height <- c(d_keep[1], diff(d_keep))
      
      if (all(is.na(Wvals)) || sum(Height, na.rm = TRUE) == 0) {
        df$Wz[ri] <- NA_real_
      } else {
        df$Wz[ri] <- sum(Wvals * Height, na.rm = TRUE) / sum(Height, na.rm = TRUE)
      }
      
      found <- TRUE
      break  # handled this row; go to next row
    }
    
    if (!found && verbose) {
      message("No matching date in any file for row ", ri, " (", dt_i, ").")
    }
    
    if (export_step && (ri %% 50 == 0 || ri == nrow(df))) {
      utils::write.csv(df[, -which(names(df) == "aux.date")], paste0(export_path, ".csv"), row.names = FALSE)
    }
    
    if (verbose) setTxtProgressBar(pb, ri)
  }
  
  if (verbose) close(pb)
  
  if (!keep_nc_files)
    invisible(file.remove(folder_name))
  return(df)
}




create_monthly_rasters_from_files_CM <- function(nc_folder, desired_depth, reference_raster = NULL, output_path, output_filename, start_date = NULL) {
  
  # Function to process auxiliary data from nc file
  process_aux_data_CM <- function(nc_file) {
    # Load the NetCDF file
    nc_CM <- terra::rast(nc_file)
    gc()  # Garbage collection to free up memory
    
    # Extract and process variable names
    aux_names <- stringr::str_remove(names(nc_CM), pattern = "wo_depth=")
    aux_names <- stringr::str_split(aux_names, pattern = "_")
    
    aux_depth <- NULL
    aux_time <- NULL
    for (i in 1:length(aux_names)) {
      aux_depth <- c(aux_depth, aux_names[[i]][1])
      aux_time <- c(aux_time, aux_names[[i]][2])
    }
    
    aux_names_df <- data.frame(Depth = as.numeric(aux_depth), Time = as.numeric(aux_time))
    
    # start_date handling
    # If start_date is NULL, keep original origin of "2024-01-16"
    origin <- if (is.null(start_date)) as.Date("2024-01-16", tz = "UTC") else as.Date(start_date, tz = "UTC")
    idx <- as.integer(aux_names_df$Time)  
    date_seq <- seq.Date(from = origin, by = "1 month", length.out = max(idx, na.rm = TRUE))
    aux_names_df$Time <- date_seq[idx]
    # ----------------------------------------------------------------
    
    gc()  
    return(aux_names_df)
  }
  
  # List all NetCDF files in the directory
  nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    stop("No NetCDF files found in the specified folder.")
  }
  
  all_layers <- list()
  available_depths <- NULL
  
  for (nc_file in nc_files) {
    message(paste("Processing file:", nc_file))
    aux_data <- process_aux_data_CM(nc_file)
    available_depths <- unique(c(available_depths, aux_data$Depth))
  }
  
  # Find the nearest available depth
  differences <- abs(available_depths - desired_depth)
  nearest_depth <- available_depths[which.min(differences)]
  message(paste("Nearest available depth layer selected:", nearest_depth))
  
  for (nc_file in nc_files) {
    aux_data <- process_aux_data_CM(nc_file)
    
    # Filter by the nearest depth layer
    depth_filter <- aux_data$Depth == nearest_depth
    if (!any(depth_filter)) {
      message(paste("No matching depth layer found in file:", nc_file))
      next  # Skip this file if no layers match the specified depth
    }
    
    filtered_data <- aux_data[depth_filter, ]
    
    # Load the NetCDF file and filter by the selected layers
    nc_CM <- terra::rast(nc_file)
    filtered_layers <- nc_CM[[which(depth_filter)]]
    gc()  
    
    for (i in 1:nrow(filtered_data)) {
      current_layer <- filtered_layers[[i]]
      names(current_layer) <- as.character(filtered_data$Time[i])
      all_layers <- c(all_layers, list(current_layer))
      gc() 
    }
  }
  
  if (length(all_layers) == 0) {
    stop("No layers were extracted. Please check the depth layer and NetCDF files.")
  }
  
  combined_raster <- rast(all_layers)
  
  # Resample and crop the combined raster 
  if (!is.null(reference_raster)) {
    combined_raster <- terra::resample(combined_raster, reference_raster)
    combined_raster <- terra::crop(combined_raster, reference_raster)
  }
  
  # Construct the full output file path 
  output_file <- file.path(output_path, paste0(output_filename, ".tif"))
  
  # Save
  terra::writeRaster(combined_raster, output_file, overwrite = TRUE)
  message(paste("Saved combined raster stack to:", output_file))
  
  assign(output_filename, terra::rast(output_file), envir = .GlobalEnv)
  
  gc()
}












# Downlaod MUR SST monthly global data   ----------------------------------


download_mursst41_monthly <- function(from_year = 2005L,
                                      out_dir,
                                      until_yearmon = NULL,   # e.g., "2025-08"; NULL = last full month
                                      retries = 4L,
                                      backoff_sec = 2,
                                      progress = c("cli","text","none"),
                                      write_log = TRUE) {
  if (base::missing(out_dir) || base::is.null(out_dir)) {
    base::stop("Please provide 'out_dir'.", call. = FALSE)
  }
  if (!base::requireNamespace("curl", quietly = TRUE)) {
    base::stop("Please install.packages('curl') first.", call. = FALSE)
  }
  
  progress <- base::match.arg(progress)
  
  has_cli <- base::requireNamespace("cli", quietly = TRUE)
  if (progress == "cli" && !has_cli) progress <- "text"
  
  if (!base::dir.exists(out_dir)) base::dir.create(out_dir, recursive = TRUE)
  
  # Determine last full month
  last_day <- if (base::is.null(until_yearmon)) {
    first_this_month <- base::as.Date(base::format(base::Sys.Date(), "%Y-%m-01"))
    first_this_month - 1L
  } else {
    ym <- base::as.Date(paste0(until_yearmon, "-01"))
    seq_to_next <- base::seq.Date(ym, by = "month", length.out = 2L)
    seq_to_next[2L] - 1L
  }
  
  # Month sequence
  start_date  <- base::as.Date(sprintf("%d-01-01", from_year))
  first_month <- base::as.Date(base::format(start_date, "%Y-%m-01"))
  last_month  <- base::as.Date(base::format(last_day,   "%Y-%m-01"))
  months_seq  <- base::seq.Date(from = first_month, to = last_month, by = "month")
  if (base::length(months_seq) == 0L) {
    base::message("No months to download for the requested range.")
    return(invisible(base::data.frame()))
  }
  
  # Start/end per month
  next_months <- base::seq.Date(from = months_seq[1L], by = "month", length.out = base::length(months_seq) + 1L)
  month_start <- months_seq
  month_end   <- next_months[-1L] - 1L
  
  # Filenames + URLs
  make_fname <- function(s, e) {
    base::sprintf("%s%s-GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc",
                  base::format(s, "%Y%m%d"), base::format(e, "%Y%m%d"))
  }
  fnames  <- base::mapply(make_fname, month_start, month_end, USE.NAMES = FALSE)
  base_url <- "https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/"
  urls    <- base::paste0(base_url, fnames)
  dest    <- base::file.path(out_dir, fnames)
  
  # Skip existing files
  need <- !base::file.exists(dest)
  urls <- urls[need]; dest <- dest[need]
  month_start <- month_start[need]; month_end <- month_end[need]
  fnames <- fnames[need]
  
  if (base::length(urls) == 0L) {
    base::message("All files already present.")
    return(invisible(base::data.frame()))
  }
  
  # Helpers
  fmt_bytes <- function(n) {
    if (!base::is.finite(n) || base::is.na(n)) return("NA")
    units <- c("B","kB","MB","GB","TB")
    i <- base::floor(base::log(max(n, 1), 1024))
    i <- max(0, min(i, base::length(units)-1))
    base::sprintf("%.1f %s", n / (1024^i), units[i+1])
  }
  
  n <- base::length(urls)
  results <- base::data.frame(
    file = dest,
    month_start = month_start,
    month_end = month_end,
    status = base::rep(NA_character_, n),
    size_bytes = base::rep(NA_real_, n),
    attempt = base::rep(NA_integer_, n),
    stringsAsFactors = FALSE
  )
  
  # Progress
  if (progress == "cli") {
    cli::cli_progress_bar("Downloading MUR monthly files",
                          total = n, clear = FALSE, format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} {cli::pb_eta}")
  }
  
  # Download loop with retries + friendly status lines
  for (i in base::seq_len(n)) {
    ok <- FALSE
    last_err <- ""
    for (attempt in 0L:retries) {
      res <- base::tryCatch({
        curl::curl_download(url = urls[i], destfile = dest[i], mode = "wb")
        TRUE
      }, error = function(e) {
        last_err <<- base::conditionMessage(e)
        FALSE
      })
      if (res) { ok <- TRUE; results$attempt[i] <- attempt + 1L; break }
      base::Sys.sleep(backoff_sec * (2 ^ attempt))
    }
    
    # size and status
    if (ok) {
      sz <- base::file.info(dest[i])$size
      results$status[i] <- "ok"
      results$size_bytes[i] <- if (base::is.na(sz)) NA_real_ else sz
      line <- base::sprintf("[%d/%d] %s : ok (%s)",
                            i, n, base::basename(dest[i]), fmt_bytes(results$size_bytes[i]))
    } else {
      results$status[i] <- paste0("error: ", last_err)
      line <- base::sprintf("[%d/%d] %s : FAILED (%s)",
                            i, n, base::basename(dest[i]), last_err)
    }
    
    if (progress == "cli") {
      cli::cli_progress_update()
      # also echo the line so you get the nice per-file status you liked
      cli::cli_inform(line)
    } else if (progress == "text") {
      base::message(line)
    }
  }
  
  if (progress == "cli") cli::cli_progress_done()
  
  # Final summary
  n_ok <- base::sum(results$status == "ok", na.rm = TRUE)
  n_fail <- base::sum(!base::startsWith(results$status, "ok"), na.rm = TRUE)
  total_sz <- fmt_bytes(base::sum(results$size_bytes[results$status == "ok"], na.rm = TRUE))
  
  base::message("")
  base::message("===== MUR monthly download summary =====")
  base::message(base::sprintf("Saved to: %s", out_dir))
  base::message(base::sprintf("OK: %d | Failed: %d | Total size: %s", n_ok, n_fail, total_sz))
  if (n_fail > 0) {
    base::message("Failed files:")
    base::message(paste0("  - ", base::basename(results$file[!base::startsWith(results$status, "ok")])), sep = "\n")
  }
  
  # Optional CSV log
  if (write_log) {
    ts <- base::format(base::Sys.time(), "%Y%m%d_%H%M%S")
    log_path <- base::file.path(out_dir, paste0("mursst_download_log_", ts, ".csv"))
    utils::write.csv(results, log_path, row.names = FALSE)
    base::message(base::sprintf("Log written: %s", log_path))
  }
  
  results
}







# Predictor Raster Functions  ---------------------------------------------



# Extract the month from the layer names
extract_month_from_names <- function(spat_raster) {
  # Extract the names of the layers
  layer_names <- names(spat_raster)
  
  # Extract the month part assuming the format "YYYY-MM-DD"
  months <- as.integer(substring(layer_names, 6, 7))
  
  # Handle potential NAs and return the result
  if (any(is.na(months))) {
    warning("Some month values could not be coerced to integers. Check the layer names format.")
  }
  
  return(months)
}




## make monthly raster stacks from CMEMS Wo dowanloads for Wo variable 

create_monthly_rasters_from_files_CM <- function(nc_folder, desired_depth, reference_raster = NULL, output_path, output_filename, start_date = NULL) {
  
  # Function to process auxiliary data from nc file
  process_aux_data_CM <- function(nc_file) {
    # Load the NetCDF file
    nc_CM <- terra::rast(nc_file)
    gc()  # Garbage collection to free up memory
    
    # Extract and process variable names
    aux_names <- stringr::str_remove(names(nc_CM), pattern = "wo_depth=")
    aux_names <- stringr::str_split(aux_names, pattern = "_")
    
    aux_depth <- NULL
    aux_time <- NULL
    for (i in 1:length(aux_names)) {
      aux_depth <- c(aux_depth, aux_names[[i]][1])
      aux_time <- c(aux_time, aux_names[[i]][2])
    }
    
    aux_names_df <- data.frame(Depth = as.numeric(aux_depth), Time = as.numeric(aux_time))
    
    # start_date handling (no structural changes elsewhere) 
    # If start_date is NULL, keep original origin of "2024-01-16"
    origin <- if (is.null(start_date)) as.Date("2024-01-16", tz = "UTC") else as.Date(start_date, tz = "UTC")
    idx <- as.integer(aux_names_df$Time)  # 1,2,3,... month indices
    date_seq <- seq.Date(from = origin, by = "1 month", length.out = max(idx, na.rm = TRUE))
    aux_names_df$Time <- date_seq[idx]
    
    
    gc()  
    return(aux_names_df)
  }
  
  # List all NetCDF files in the directory
  nc_files <- list.files(nc_folder, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    stop("No NetCDF files found in the specified folder.")
  }
  
  all_layers <- list()
  available_depths <- NULL
  
  for (nc_file in nc_files) {
    message(paste("Processing file:", nc_file))
    aux_data <- process_aux_data_CM(nc_file)
    available_depths <- unique(c(available_depths, aux_data$Depth))
  }
  
  # Find the nearest available depth
  differences <- abs(available_depths - desired_depth)
  nearest_depth <- available_depths[which.min(differences)]
  message(paste("Nearest available depth layer selected:", nearest_depth))
  
  for (nc_file in nc_files) {
    aux_data <- process_aux_data_CM(nc_file)
    
    # Filter by the nearest depth layer
    depth_filter <- aux_data$Depth == nearest_depth
    if (!any(depth_filter)) {
      message(paste("No matching depth layer found in file:", nc_file))
      next  # Skip this file if no layers match the specified depth
    }
    
    filtered_data <- aux_data[depth_filter, ]
    
    # Load the NetCDF file and filter by the selected layers
    nc_CM <- terra::rast(nc_file)
    filtered_layers <- nc_CM[[which(depth_filter)]]
    gc()  
    
    for (i in 1:nrow(filtered_data)) {
      current_layer <- filtered_layers[[i]]
      names(current_layer) <- as.character(filtered_data$Time[i])
      all_layers <- c(all_layers, list(current_layer))
      gc()  
    }
  }
  
  if (length(all_layers) == 0) {
    stop("No layers were extracted. Please check the depth layer and NetCDF files.")
  }
  
  combined_raster <- rast(all_layers)
  
  # Resample and crop the combined
  if (!is.null(reference_raster)) {
    combined_raster <- terra::resample(combined_raster, reference_raster)
    combined_raster <- terra::crop(combined_raster, reference_raster)
  }
  
  # Construct the full output file path with .tif extension
  output_file <- file.path(output_path, paste0(output_filename, ".tif"))
  
  # Save the combined raster stack to the specified output file
  terra::writeRaster(combined_raster, output_file, overwrite = TRUE)
  message(paste("Saved combined raster stack to:", output_file))
  
  assign(output_filename, terra::rast(output_file), envir = .GlobalEnv)
  
  gc()
}



## Calculate mean velocity:

calculate_uv <- function(stack) {
  layers <- nlyr(stack) / 2
  result_stack <- terra::rast(ncol=ncol(stack), nrow=nrow(stack), extent=ext(stack), crs=crs(stack))
  
  for (i in 1:layers) {
    uo <- stack[[i]]
    vo <- stack[[i + layers]]
    uv_layer <- sqrt(uo^2 + vo^2)
    names(uv_layer) <- names(stack)[i]  # Preserve the original date name
    result_stack <- c(result_stack, uv_layer)
    gc() 
  }
  return(result_stack)
}






# SDMtune Predictions -----------------------------------------------------
## this helper enabes to predict overe a stack of monhtly predictor stack to perform dynamic predictions rather than a single static one 

sdmtune_predict_monthly <- function(model,
                                    monthly_list,       # list of SpatRaster (one per month)
                                    dates = NULL,       # Date vector; if NULL we try names(monthly_list)
                                    type = "cloglog",
                                    extent = NULL,      # crop to this extent before predicting
                                    verbose = TRUE) {
  
  stopifnot(inherits(model, "SDMmodel"))
  
  # Required variables (used to subset each stack)
  req_vars <- colnames(model@data@data)
  if (is.null(req_vars) || !length(req_vars)) {
    stop("Couldn't detect model variables from model@data@data.")
  }
  
  # Dates
  if (is.null(dates)) {
    nm <- names(monthly_list)
    if (is.null(nm)) stop("Provide 'dates' or set names(monthly_list) to ISO dates.")
    dates <- as.Date(nm)
    if (any(is.na(dates))) stop("names(monthly_list) must be parseable as Date; otherwise pass 'dates='.")
  }
  if (length(dates) != length(monthly_list))
    stop("Length of 'dates' must match length of 'monthly_list'.")
  
  # Coerce 'extent' to a SpatExtent if provided
  ex <- NULL
  if (!is.null(extent)) {
    if (inherits(extent, "SpatExtent")) {
      ex <- extent
    } else if (inherits(extent, "SpatRaster") || inherits(extent, "SpatVector")) {
      ex <- terra::ext(extent)
    } else if (is.numeric(extent) && length(extent) == 4) {
      ex <- terra::ext(extent)
    } else {
      stop("`extent` must be a SpatExtent, SpatRaster/SpatVector, or numeric c(xmin, xmax, ymin, ymax).")
    }
  }
  
  N <- length(monthly_list)
  preds <- vector("list", N)
  
  for (i in seq_len(N)) {
    r <- monthly_list[[i]]
    if (!inherits(r, "SpatRaster")) stop("monthly_list[[", i, "]] is not a SpatRaster.")
    
    # Crop to extent if requested
    if (!is.null(ex)) {
      r <- terra::crop(r, ex)
      if (terra::ncell(r) == 0) {
        stop(sprintf("Cropped raster for %s has zero cells; check the extent/CRS.",
                     format(dates[i], "%Y-%m-%d")))
      }
    }
    
    # Make sure all required variables are present
    missing <- setdiff(req_vars, names(r))
    if (length(missing)) {
      stop(sprintf("Month %s is missing layers: %s",
                   format(dates[i], "%Y-%m-%d"), paste(missing, collapse = ", ")))
    }
    
    # Subset and keep model order
    r <- r[[req_vars]]
    
    if (verbose) cat("Predicting", format(dates[i], "%Y-%m"), sprintf("(%d/%d)\n", i, N))
    
    pr <- SDMtune::predict(model, data = r, type = type)
    names(pr) <- as.character(dates[i])   # e.g., "2010-01-01"
    preds[[i]] <- pr
  }
  
  out <- terra::rast(preds)
  names(out) <- as.character(dates)
  out
}


# calculate how many models will be run according to tune settings 
expected_fits <- function(pop = 20, gen = 5, keep_best = 0.4, keep_random = 0.2,
                          rounding = c("floor", "round", "ceiling")) {
  rounding <- match.arg(rounding)
  
  # basic checks
  stopifnot(pop > 0, gen >= 0,
            keep_best >= 0, keep_random >= 0,
            keep_best + keep_random <= 1)
  
  rfun <- switch(rounding,
                 floor   = base::floor,
                 round   = base::round,
                 ceiling = base::ceiling)
  
  kept_per_gen <- rfun(pop * keep_best) + rfun(pop * keep_random)
  new_per_gen  <- pop - kept_per_gen
  total_fits   <- pop + gen * new_per_gen
  
  as.integer(total_fits)
}





# define a patched trainBRT that reads weights from data@data$case_w if present ---
patched_trainBRT <- function(data,
                             distribution = "bernoulli",
                             n.trees = 100,
                             interaction.depth = 1,
                             shrinkage = 0.1,
                             bag.fraction = 0.5) {
  # Build the SDMmodel shell
  result <- SDMtune::SDMmodel(data = data)
  
  # Predictors 
  df <- base::cbind(pa = data@pa, data@data)
  
  # weight subsetting for CV
  w <- NULL
  
  if ("month" %in% base::names(data@data)) {
    # Month-wise balancing: in each month, sum(w_pres) ≈ sum(w_abs)
    w <- base::numeric(nrow(df))
    mlev <- base::unique(data@data$month)
    for (m in mlev) {
      idx  <- which(data@data$month == m)
      if (length(idx)) {
        pres <- base::sum(data@pa[idx] == 1L)
        absn <- base::sum(data@pa[idx] == 0L)
        w[idx] <- if (absn > 0) ifelse(data@pa[idx] == 1L, 1, pres/absn) else 1
      }
    }
  } else {
    # Global balancing: sum(w_pres) ≈ sum(w_abs)
    pres <- base::sum(data@pa == 1L)
    absn <- base::sum(data@pa == 0L)
    w <- if (absn > 0) ifelse(data@pa == 1L, 1, pres/absn) else base::rep(1, nrow(df))
  }
  
 
  if (base::anyNA(w)) {
    w[base::is.na(w)] <- 1
  }
  w[w < 0] <- 0
  
  
  # fit gbm with weights
  gbm_fit <- gbm::gbm(
    pa ~ .,
    data = df,
    distribution = distribution,
    n.trees = n.trees,
    interaction.depth = interaction.depth,
    shrinkage = shrinkage,
    bag.fraction = bag.fraction,
    weights = w
  )
  
  # Wrap back into SDMtune classes
  brt_obj <- SDMtune::BRT(
    n.trees = n.trees,
    distribution = distribution,
    interaction.depth = interaction.depth,
    shrinkage = shrinkage,
    bag.fraction = bag.fraction,
    model = gbm_fit
  )
  result@model <- brt_obj
  result
}





gam_cv <- function(dat,
                   fid,
                   fml,
                   family     = binomial(link = "cloglog"),
                   base_preds = NULL,        
                   re_term    = "s(id)") {   
  
  kvals <- sort(unique(fid))
  
  tss_at <- function(y, p, thr){
    y <- as.integer(y == 1)
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & y == 1); tn <- sum(pred == 0 & y == 0)
    fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
    if ((tp + fn) == 0 || (tn + fp) == 0) return(NA_real_)
    sens <- tp / (tp + fn); spec <- tn / (tn + fp)
    sens + spec - 1
  }
  
  ## omission rate at a given threshold
  omission_rate <- function(y, p, thr){
    y <- as.integer(y == 1)
    pred <- as.integer(p >= thr)
    tp <- sum(pred == 1 & y == 1); fn <- sum(pred == 0 & y == 1)
    if ((tp + fn) == 0) return(NA_real_)
    fn / (tp + fn)
  }
  
  ## continuous Boyce index (Hirzel-style)
  boyce_index <- function(pred, pres, nbins = 20){
    pred <- pred[is.finite(pred)]
    pres <- pres[is.finite(pres)]
    if (length(pred) < 2L || length(pres) < 2L) return(NA_real_)
    
    bks <- stats::quantile(pred, probs = seq(0, 1, length.out = nbins + 1),
                           na.rm = TRUE)
    bks[1] <- -Inf
    bks[length(bks)] <- Inf
    mids <- (bks[-1] + bks[-length(bks)]) / 2
    
    h_pres <- graphics::hist(pres, breaks = bks, plot = FALSE)
    h_all  <- graphics::hist(pred, breaks = bks, plot = FALSE)
    
    P <- h_pres$counts
    E <- h_all$counts
    keep <- P > 0 & E > 0
    if (sum(keep) < 2L) return(NA_real_)
    
    P <- P[keep] / sum(P[keep])
    E <- E[keep] / sum(E[keep])
    ratio <- P / E
    mids  <- mids[keep]
    
    stats::cor(mids, ratio, method = "spearman")
  }
  
  models  <- vector("list", length(kvals))
  metrics <- vector("list", length(kvals))
  
  ## container for per-fold predictor importance (percent)
  do_varimp <- !is.null(base_preds) && length(base_preds) > 0L
  if (do_varimp) {
    vi_mat <- matrix(NA_real_,
                     nrow = length(kvals),
                     ncol = length(base_preds),
                     dimnames = list(kvals, base_preds))
  }
  
  for (i in seq_along(kvals)){
    k <- kvals[i]
    tr_idx <- fid != k; te_idx <- fid == k
    train  <- dat[tr_idx, , drop = FALSE]
    test   <- dat[te_idx, , drop = FALSE]
    
    # class-balanced weights on TRAIN subset
    train$w <- ifelse(train$PA == 1, 1, 0.1)
    
    m <- mgcv::bam(
      fml,
      data    = train,
      family  = family,
      method  = "fREML",
      control = mgcv::gam.control(maxit = 500, epsilon = 1e-5, nthreads = 6),
      discrete = TRUE,
      weights  = w,
      na.action = na.fail,
      select    = TRUE,
      gamma     = 1
    )
    
    # predictions
    p_tr <- stats::predict(m, newdata = train, type = "response")
    p_te <- stats::predict(m, newdata = test,  type = "response")
    
    
    # Miller calibration (train and test)
    mc_tr <- modEvA::MillerCalib(
      obs  = train$PA,
      pred = p_tr,
      plot = FALSE
    )
    
    mc_te <- modEvA::MillerCalib(
      obs  = test$PA,
      pred = p_te,
      plot = FALSE
    )
    
    
    # TRAIN metrics 
    tr_auc <- tryCatch(
      as.numeric(pROC::roc(train$PA, p_tr, quiet = TRUE)$auc),
      error = function(e) NA_real_
    )
    co_tr  <- tryCatch(
      pROC::coords(
        pROC::roc(train$PA, p_tr, quiet = TRUE),
        x           = "best",
        best.method = "youden",
        ret         = c("threshold","sensitivity","specificity")
      ),
      error = function(e) c(threshold = NA, sensitivity = NA, specificity = NA)
    )
    thr_tr <- as.numeric(co_tr["threshold"])
    tr_tss <- if (is.na(thr_tr)) NA_real_ else tss_at(train$PA, p_tr, thr_tr)
    
    ## omission rate 
    tr_omr <- if (is.na(thr_tr)) NA_real_ else
      omission_rate(train$PA, p_tr, thr_tr)
    
    ## Boyce 
    tr_boyce <- boyce_index(pred = p_tr, pres = p_tr[train$PA == 1L])
    
    # TEST metrics
    te_auc <- tryCatch(
      as.numeric(pROC::roc(test$PA, p_te, quiet = TRUE)$auc),
      error = function(e) NA_real_
    )
    # (i) Test TSS at TRAIN-chosen threshold (proper CV) 
    te_tss_at_tr <- if (is.na(thr_tr)) NA_real_ else tss_at(test$PA, p_te, thr_tr)
    # (ii) Test-optimal TSS
    co_te  <- tryCatch(
      pROC::coords(
        pROC::roc(test$PA, p_te, quiet = TRUE),
        x           = "best",
        best.method = "youden",
        ret         = c("threshold","sensitivity","specificity")
      ),
      error = function(e) c(threshold = NA, sensitivity = NA, specificity = NA)
    )
    te_tss_opt <- if (any(is.na(co_te))) NA_real_ else
      as.numeric(co_te["sensitivity"] + co_te["specificity"] - 1)
    
    ## omission rate (test) at *test*-optimal threshold (consistent with test_tss)
    te_omr <- if (any(is.na(co_te))) NA_real_ else
      1 - as.numeric(co_te["sensitivity"])
    
    ## Boyce (test)
    te_boyce <- boyce_index(pred = p_te, pres = p_te[test$PA == 1L])
    
    models[[i]] <- m
    metrics[[i]] <- data.frame(
      fold        = k,
      train_auc   = tr_auc,
      test_auc    = te_auc,
      train_tss   = tr_tss,
      test_tss    = te_tss_opt,   # test-optimised TSS
      thr_train   = thr_tr,
      train_omr   = tr_omr,
      test_omr    = te_omr,
      train_boyce = tr_boyce,
      test_boyce  = te_boyce,
      # Miller calibration:
      train_cal_intercept = mc_tr$intercept,
      train_cal_slope     = mc_tr$slope,
      test_cal_intercept  = mc_te$intercept,
      test_cal_slope      = mc_te$slope,
      test_cal_slopeDiff  = mc_te$slopeDiff
    )
    
    ## link-scale CV
    if (do_varimp) {
      model_data <- train
      
      terms_link <- tryCatch(
        stats::predict(
          object  = m,
          newdata = model_data,
          type    = "terms",
          exclude = re_term
        ),
        error = function(e) {
          stats::predict(
            object  = m,
            newdata = model_data,
            type    = "terms"
          )
        }
      )
      
      term_imp <- colMeans(abs(terms_link), na.rm = TRUE)
      
      pattern <- paste0("\\b(", paste(base_preds, collapse = "|"), ")\\b")
      pred_imp <- stats::setNames(numeric(length(base_preds)), base_preds)
      
      for (j in seq_along(term_imp)) {
        term_name    <- names(term_imp)[j]
        vars_in_term <- stringr::str_extract_all(term_name, pattern)[[1]] |> unique()
        if (length(vars_in_term) > 0L) {
          share <- term_imp[[j]] / length(vars_in_term)
          for (v in vars_in_term) {
            pred_imp[[v]] <- pred_imp[[v]] + share
          }
        }
      }
      
      if (sum(pred_imp) > 0) {
        pred_pct <- 100 * pred_imp / sum(pred_imp)
      } else {
        pred_pct <- pred_imp * NA_real_
      }
      
      vi_mat[i, ] <- pred_pct
    }
  }
  
  metrics <- do.call(rbind, metrics)
  
  summary <- within(data.frame(
    train_auc     = mean(metrics$train_auc,   na.rm = TRUE),
    test_auc      = mean(metrics$test_auc,    na.rm = TRUE),
    train_tss     = mean(metrics$train_tss,   na.rm = TRUE),
    test_tss      = mean(metrics$test_tss,    na.rm = TRUE),
    train_omr     = mean(metrics$train_omr,   na.rm = TRUE),
    test_omr      = mean(metrics$test_omr,    na.rm = TRUE),
    train_boyce   = mean(metrics$train_boyce, na.rm = TRUE),
    test_boyce    = mean(metrics$test_boyce,  na.rm = TRUE),
    train_auc_sd  = stats::sd(metrics$train_auc,   na.rm = TRUE),
    test_auc_sd   = stats::sd(metrics$test_auc,    na.rm = TRUE),
    train_tss_sd  = stats::sd(metrics$train_tss,   na.rm = TRUE),
    test_tss_sd   = stats::sd(metrics$test_tss,    na.rm = TRUE),
    train_omr_sd  = stats::sd(metrics$train_omr,   na.rm = TRUE),
    test_omr_sd   = stats::sd(metrics$test_omr,    na.rm = TRUE),
    train_boyce_sd= stats::sd(metrics$train_boyce, na.rm = TRUE),
    test_boyce_sd = stats::sd(metrics$test_boyce,  na.rm = TRUE),
    # Miller calibration means:
    train_cal_intercept = mean(metrics$train_cal_intercept, na.rm = TRUE),
    train_cal_slope     = mean(metrics$train_cal_slope,     na.rm = TRUE),
    test_cal_intercept  = mean(metrics$test_cal_intercept,  na.rm = TRUE),
    test_cal_slope      = mean(metrics$test_cal_slope,      na.rm = TRUE),
    test_cal_slopeDiff  = mean(metrics$test_cal_slopeDiff,  na.rm = TRUE)
  ), diff_tss <- train_tss - test_tss)
  
  out <- list(models = models, metrics = metrics, summary = summary)
  
  ## aggregate VI across folds: mean % and SD of %
  if (do_varimp) {
    mean_imp <- colMeans(vi_mat, na.rm = TRUE)
    sd_imp   <- apply(vi_mat, 2, stats::sd, na.rm = TRUE)
    
    out$varimp <- data.frame(
      variable   = names(mean_imp),
      importance = as.numeric(mean_imp),   
      sd         = as.numeric(sd_imp),     
      row.names  = NULL
    )
  }
  
  out
}

