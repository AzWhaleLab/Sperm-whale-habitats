              ##EFunctions for computation of spatial mean and SDs for environmental data from NETCDF files ##
##Developed by Mareike Duffing Romero. However, some codes for different sections were updated or used as a reference for building different functions and these were received from Dr. Elliott Hazen and Dr. Sergi Jorge-Perez and further updated and adapted to build these functions for the application of Mareike's PhD. Other contributing collaborator is M.Ines Pinhero for the code on extracting data from multiple netcdf files within a folder
## Workflow: ----
###Load input nc file 
###Break array by dimension (time, lat/long, depth, var). 

## Version Control ----
## Version 0.1
  ## Fxn_SpatialCalcMean: load files and compute the mean and standard deviation by time stamp (in this case day)
  ##Fxn_Extract variables: keeps similar workflow as previous workflow, but can keep building things further
##v0.2
  ##Fxn_ExtractVarCacl: change to input_folder from the previous document that we created as the output folder where netcdf are saved into.
        ##add variable of VariableShortN to include in the function to automate the function and create an input file inside the function to open up the netcdf file of interest.
##v0.3 added memory use and time stamps for extraction functions. Adds new computed variable into dataframe as cbind
##v.04 Updated both extraction functions to simplify the loop for the indices portion and remove the large data numeric item from global environment. Yet, added a code that pulls out the original name of the target dataframe and adds the new column of interest to the dataframe and pulls out the original dataframe with the new columns without needing to recall the output dataframe as the original dataframe. 
  ##fixed Fxn_SpatialCalcMean that was saving the wrong SD raster file and put the right code in. 
  ## added new function to create boxplots to compare differences in observed vs spatial mean and SD  
##v.05 added function for netcdf files whose variables may only have 3dimension rather than 4 (aka no depth) Fxn_Extract_EnviroRaw_3D
##v0.6 added new function for chlorophyll lag 1 month and 2 month, but same function can be applied to other variables as needed 
##v0.7: For all extraction functions add a slot of code that allows to pull out latitude, longitude, depth with different possible names (ie lat, latitude, long, Depth) and still pull out the variable of interest, rather than manually change it.
##v0.7 add a 3D version of lag function (without depth)
##v0.8 Fxn_Extract_EnviroRaw_Lag2Months to match the exact dates and not have any errors when extracting data 
##v0.9 Fxn_Extract_EnviroRaw_MultiNC: a fxn to extract your variable if you have a folder containing multiple nc files with your variable and file name has a "date" as part of the name.
##v1.0 updated Fxn_ExtractRaw to support 4D or 3D vars
##Load Packages ----
require(ncdf4)
require(lubridate)

##Let's build the functions ----
###Available functions: ----
    ###Fxn_Extract_EnviroRaw: extract variable from one single netcdf file that has a 4D (ie time/date, lat/long,depth, variable). 
    ###Fxn_Extract_EnviroRaw_3D: extract variable from one single data file that has 3D, meaning same as above minus depth. 
    ###Fxn_Extract_EnviroRaw_MultiNC: function to extract variable based on a folder that contains multiple netcdf files, specifically for cases for when you have a multiple netcdf whose file names contain a date eg. azores_netcdf_salt_1may2025.nc, etc. The function has been created in a way that it can read various dates naming conventions and pulls out the specific dates you need. It has been adapted based on a code that M.Ines Pinheiro provided.
    ###Fxn_Extract_EnviroRaw_Lag2Months: allows to create a 2 month lag to pull out a variable, this is mostly used for chlorophyll lags, but could be used for other variables 
    ###Fxn_Extract_EnviroCalc: this function still extracts variables, but based on previously built netcdf files whose means and SD has been spatially calculated at a 3x3 sqaure around a point. 
    ###Fxn_CalcSpatMeanSD: function to calculate mean and SD by a 3 x 3 grid within a netcdf file, it then creates a new netcdf file. 

###Items for the function ----
## Need input files and maintain the similar workflow as we had, but make sure the last variable is the one you choose. Have a If for Mean, SD, Original variables. 
##Input Folder
##variable_name,
##Underscore variable
## dataframe of interestt to merge to.
##Calculation=Mean, SD
##Calc_Variable = the one to enter for the IF function 
###Items for a function examples----
### Define the input NetCDF file and variable
# input_nc_file #<- "input_file.nc"  # Replace with your file path
# output_nc_file #<- "output_file.nc" # File to save output (can overwrite input if desired)
# variable_name #<- "your_variable"  # Replace with the variable name in your NetCDF ex "thetao", "so", "uo", "vo", etc
# output_folder # our output folder 
# VariableShortN # short variable name ie"SurfaceTemp" that can be used for saving output
# VariableLongN # longer variable "Surface Temperature"


##Create Functions for Extraction ----
###Function 1 Extraction Raw----

Fxn_Extract_EnviroRaw=function(input_file, ##input netcdf file
                               variable_name, ##original variable to pull out of the netcdf file
                               target_df, ## our dataframe containing the data that we'll be downloading enviro data to
                               new_var_name ##Variable name to be put in the data frame (ie. surface_temperature)
                               ) 
  {
  library(ncdf4); require(terra);require(lubridate); require(furrr)
  ##Console cleaning, saving memory, setting up times and setting up the field
  gc()
  
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  
  
  require(hms) #V0.3
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ##Begin the function
  print(paste0("Loading files and start extracting ", new_var_name)) #v0.2
  
  ## Pull out original name of target_df v0.4
  target_df_name <- deparse(substitute(target_df))
  
  ## Open NetCDF file
  nc <- nc_open(input_file)
  print(nc)
  
  ## detect variable dimensionality  #v1.0
  var_obj <- nc$var[[variable_name]]
  ndims_var <- var_obj$ndims
  
  # Helper to flexibly get variable by alternative names #v0.7
  get_nc_variable <- function(nc, possible_names, return_name = FALSE) {
    for (name in possible_names) {
      if (name %in% names(nc$var) || name %in% names(nc$dim)) {
        if (return_name) return(name)
        return(ncvar_get(nc, name))
      }
    }
    stop(paste("None of the variable names found:", paste(possible_names, collapse = ", ")))
  }
  
  ## Get latitude, longitude, depth, and time variables with potential variables naming conventions #v0.7
  latitudes  <- get_nc_variable(nc, c("latitude", "lat"))
  longitudes <- get_nc_variable(nc, c("longitude", "lon", "long"))
  depths     <- get_nc_variable(nc, c("depth", "Depth"))
  # print(paste0("depth:", depths))
  
  ## Get actual name of the time variable #v0.7
  time_varname <- get_nc_variable(nc, c("time", "Time"), return_name = TRUE);print(time_varname)
  
  ## Use it to get the data and attribute #v0.7
  time_var   <- ncvar_get(nc, time_varname)
  time_units <- ncatt_get(nc, time_varname, "units")$value

  ## Convert time variable to dates
  time_origin <- sub("seconds since ", "", time_units)
  time_dates <- as.POSIXct(time_var, origin = time_origin, tz = "UTC")
  
  time_slots_df <- as.data.frame(time_var)
  time_slots_df$time_dates <- time_dates
  
  ## Find the indices for target dates, latitudes, and longitudes
  time_dates <- as.Date(time_dates)
  find_indices <- function(x, value) which(x == value)
  time_indices <- sapply(target_df$date, find_indices, value = time_dates)
  lat_indices <- sapply(target_df$lat, function(lat) which.min(abs(latitudes - lat)))
  lon_indices <- sapply(target_df$lon, function(lon) which.min(abs(longitudes - lon)))
  
  ## Assume depth index is always 1 (modify as needed)
  depth_indices <- rep(1, length(lat_indices))
  
  ## Extract data updated #v0.4 
  ## Update extract value dimension of 4D or 3D all in one for the loop 
  extracted_values <- numeric(length(lat_indices))
  
  for (i in seq_along(lat_indices)) {
    
    lat  <- lat_indices[i]
    lon  <- lon_indices[i]
    time <- time_indices[i]
    
    if (ndims_var == 4) {
      
      depth <- depth_indices[i]
      
      extracted_values[i] <- ncvar_get(
        nc,
        variable_name,
        start = c(lon, lat, depth, time),
        count = c(1, 1, 1, 1)
      ) ##original v0.4
      
    } else if (ndims_var == 3) {
      
      extracted_values[i] <- ncvar_get(
        nc,
        variable_name,
        start = c(lon, lat, time),
        count = c(1, 1, 1)
      ) #new one v1.0 to support 3D
      
    } else {
      stop(paste("Unsupported number of dimensions:", ndims_var))
    }
  }
  
  
  ## Close NetCDF file
  nc_close(nc)
  
  ## Ensure multiple runs keep original columns and append new_var_name only if it doesn't already exist #v0.3 and 0.4
  if (!(new_var_name %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name))
  } else {
    target_df[[new_var_name]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  future::plan(future::sequential)  #return to single core
  
  print(paste0("Finish extracting ", new_var_name)) #v0.2
  
  require(hms) #V2
  EndTime=Sys.time()
  RunTimeAll=EndTime - StartTime;
  TimeString_Run=paste0("Start Time: ",
                        StartTime,
                        ";  End Time: ",
                        EndTime,
                        ";  Duration ",
                        as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_Run)
  require(pryr)
  loadMemory<-mem_used() 
  print(paste0("Memory used: ")); print(loadMemory) #v0.3
  gc()
  return(target_df)
} ##End of Fxn_Extract raw variables


###Function 1.2 Extraction Raw by multiple files----
Fxn_Extract_EnviroRaw_MultiNC <- function(
    dir,                ## directory with many nc files
    variable_name,      ## original nc variable (thetao, sob, etc.)
    target_df,          ## dataframe with Date, lat, lon
    new_var_name         ## name to add to dataframe
) {
  
  require(ncdf4)
  require(stringr)
  require(sp)
  require(hms)
  require(raster)
  
  gc()
  StartTime <- Sys.time()
  print(paste0("Start Time: ", StartTime))
  print(paste0("Extracting ", new_var_name))
  
  ## Keep original df name
  target_df_name <- deparse(substitute(target_df))
  
  ## --- index NetCDF files by date ---
  nc_files <- list.files(dir, pattern = "\\.nc$", full.names = TRUE)
  
  date_str <- str_extract(
    basename(nc_files),
    "\\d{8}-\\d{8}|\\d{8}|\\d{4}-\\d{2}-\\d{2}|\\d{2}-\\d{2}-\\d{4}|\\d{4}_\\d{2}_\\d{2}"
  )
  
  date_str <- ifelse(
    grepl("\\d{8}-\\d{8}", date_str),
    substr(date_str, 1, 8),
    date_str
  )
  
  file_dates <- as.Date(NA)
  file_dates <- ifelse(grepl("^\\d{8}$", date_str), as.Date(date_str, "%Y%m%d"),
                       ifelse(grepl("^\\d{4}-\\d{2}-\\d{2}$", date_str), as.Date(date_str, "%Y-%m-%d"),
                              ifelse(grepl("^\\d{2}-\\d{2}-\\d{4}$", date_str), as.Date(date_str, "%d-%m-%Y"),
                                     ifelse(grepl("^\\d{4}_\\d{2}_\\d{2}$", date_str), as.Date(date_str, "%Y_%m_%d"),
                                            NA))))
  
  nc_index <- data.frame(Date = file_dates, file = nc_files, stringsAsFactors = FALSE)
  
  ## --- prepare output ---
  out_values <- rep(NA_real_, nrow(target_df))
  
  target_df$Date <- as.Date(target_df$date)
  unique_dates <- sort(unique(target_df$Date))
  
  ## --- loop per date ---
  for (d in unique_dates) {
    
    message("Processing date: ", as.character(d))
    
    rows <- which(target_df$Date == d)
    hh4  <- target_df[rows, ]
    
    coordinates(hh4) <- ~lon + lat
    
    ncname <- nc_index$file[nc_index$Date == d]
    
    if (length(ncname) == 0) {
      warning("No NetCDF file found for date: ", d)
      next
    }
    
    ## read raster directly (no nc_open needed here)
    r <- raster(ncname[1], varname = variable_name)
    
    out_values[rows] <- raster::extract(r, hh4)
  }
  
  ## --- append output ---
  if (!(new_var_name %in% colnames(target_df))) {
    target_df[[new_var_name]] <- out_values
  } else {
    target_df[[new_var_name]] <- out_values
  }
  
  ## MDR test to remove Date 
  target_df<-subset(target_df, select=-c(Date))

  ## --- return df to global env (your convention) ---
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  EndTime <- Sys.time()
  RunTimeAll <- EndTime - StartTime
  
  print(paste0(
    "End Time: ", EndTime,
    "; Duration ", as_hms(round(RunTimeAll, 2))
  ))
  
  gc()
  return(target_df)
}


###Function 1b Extraction Raw 3D----

Fxn_Extract_EnviroRaw_3D=function(input_file, ##input netcdf file
                               variable_name, ##original variable to pull out of the netcdf file
                               target_df, ## our dataframe containing the data that we'll be downloading enviro data to
                               new_var_name ##Variable name to be put in the data frame (ie. surface_temperature)
                               ) 
  {
  library(ncdf4); require(terra);require(lubridate); require(furrr)
  ##Console cleaning, saving memory, setting up times and setting up the field
  gc()
  
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  
  
  require(hms) #V0.3
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ##Begin the function
  print(paste0("Loading files and start extracting ", new_var_name)) #v0.2
  
  ## Pull out original name of target_df v0.4
  target_df_name <- deparse(substitute(target_df))
  
  ## Open NetCDF file
  nc <- nc_open(input_file)

  ##Helper to flexibly get variable by alternative names #v0.7
  get_nc_variable <- function(nc, possible_names, return_name = FALSE) {
    for (name in possible_names) {
      if (name %in% names(nc$var) || name %in% names(nc$dim)) {
        if (return_name) return(name)
        return(ncvar_get(nc, name))
      }
    }
    stop(paste("None of the variable names found:", paste(possible_names, collapse = ", ")))
  }
  
  ## Get latitude, longitude, depth, and time variables with potential variables naming conventions #v0.7
  latitudes  <- get_nc_variable(nc, c("latitude", "lat"))
  longitudes <- get_nc_variable(nc, c("longitude", "lon", "long"))
  # depths     <- get_nc_variable(nc, c("depth", "Depth"))
  # print(paste0("depth:", depths))
  
  ## Get actual name of the time variable #v0.7
  time_varname <- get_nc_variable(nc, c("time", "Time"), return_name = TRUE);print(time_varname)
  
  ## Use it to get the data and attribute #v0.7
  time_var   <- ncvar_get(nc, time_varname)
  time_units <- ncatt_get(nc, time_varname, "units")$value
  
  ## Convert time variable to dates
  time_origin <- sub("seconds since ", "", time_units)
  time_dates <- as.POSIXct(time_var, origin = time_origin, tz = "UTC")
  
  time_slots_df <- as.data.frame(time_var)
  time_slots_df$time_dates <- time_dates
  
  ## Find the indices for target dates, latitudes, and longitudes
  time_dates <- as.Date(time_dates)
  find_indices <- function(x, value) which(x == value)
  time_indices <- sapply(target_df$date, find_indices, value = time_dates)
  lat_indices <- sapply(target_df$lat, function(lat) which.min(abs(latitudes - lat)))
  lon_indices <- sapply(target_df$lon, function(lon) which.min(abs(longitudes - lon)))
  
  ## Assume depth index is always 1 (modify as needed)
  depth_indices <- rep(1, length(lat_indices))
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    # depth <- depth_indices[i]
    time <- time_indices[i]
    # extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, depth, time), count = c(1, 1, 1, 1))
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, time), count = c(1, 1, 1))
  }
  
  ## Close NetCDF file
  nc_close(nc)
  
  ## Ensure multiple runs keep original columns and append new_var_name only if it doesn't already exist #v0.3 and 0.4
  if (!(new_var_name %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name))
  } else {
    target_df[[new_var_name]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  future::plan(future::sequential)  #return to single core
  
  print(paste0("Finish extracting ", new_var_name)) #v0.2
  
  require(hms) #V2
  EndTime=Sys.time()
  RunTimeAll=EndTime - StartTime;
  TimeString_Run=paste0("Start Time: ",
                        StartTime,
                        ";  End Time: ",
                        EndTime,
                        ";  Duration ",
                        as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_Run)
  require(pryr)
  loadMemory<-mem_used() 
  print(paste0("Memory used: ")); print(loadMemory) #v0.3
  gc()
  return(target_df)
} ##End of Fxn_Extract raw variables



###Function 1c Extraction Raw with time lag----

Fxn_Extract_EnviroRaw_Lag2Months=function(input_file, ##input netcdf file
                               variable_name, ##original variable to pull out of the netcdf file
                               target_df, ## our dataframe containing the data that we'll be downloading enviro data to
                               new_var_name ##Variable name to be put in the data frame (ie. surface_temperature)
) 
{
  ##Note make sure you create a lag time before hand
  library(ncdf4); require(terra);require(lubridate); require(furrr)
  ##Console cleaning, saving memory, setting up times and setting up the field
  gc()
  
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  
  
  require(hms) #V0.3
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ##Begin the function
  print(paste0("Loading files and start extracting ", new_var_name)) #v0.2
  
  ## Pull out original name of target_df v0.4
  target_df_name <- deparse(substitute(target_df))
  
  ## Open NetCDF file
  nc <- nc_open(input_file)
  print(nc)
  ##Helper to flexibly get variable by alternative names #v0.7
  get_nc_variable <- function(nc, possible_names, return_name = FALSE) {
    for (name in possible_names) {
      if (name %in% names(nc$var) || name %in% names(nc$dim)) {
        if (return_name) return(name)
        return(ncvar_get(nc, name))
      }
    }
    stop(paste("None of the variable names found:", paste(possible_names, collapse = ", ")))
  }
  
  ## Get latitude, longitude, depth, and time variables with potential variables naming conventions #v0.7
  latitudes  <- get_nc_variable(nc, c("latitude", "lat"))
  longitudes <- get_nc_variable(nc, c("longitude", "lon", "long"))
  depths     <- get_nc_variable(nc, c("depth", "Depth"))
  # print(paste0("depth:", depths))
  
  ## Get actual name of the time variable #v0.7
  time_varname <- get_nc_variable(nc, c("time", "Time"), return_name = TRUE);print(time_varname)
  
  ## Use it to get the data and attribute #v0.7
  time_var   <- ncvar_get(nc, time_varname)
  time_units <- ncatt_get(nc, time_varname, "units")$value
  
  ## Convert time variable to dates
  time_origin <- sub("seconds since ", "", time_units)
  time_dates <- as.POSIXct(time_var, origin = time_origin, tz = "UTC")
  time_dates <- as.Date(time_dates) ##moved up

  time_slots_df <- as.data.frame(time_var)
  time_slots_df$time_dates <- time_dates
  
  ##Create lag months in target dataframe #v0.6
  print("Creating lag months in target df")
  target_df <- target_df %>%
    mutate(
      date_minus_1m = date %m-% months(1),
      date_minus_2m = date %m-% months(2)
    )
  
  ##Lag 1 v06 ----
  print("Extracting variables for lag minus 1 months")
  
  ## Find the indices for target dates, latitudes, and longitudes within df
  find_indices <- function(x, value) which(x == value)
  # time_indices <- sapply(target_df$date_minus_1m, find_indices, value = time_dates) ##For lag 1  #v0.6 subset
  time_indices <- match(as.Date(target_df$date_minus_1m), time_dates) ##for lag1 #v0.8 updated to match dates
  lat_indices <- sapply(target_df$lat, function(lat) which.min(abs(latitudes - lat)))
  lon_indices <- sapply(target_df$lon, function(lon) which.min(abs(longitudes - lon)))
  
  ## Assume depth index is always 1 (modify as needed)
  depth_indices <- rep(1, length(lat_indices))
  
  ##v0.8 add empty space if needed
  extracted_values <- rep(NA_real_, length(lat_indices))
  valid_idx <- which(!is.na(time_indices))
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    depth <- depth_indices[i]
    time <- time_indices[i]
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, depth, time), count = c(1, 1, 1, 1))
  }
  
  
  ## Ensure multiple runs keep original columns and append new_var_name only if it doesn't already exist #v0.3 and 0.4
  ##For lag 1 #v0.6
  new_var_name_lag1=paste0(new_var_name, "_lag1Mo")
  if (!(new_var_name_lag1 %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name_lag1))
  } else {
    target_df[[new_var_name_lag1]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  ##For Lag 2 #v0.6 ----
  print("Extracting variables for lag minus 2 months")
  # time_indices <- sapply(target_df$date_minus_2m, find_indices, value = time_dates) 
  time_indices <- match(as.Date(target_df$date_minus_2m), time_dates) ##v0.8 updated to match propery
  
  ##v0.8 add empty space if needed
  extracted_values <- rep(NA_real_, length(lat_indices))
  valid_idx <- which(!is.na(time_indices))
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    depth <- depth_indices[i]
    time <- time_indices[i]
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, depth, time), count = c(1, 1, 1, 1))
  }
  
  ##For lag 1
  new_var_name_lag2=paste0(new_var_name, "_lag2Mo")
  if (!(new_var_name_lag2 %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name_lag2))
  } else {
    target_df[[new_var_name_lag2]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  ## Close NetCDF file
  nc_close(nc)
  
  future::plan(future::sequential)  #return to single core
  
  print(paste0("Finish extracting 2 months of lag for: ", new_var_name)) #v0.2
  
  ##Remove the date months
  # target_df<-target_df %>% select(!c(date_minus_1m,date_minus_2m))
  target_df<-subset(target_df, select=-c(date_minus_1m, date_minus_2m))
  require(hms) #V2
  EndTime=Sys.time()
  RunTimeAll=EndTime - StartTime;
  TimeString_Run=paste0("Start Time: ",
                        StartTime,
                        ";  End Time: ",
                        EndTime,
                        ";  Duration ",
                        as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_Run)
  require(pryr)
  loadMemory<-mem_used() 
  print(paste0("Memory used: ")); print(loadMemory) #v0.3
  gc()
  return(target_df)
} ##End of Fxn_Extract raw variables with time lag 



###Function 1d Extraction Raw with time lag for 3D----

Fxn_Extract_EnviroRaw_Lag2Months_3D=function(input_file, ##input netcdf file
                               variable_name, ##original variable to pull out of the netcdf file
                               target_df, ## our dataframe containing the data that we'll be downloading enviro data to
                               new_var_name ##Variable name to be put in the data frame (ie. surface_temperature)
) 
{
  ##Note make sure you create a lag time before hand
  library(ncdf4); require(terra);require(lubridate); require(furrr)
  ##Console cleaning, saving memory, setting up times and setting up the field
  gc()
  
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  
  
  require(hms) #V0.3
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ##Begin the function
  print(paste0("Loading files and start extracting ", new_var_name)) #v0.2
  
  ## Pull out original name of target_df v0.4
  target_df_name <- deparse(substitute(target_df))
  
  ## Open NetCDF file
  nc <- nc_open(input_file)
  print(nc)
  ##Helper to flexibly get variable by alternative names #v0.7
  get_nc_variable <- function(nc, possible_names, return_name = FALSE) {
    for (name in possible_names) {
      if (name %in% names(nc$var) || name %in% names(nc$dim)) {
        if (return_name) return(name)
        return(ncvar_get(nc, name))
      }
    }
    stop(paste("None of the variable names found:", paste(possible_names, collapse = ", ")))
  }
  
  ## Get latitude, longitude, depth, and time variables with potential variables naming conventions #v0.7
  latitudes  <- get_nc_variable(nc, c("latitude", "lat"))
  longitudes <- get_nc_variable(nc, c("longitude", "lon", "long"))
  # #depths     <- get_nc_variable(nc, c("depth", "Depth")) #not needed
  # print(paste0("depth:", depths))
  
  ## Get actual name of the time variable #v0.7
  time_varname <- get_nc_variable(nc, c("time", "Time"), return_name = TRUE);print(time_varname)
  
  ## Use it to get the data and attribute #v0.7
  time_var   <- ncvar_get(nc, time_varname)
  time_units <- ncatt_get(nc, time_varname, "units")$value
  
  ## Convert time variable to dates
  time_origin <- sub("seconds since ", "", time_units)
  time_dates <- as.POSIXct(time_var, origin = time_origin, tz = "UTC")
  
  time_slots_df <- as.data.frame(time_var)
  time_slots_df$time_dates <- time_dates
  
  ##Create lag months in target dataframe #v0.6
  print("Creating lag months in target df")
  target_df <- target_df %>%
    mutate(
      date_minus_1m = date %m-% months(1),
      date_minus_2m = date %m-% months(2)
    )
  
  print("Extracting variables for lag minus 1 months")
  
  ## Find the indices for target dates, latitudes, and longitudes
  time_dates <- as.Date(time_dates)
  find_indices <- function(x, value) which(x == value)
  
  time_indices <- sapply(target_df$date_minus_1m, find_indices, value = time_dates) ##For lag 1  #v0.6 subset
  lat_indices <- sapply(target_df$lat, function(lat) which.min(abs(latitudes - lat)))
  lon_indices <- sapply(target_df$lon, function(lon) which.min(abs(longitudes - lon)))
  
  ## Assume depth index is always 1 (modify as needed)
  depth_indices <- rep(1, length(lat_indices))
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    # #depth <- depth_indices[i] #not needed
    time <- time_indices[i]
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, time), count = c(1, 1, 1))
  } ##end of fxn extract lag 3d
  
  
  ## Ensure multiple runs keep original columns and append new_var_name only if it doesn't already exist #v0.3 and 0.4
  ##For lag 1 #v0.6
  new_var_name_lag1=paste0(new_var_name, "_lag1Mo")
  if (!(new_var_name_lag1 %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name_lag1))
  } else {
    target_df[[new_var_name_lag1]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  ##For Lag 2 #v0.6
  print("Extracting variables for lag minus 2 months")
  time_indices <- sapply(target_df$date_minus_2m, find_indices, value = time_dates) 
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    # #depth <- depth_indices[i] #not needed
    time <- time_indices[i]
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, time), count = c(1, 1, 1))
  }
  
  ##For lag 1
  new_var_name_lag2=paste0(new_var_name, "_lag2Mo")
  if (!(new_var_name_lag2 %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), new_var_name_lag2))
  } else {
    target_df[[new_var_name_lag2]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment ##v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  ## Close NetCDF file
  nc_close(nc)
  
  future::plan(future::sequential)  #return to single core
  
  print(paste0("Finish extracting 2 months of lag for: ", new_var_name)) #v0.2
  
  ##Remove the date months
  # target_df<-target_df %>% select(!c(date_minus_1m,date_minus_2m))
  target_df<-subset(target_df, select=-c(date_minus_1m, date_minus_2m))
  require(hms) #V2
  EndTime=Sys.time()
  RunTimeAll=EndTime - StartTime;
  TimeString_Run=paste0("Start Time: ",
                        StartTime,
                        ";  End Time: ",
                        EndTime,
                        ";  Duration ",
                        as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_Run)
  require(pryr)
  loadMemory<-mem_used() 
  print(paste0("Memory used: ")); print(loadMemory) #v0.3
  gc()
  return(target_df)
} ##End of Fxn_Extract raw variables with time lag 



###Function 2: Extract Calculated Variables ----
##
Fxn_Extract_EnviroCalc = function(input_folder, ##input folder that has the files v0.2 changed from input file to folder
                                  variable_name, ##original variable to pull out of the netcdf file
                                  target_df, ## our dataframe containing the data that we'll be downloading enviro data to
                                  new_var_name, ##Variable name to be put in the data frame (ie. surface_temperature, 
                                  calc, ##what calculating was applied to it: Mean or SD?,
                                  VariableShortN ## short variable 
                                  )

{
  gc()
  require(ncdf4); require(terra); require(furrr) ##added furr to help run things in parallel
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  
  
  require(hms) #V0.3
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ## Pull out original name of target_df #v0.4
  target_df_name <- deparse(substitute(target_df))
  
  ## Open NetCDF file
  print(paste0("Loading files and start extracting calculated ", calc, " of ", VariableShortN)) #v0.2
  input_file=paste0(input_folder, "/", VariableShortN, "_Calc_",calc,".nc") ##0.2
  nc <- nc_open(input_file)
  
  ##Helper to flexibly get variable by alternative names #v0.7
  get_nc_variable <- function(nc, possible_names, return_name = FALSE) {
    for (name in possible_names) {
      if (name %in% names(nc$var) || name %in% names(nc$dim)) {
        if (return_name) return(name)
        return(ncvar_get(nc, name))
      }
    }
    stop(paste("None of the variable names found:", paste(possible_names, collapse = ", ")))
  }
  
  ## Get latitude, longitude, depth, and time variables with potential variables naming conventions #v0.7
  latitudes  <- get_nc_variable(nc, c("latitude", "lat"))
  longitudes <- get_nc_variable(nc, c("longitude", "lon", "long"))
  depths     <- get_nc_variable(nc, c("depth", "Depth"))
  # print(paste0("depth:", depths))
  
  ## Get actual name of the time variable #v0.7
  time_varname <- get_nc_variable(nc, c("time", "Time"), return_name = TRUE);print(time_varname)
  
  ## Use it to get the data and attribute #v0.7
  time_var   <- ncvar_get(nc, time_varname)
  time_units <- ncatt_get(nc, time_varname, "units")$value
  
  
  # Convert time variable to dates
  time_origin <- sub("seconds since ", "", time_units)
  time_dates <- as.POSIXct(time_var, origin = time_origin, tz = "UTC")
  
  time_slots_df <- as.data.frame(time_var)
  time_slots_df$time_dates <- time_dates
  
  ## Find the indices for target dates, latitudes, and longitudes
  time_dates <- as.Date(time_dates)
  find_indices <- function(x, value) which(x == value)
  time_indices <- sapply(target_df$date, find_indices, value = time_dates)
  lat_indices <- sapply(target_df$lat, function(lat) which.min(abs(latitudes - lat)))
  lon_indices <- sapply(target_df$lon, function(lon) which.min(abs(longitudes - lon)))
  
  ## Assume depth index is always 1 (modify as needed)
  depth_indices <- rep(1, length(lat_indices))
  
  # ## Extract data
  # for (i in 1:length(lat_indices)) {
  #   lat <- lat_indices[i]
  #   lon <- lon_indices[i]
  #   # depth <- depth_indices[i]
  #   time <- time_indices[i]
  #   data <- ncvar_get(nc, variable_name, start = c(lon, lat, time), count = c(1, 1, 1))
  #   hh5 <- as.data.frame(data)
  #   if (i == 1) {
  #     hh6 <- hh5
  #   } else {
  #     hh7 <- rbind(hh6, hh5)
  #     hh6 <- hh7
  #   }
  # }
  
  ## Extract data updated #v0.4 
  extracted_values <- numeric(length(lat_indices))
  for (i in seq_along(lat_indices)) {
    lat <- lat_indices[i]
    lon <- lon_indices[i]
    time <- time_indices[i]
    extracted_values[i] <- ncvar_get(nc, variable_name, start = c(lon, lat, time), count = c(1, 1, 1))
  }
  
  
  ## Close NetCDF file
  nc_close(nc)
  
  ## Rename extracted variable and add to dataframe. But first rename the variable to include the mean or sd epiphet
  CalcVar_name=paste0(new_var_name,  "_", calc)

  # ## Ensure multiple runs keep original columns and append CalcVar_name
  # if (!(CalcVar_name %in% colnames(target_df))) {
  #   target_df <- cbind(target_df, setNames(as.data.frame(hh7$data), CalcVar_name))
  # } else {
  #   target_df[[CalcVar_name]] <- hh7$data
  # }
  
  # Ensure multiple runs keep original columns and append new_var_name only if it doesn't already exist
  if (!(CalcVar_name %in% colnames(target_df))) {
    target_df <- cbind(target_df, setNames(as.data.frame(extracted_values), CalcVar_name))
  } else {
    target_df[[CalcVar_name]] <- extracted_values
  }
  
  ## Assign updated dataframe back to its original name in the global environment v0.4
  assign(target_df_name, target_df, envir = .GlobalEnv)
  
  future::plan(future::sequential)  #return to single core
  
  print(paste0("Finish extracting calculated ", calc, " of ", VariableShortN)) #v0.2
  
  require(hms) #V2
  EndTime=Sys.time()
  RunTimeAll=EndTime - StartTime;
  TimeString_Run=paste0("Start Time: ",
                        StartTime,
                        ";  End Time: ",
                        EndTime,
                        ";  Duration ",
                        as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_Run)
  require(pryr)
  loadMemory<-mem_used() 
  print(paste0("Memory used: ")); print(loadMemory) #v0.3
  gc()
  
  return(target_df)
} #end fxn_for extracting calculated variables


##Function for Computing Spatial means and SDs -----
Fxn_CalcSpatMeanSD=function(input_nc_file, ##input folder
                            output_folder, ##Output folder path
                            variable_name, ##variable name to be pulled out of netcdf file
                            VariableShortN, ## short variable 
                            VariableLongN
){
  require(terra); require(ncdf4);require(lubridate)
  require(furrr) ##added furr to help run things in parallel and reduce the load as needed
  
  ##These 2 lines of code are for setting up your computer and make sure you can run the code.
  progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
  future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel
  
  require(hms)
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ## Load raster from NetCDF
  print(paste0("Loading input raster for:", VariableLongN))
  var_data <- rast(input_nc_file, subds = variable_name)  
  
  ## Get number of time steps
  time_steps <- nlyr(var_data)
  
  ## Apply focal function for mean and standard deviation at each time step
  mean_list <- list()
  sd_list <- list()
  
  print("Calculating Means & SD")
  for (i in 1:time_steps) {
    mean_list[[i]] <- focal(var_data[[i]], w = matrix(1, 3, 3), fun = mean, na.rm = TRUE)
    sd_list[[i]] <- focal(var_data[[i]], w = matrix(1, 3, 3), fun = sd, na.rm = TRUE)
  }
  
  ## Convert lists to rasters
  mean_stack <- rast(mean_list)
  sd_stack <- rast(sd_list)
  
  # Rename the varnames in the stacks
  mean_varname <- paste0("mean_", variable_name)
  sd_varname <- paste0("sd_", variable_name)
  
  varnames(mean_stack) <- rep(mean_varname, nlyr(mean_stack))
  varnames(sd_stack) <- rep(sd_varname, nlyr(sd_stack))
  
  ## Merge both into a single raster
  # focal_raster <- c(mean_stack, sd_stack) #hold off for a while v0.2
  
  ## Assign a simple CRS (EPSG:4326) or remove CRS
  # crs(focal_raster) <- "EPSG:4326"  # Standard WGS84
  crs(mean_stack) <- "EPSG:4326"  # Standard WGS84
  crs(sd_stack) <- "EPSG:4326"  # Standard WGS84
  
  ##Extract original units
  nc <- nc_open(input_nc_file)
  original_units <- ncatt_get(nc, variable_name, "units")$value
  nc_close(nc)
  
  ##Save output
  print("Saving Output")
  var_calc_name=paste0("mean_", variable_name) ##Maybe
  # OutputFile_Mean=paste0(NC_2023_WD, "/CalcProducts_2023-24/", VariableShortN, "_CalcMean.nc")
  OutputFile_Mean=paste0(output_folder,"/",VariableShortN, "_Calc_Mean.nc")
  writeCDF(mean_stack,
           filename = OutputFile_Mean, 
           overwrite=TRUE, 
           varname=variable_name
           # varname=var_calc_name
  )
  
  OutputFile_SD=paste0(output_folder,"/",VariableShortN, "_Calc_SD.nc")
  writeCDF(sd_stack, #v0.4 updated to the proper item of sd 
           filename = OutputFile_SD, 
           overwrite=TRUE, 
           varname=variable_name
           # varname=var_calc_name
  )
  
  require(pryr)
  loadMemory<-mem_used()
  print(paste0("memory used:", loadMemory))
  
  Results<-list(mean_stack, sd_stack) #, focal_raster ##hold off for a while v0.2
  print(paste0("Complete calculation of mean and sd for: ", VariableLongN))
  
  EndSave=Sys.time(); 
  
  RunTimeAll=EndSave - StartTime;
  TimeString_Run2=paste0("Calculation Start Time: ", 
                         StartTime,
                         ";  End Time after saving files: ",
                         EndSave,
                         ";  Duration ",
                         as_hms(round(RunTimeAll,2)))
  future::plan(future::sequential)  #return to single core
  gc()
  return(Results); 
}


##Function for plotting boxplots that compare  observed and calculated variables ----
Fxn_CompObsCalc_Box=function(df, #dataframe
                             var_small, #variable in underscore
                             VarName, #clean variable
                             OutputFolder
                             
){
  require(tidyr); require(tidyselect); require(dplyr)
  ##Start
  df$type<-as.factor(df$type)
  var.sel<-c(var_small)
  var.dyn_Mean<-c(paste0(var_small,"_Mean"))
  var.dyn_SD<-c(paste0(var_small,"_SD"))
  
  ##Select variables
  df_sel<-df %>% select(all_of(c("type","month",var.sel)))
  df_sel_Mean<-df %>% select(all_of(c("type","month",var.dyn_Mean)))
  df_sel_SD<-df %>% select(all_of(c("type","month",var.dyn_SD)))
  
  ##rename variables had as var_small
  colnames(df_sel)[3] <- "variable"
  colnames(df_sel_Mean)[3] <- "variable"
  colnames(df_sel_SD)[3] <- "variable"
  
  ##Add variable of interest
  df_sel$calc<-"Observed"
  df_sel_Mean$calc<-"Spatial Mean"
  df_sel_SD$calc<-"Spatial SD"
  
  df_sel$VarName<-var_small
  df_sel_Mean$VarName<-var_small
  df_sel_SD$VarName<-var_small
  
  ##Create new dataframe 
  var_df_long<-rbind(df_sel, df_sel_Mean, df_sel_Mean)
  var_df_ObsMean<-rbind(df_sel, df_sel_Mean)
  
  require(ggplot2); require(cowplot)
  
  box<-ggplot()+
    geom_boxplot(data = var_df_ObsMean, aes(x=calc,y=variable, fill = type), width=0.5)+
    theme_classic()+
    theme(legend.position = "bottom")
  box
  
  box1<-box+
    theme(legend.position = "none",
          axis.title = element_blank())
  box1
  
  box2<-ggplot()+
    geom_boxplot(data = df_sel_SD, aes(x=calc,y=variable, fill = type), width=0.5)+
    theme_classic()+
    theme(legend.position = "none",
          axis.title = element_blank())
  box2
  plot<-plot_grid(box1, box2, ncol = 2, 
                  rel_widths = c(1, 0.8)
  )
  plot
  
  # Extract the legend using ggplotGrob()
  g <- ggplotGrob(box)
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  
  ##Get axes and start building
  y_axis <- ggdraw() + draw_label(paste0(VarName), angle = 90, size = 9) #, fontface = "bold"
  
  final_plot_with_axes <- plot_grid(y_axis, plot, ncol = 2, rel_widths = c(0.1, 1))
  
  final_plot_with_legend <- plot_grid(final_plot_with_axes, legend, ncol = 1, rel_heights = c(1,0.1)) #0.1,NULL,
  
  title <- ggdraw() + draw_label(paste0("Differences between observed and \nspatial calculation of ", VarName), 
                                 fontface = "bold", size = 10.5)
  
  final_arranged_plot <- plot_grid(title,NULL,  final_plot_with_legend, ncol = 1, rel_heights = c(0.1, 0.02,1))
  
  # Print the final plot
  # print(final_arranged_plot)
  
  ggsave(final_arranged_plot, filename = paste0(OutputFolder, "/", var_small, "_DifObsCalc_boxplot.tiff"), width = 5, height = 2.8, units = "in", dpi = 200, bg="white", compression="lzw") #width = 11, height = 7.5,
  return(final_arranged_plot)
} ##End of Fxn_CompObsCalc_Box #v0.4 after the fact new function 
