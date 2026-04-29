                  ### Part 2b: download environmental data for Prediction code

## Script for downloading environmental data for creating prediction code.
## See notes from Part 2 downloading environmental data. 
##Here we follow similar workflow as in that script, however, we create a grid dataframe containing the temporal scale we want to help extract the specific data from the netcdf. Im this case we used the first Day of each month. 
## See Methods in publication to understand better the different sources of variables used for this.
## Need to upload functions below to help extract variables
## Most of the script follows the same workflow as Hazen et al (2017) and Pérez-Jorge et al 2021. However, some items were turned into functions to help facilitate the extraction process, and many of the codes were updated to reflect updates to packages or transitions into new spatial packages from retired spatial packages. 
##In our workflow given that we have 2018, 2023, 2024, the workflow for extracting dynamic variables was split into 2018 data and 2023-2024 data. Then the data was combined to extract static variables. 
##Authors: Mareike D. Duffing Romero, Elliott Hazen, Sergi Pérez-Jorge. 
##There are notes and clarification throughout the document for seperate sections as needed, mostly defined as "##". Similarly, you can close/open different sections as needed using the navigation pane.

##Load files ----
library(ncdf4)
library(lubridate)

rm(list = ls())

##Setting working directories 
SDM_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD) ## SDM Directory
Results<-"~//Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/PredictionData" 
Sergi_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi" ##for spatial files
Version="v1.0_beta.R"

source(paste0(SDM_HD, "/Functions/Function_SpatialMeansSD_NetCDF_", Version))

##Create grid for data extraction ----
library(sf); library(sp); library(lubridate)
grid_base_pre <- st_read(paste0(Sergi_HD,"/Grid/Centroids/Grid3/Grid3_Centroids.shp")) 
grid_base<-grid_base_pre

# grid_base<-st_transform(grid_base, grid_base_pre)
coordinates<-st_coordinates(grid_base)
centroids<-as.data.frame(coordinates)
colnames(centroids)[1] <- "lon"
colnames(centroids)[2] <- "lat"

##Create target dates
Feb_2018<-rep("2018-02-01", each = 1, len = nrow(centroids)) 
March_2018<-rep("2018-03-01", each = 1, len = nrow(centroids)) 
April_2018<-rep("2018-04-01", each = 1, len = nrow(centroids)) 

May_2018<-rep("2018-05-01", each = 1, len = nrow(centroids)) 
June_2018<-rep("2018-06-01", each = 1, len = nrow(centroids)) 
July_2018<-rep("2018-07-01", each = 1, len = nrow(centroids)) 
August_2018<-rep("2018-08-01", each = 1, len = nrow(centroids)) 
September_2018<-rep("2018-09-01", each = 1, len = nrow(centroids)) 
October_2018<-rep("2018-10-01", each = 1, len = nrow(centroids)) 

date2018<-rbind(Feb_2018, March_2018, April_2018,#v3
                May_2018,June_2018, July_2018, August_2018, September_2018, ##old
                October_2018 #v3
) #n=5-->n=9 for v3

##v3 (Jan 2026): Added Feb-June 2023 (n=5 more), total 11 months (n=11)
Feb_2023<-rep("2023-02-01", each = 1, len = nrow(centroids)) 
March_2023<-rep("2023-03-01", each = 1, len = nrow(centroids)) 
April_2023<-rep("2023-04-01", each = 1, len = nrow(centroids)) 
May_2023<-rep("2023-05-01", each = 1, len = nrow(centroids)) 
June_2023<-rep("2023-06-01", each = 1, len = nrow(centroids)) 

July_2023<-rep("2023-07-01", each = 1, len = nrow(centroids)) 
August_2023<-rep("2023-08-01", each = 1, len = nrow(centroids)) 
September_2023<-rep("2023-09-01", each = 1, len = nrow(centroids)) 
October_2023<-rep("2023-10-01", each = 1, len = nrow(centroids))
November_2023<-rep("2023-11-01", each = 1, len = nrow(centroids))
December_2023<-rep("2023-12-01", each = 1, len = nrow(centroids))

January_2024<-rep("2024-01-01", each = 1, len = nrow(centroids)) 
February_2024<-rep("2024-02-01", each = 1, len = nrow(centroids)) 
March_2024<-rep("2024-03-01", each = 1, len = nrow(centroids)) 
April_2024<-rep("2024-04-01", each = 1, len = nrow(centroids)) 
May_2024<-rep("2024-05-01", each = 1, len = nrow(centroids)) 
June_2024<-rep("2024-06-01", each = 1, len = nrow(centroids)) 
July_2024<-rep("2024-07-01", each = 1, len = nrow(centroids)) 
August_2024<-rep("2024-08-01", each = 1, len = nrow(centroids)) #n=14 (Aug2024 data); 
September_2024<-rep("2024-09-01", each = 1, len = nrow(centroids))
October_2024<-rep("2024-10-01", each = 1, len = nrow(centroids))  #n=16

#v3: new date 202324 has extra months for 2023 total is 21 months
date202324<-rbind(Feb_2023, March_2023, April_2023, May_2023, June_2023, #v3
                  July_2023,August_2023,September_2023,October_2023,November_2023,December_2023, January_2024,February_2024,March_2024,April_2024,May_2024,June_2024, July_2024, August_2024, September_2024, October_2024)

date<-rbind(date2018, date202324)
date<-as.Date(date)

date_final <- date[order(date)]

##Merge, need 19 centroids (Aug2024 data); for October 2024 need 21 ##old v2
##Merge, need 9 more, total 30 centroids
centroids_final<-rbind(centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids, 
                       centroids,
                       centroids,
                       centroids,
                       centroids,
                       centroids, 
                       centroids, 
                       centroids)

centroids_new<-rbind(centroids,
                     centroids,
                     centroids,
                     centroids,
                     centroids,
                     centroids,
                     centroids,
                     centroids,
                     centroids)
centroids_final<-rbind(centroids_final, centroids_new)

pred_data<-cbind(centroids_final,date_final)
pred_data$date<-pred_data$date_final
pred_data$year<-year(pred_data$date_final)

## Download environmental data for Predictions----
### Data 2023-2024 ----
library(dplyr)
all_tracks_a<-filter(all_tracks, TagYear %in% c("2023", "2024"))
# all_tracks_a<-filter(all_tracks_a, sen_rep %in% c("0","1")) ##if you want to test a subset of data

require(hms) #V2
StartTime=Sys.time() ; StartTime
print(paste0("Start Time:", StartTime))

StartTime_Dyn23=Sys.time() ; StartTime
print(paste0("Start Time Dynamic 2023-24:", StartTime))

##Here you can update to reflect your directory were you will have your NC files 
NC_2023_WD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Environmental_variables/Products_2023_2024" ##New data Nov2024

NC_2023_WD_b<-"~/Projects/SpermWhale/3_SpeciesDistMod/Environmental_variables/Products_2023_2024_Aug2024" ##old data Aug2024

#### potential temperature (sea surface temperature) ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"), 
  target_df=all_tracks_a,
  variable_name="thetao",
  new_var_name = "surface_temperature"
)

#### potential temperature (mid temperature) ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_mid.nc"), 
  target_df=all_tracks_a,
  variable_name="thetao",
  new_var_name = "mid_temperature"
)

#### potential temperature (deep temperature) ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_deep.nc"), 
  target_df=all_tracks_a,
  variable_name="thetao",
  new_var_name = "deep_temperature"
)

#### shallow_sea_water_salinity  ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"), 
  target_df=all_tracks_a,
  variable_name="so",
  new_var_name = "shallow_sea_water_salinity"
)

#### mid_sea_water_salinity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_mid.nc"), 
  target_df=all_tracks_a,
  variable_name="so",
  new_var_name = "mid_sea_water_salinity"
)

#### deep_sea_water_salinity  ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_deep.nc"), 
  target_df=all_tracks_a,
  variable_name="so",
  new_var_name = "deep_sea_water_salinity"
)

#### shallow_eastward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"), 
  target_df=all_tracks_a,
  variable_name="uo",
  new_var_name = "shallow_eastward_sea_water_velocity"
)

#### shallow_northward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"), 
  target_df=all_tracks_a,
  variable_name="vo",
  new_var_name = "shallow_northward_sea_water_velocity"
)
#### mid_eastward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_mid.nc"), 
  target_df=all_tracks_a,
  variable_name="uo",
  new_var_name = "mid_eastward_sea_water_velocity"
)

#### mid_northward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_mid.nc"), 
  target_df=all_tracks_a,
  variable_name="vo",
  new_var_name = "mid_northward_sea_water_velocity"
)

#### deep_eastward_sea_water_velocity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_deep.nc"), 
  target_df=all_tracks_a,
  variable_name="uo",
  new_var_name = "deep_eastward_sea_water_velocity"
)

#### deep_northward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_deep.nc"), 
  target_df=all_tracks_a,
  variable_name="vo",
  new_var_name = "deep_northward_sea_water_velocity"
)

#### sea_bottom_temperature ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"), 
  target_df=all_tracks_a,
  variable_name="bottomT",
  new_var_name = "sea_bottom_temperature"
)

#### Sea bottom Salinity ready-----

Run_Extract<-Fxn_Extract_EnviroRaw_3D(
  # Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_bottSal_daily.nc"), 
  target_df=all_tracks_a,
  variable_name="sob",
  new_var_name = "bottomSal"
)

#### sea_surface_height ready####
Run_Extract<-Fxn_Extract_EnviroRaw_3D(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_phy_AllYr202324_vars_daily_surf.nc"),
  target_df=all_tracks_a,
  variable_name="zos",
  new_var_name = "sea_surface_height"
)


#### chloro ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2023_WD,"/cmems_mod_glo_bgc_my_AllYrs202324_nutrients_daily_surf.nc"),
  target_df=all_tracks_a,
  variable_name="chl",
  new_var_name = "chla"
)


### Save 2023-24 enviro dataframe dynamic ----

all_tracks_a
saveRDS(all_tracks_a, file = paste0(Results, "/SpWH_AllYr2023-24_EnviroDynamic_Pred_DF_wide.rds")) ##new one that is wide and includes all months for the year.

require(hms) #V2
EndTime_Dyn23=Sys.time()

RunTimeAll_Dyn23=EndTime_Dyn23 - StartTime_Dyn23;
TimeString_Run=paste0("Dynamic 2023 Start Time: ",
                      StartTime_Dyn23,
                      ";  End Time after saving files: ",
                      EndTime_Dyn23,
                      ";  Duration ",
                      as_hms(round(RunTimeAll_Dyn23,2))) ##TimeString
print(TimeString_Run)

### Data 2018 ----
library(dplyr)
all_tracks_b<-filter(pred_data_cl, year %in% c("2018"))
# all_tracks_b<-filter(all_tracks_dyn_pre, year %in% c("2018")) ##if subsetting data to add one variable

# all_tracks_b<-filter(all_tracks_b, sen_rep %in% c("0","1"))

NC_2018_WD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Environmental_variables/Products_2018"

require(hms) #V2
StartTime_Dyn18=Sys.time() ; StartTime
print(paste0("Start Time Dynamic 2018:", StartTime_Dyn18))

##Save Final Output ----
###If using subset data or merging dynamic with static ----
# all_tracks_combo<-cbind(all_tracks_dynamic, all_tracks_stat)
# all_tracks<-all_tracks_combo

###Continue with transformations -----
##Transformations first
###Make sure these are ordered in dynamic first followed by static variables
## all_tracks <- all_tracks %>%select(-log_chloro) friendly reminder for droping cols
all_tracks$log_sst<-log10(all_tracks$surface_temperature)
all_tracks$deep_EKE<-0.5* ((all_tracks$deep_eastward_sea_water_velocity^2)+(all_tracks$deep_northward_sea_water_velocity^2))
all_tracks$shallow_EKE<-0.5* ((all_tracks$shallow_eastward_sea_water_velocity^2)+(all_tracks$shallow_northward_sea_water_velocity^2))
all_tracks$mid_EKE<-0.5* ((all_tracks$mid_eastward_sea_water_velocity^2)+(all_tracks$mid_northward_sea_water_velocity^2))

all_tracks$log_deep_EKE<-log10(all_tracks$deep_EKE)
all_tracks$log_shallow_EKE<-log10(all_tracks$shallow_EKE)
all_tracks$log_mid_EKE<-log10(all_tracks$mid_EKE)
all_tracks$log_sea_bottom_temperature<-log10(all_tracks$sea_bottom_temperature)
all_tracks$log_botSal<-log10(all_tracks$bottomSal) ##for new variable

all_tracks$log_chloro<-log10(all_tracks$chla) ##we had forgotten about this in oct 2024

all_tracks$slope<-(all_tracks$slope)+0.0001
all_tracks$log_slope<-log10(all_tracks$slope)
all_tracks$sqrt_dist_coast<-sqrt(all_tracks$distcoast_km)
all_tracks$sqrt_dist_200m<-sqrt(all_tracks$dist200m_km)
all_tracks$sqrt_dist_500m<-sqrt(all_tracks$dist500m_km)
all_tracks$sqrt_dist_1000m<-sqrt(all_tracks$dist1000m_km) ##added mdr
all_tracks$sqrt_dist_smnt<-sqrt(all_tracks$distsmnt_km)
all_tracks$sqrt_dist_smntLRG<-sqrt(all_tracks$distsmntLRG_km)
all_tracks$sqrt_dist_smntSML<-sqrt(all_tracks$distsmntLRG_km)

all_tracks_p<-all_tracks
all_tracks<-all_tracks_p

# all_tracks <- all_tracks[,-c(17, 18, 20, 21,23,24,26,27,31,32,34,35,37,38)] #deletes columns 5 and 7 ##only apply to this if we have lat, long extra columns ##based on columns in colnames

# all_tracks <- all_tracks[,-c(21,22, 24,25,27,28,30,31, 35, 36, 38, 39, 41, 42)] #v2 for Nov2024 data ##note check again because we re-arranged some things

all_tracks <- all_tracks[,-c(22,23, 25,26,28,29,31,32, 36, 37, 39, 40,42, 43)] #v2 for Jan 2026 ##note check again because we re-arranged some things


# write.csv(all_tracks, file = paste0(Results,"/pred_data_env_data.csv"))
# saveRDS(all_tracks, file = paste0(Results,"/pred_data_env_data.rds"))
# saveRDS(all_tracks_p, file = paste0(Results,"/pred_data_env_data_BackUp.rds"))

##Jan 2026
# saveRDS(all_tracks, file = paste0(Results,"/pred_data_env_data_BotSal_Data2324.rds"))
# saveRDS(all_tracks, file = paste0(Results,"/pred_data_env_data_BotSal.rds"))
saveRDS(all_tracks, file = paste0(Results,"/pred_data_env_data_AllYrs_wide.rds"))
# saveRDS(all_tracks_p, file = paste0(Results,"/pred_data_env_data_BackUp_AllYrs_wide.rds"))


# all_tracks_stat<-subset(all_tracks, select=c(1:5,21:29, 39:46)) #c(21:29, 38, 40:46))
# saveRDS(all_tracks_stat, file = paste0(Results,"/pred_data_Static.rds"))

# StartTime=StartTime_Stat
library(hms)
EndSave=Sys.time()
RunTimeAll=EndSave - StartTime;
TimeString_RunAll=paste0("Download Environmental Data Start Time: ",
                         StartTime,
                         ";  End Time after saving files: ",
                         EndSave,
                         ";  Duration ",
                         as_hms(round(RunTimeAll,2))) ##TimeString
print(TimeString_RunAll)

 