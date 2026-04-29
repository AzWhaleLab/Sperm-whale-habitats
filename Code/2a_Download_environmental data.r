                        ### Part 2: download environmental data 

## Script for downloading environmental data.
  ##Environmental data for dynamic variables is mostly extracted from netcdf files and for static variables from other shapefiles. 
  ## See Methods in publication to understand better the different sources of variables used for this.
## Need to upload functions below to help extract variables
## Most of the script follows the same workflow as Hazen et al (2017) and Pérez-Jorge et al 2021. However, some items were turned into functions to help facilitate the extraction process, and many of the codes were updated to reflect updates to packages or transitions into new spatial packages from retired spatial packages. 
##In our workflow given that we have 2018, 2023, 2024, the workflow for extracting dynamic variables was split into 2018 data and 2023-2024 data. Then the data was combined to extract static variables. 
##Authors: Mareike D. Duffing Romero, Elliott Hazen, Sergi Pérez-Jorge. 
##There are notes and clarification throughout the document for seperate sections as needed, mostly defined as "##". Similarly, you can close/open different sections as needed using the navigation pane

##Load files ----
library(ncdf4)
library(lubridate)

rm(list = ls())

##Setting working directories 
SDM_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD) 
Results<-"~/Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/EnviroResults"  
Sergi_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi" ##place where some spatial shapefiles are located

##Load functions for extracting environmental data
Version="v1.0.R"#
source(paste0(SDM_HD, "/Functions/Function_SpatialMeansSD_NetCDF_", Version))

###Original tracks ----
original_tracks_pre<- readRDS("Data/SpermWhale_AllYrsTag_Predicted_6hrSSM_13Wh_v2.rds") 
original_tracks=original_tracks_pre #n=13 whales; structure is 
original_tracks$tagid<-original_tracks$id

original_tracks<-as.data.frame(original_tracks[c('tagid', 'date', 'lon', 'lat')]) #rename id to tagid mdr

#3.- Order the database hh3 by date and assign a new ID for the whole database. Split the date column to use it for names later 
#change names to latitude and longitude to use it in the loop

##eliminate track of the migratory sperm whale 244386
original_tracks <- original_tracks[which(!original_tracks$tagid =="244386"),] #changed to tagid

##MDR removed the code of flag, but added itireation, ID, sen_rep. Note ID and iteration are based on the same simulation id.
# original_tracks$flag <-0 
original_tracks$ID <-0
original_tracks$iteration <-0
original_tracks$sen_rep <-0

original_tracks$type<-"original"

### Pseudo- CRW----
#Open pseudo absences
pseudo_CRW_pre<- readRDS("~/Projects/SpermWhale/3_SpeciesDistMod/Pseudo-absences/Results/New_tracks_original/sim_CRW_all.rds") 
pseudo_CRW=pseudo_CRW_pre

colnames(pseudo_CRW)[1] <- "lon"
colnames(pseudo_CRW)[2] <- "lat"
colnames(pseudo_CRW)[3] <- "date"

pseudo_CRW
pseudo_CRW$type<-"CRW"
pseudo_CRW <- pseudo_CRW[, c("tagid","date", "lon", "lat", "ID", "iteration", "sen_rep", "type")]

###join all datasets ----
library(openxlsx)
tracksum<-read.xlsx("~/Projects/SpermWhale/1_SatelliteData/DataExploration/RawResults/Sperm Whale_Data_SummaryTable_AllYears_b_Oct2024.xlsx", sheet = "TagYear", detectDates = TRUE) ##Any Metadata you would like to include in your data set.
tracksum$tagid<-as.factor(tracksum$tagid)

## combine original track data with pseudo absences
all_tracks<-rbind(original_tracks,pseudo_CRW)
all_tracks$tagid<-as.factor(all_tracks$tagid)
##combine data with any metadata
all_tracks<-dplyr::left_join(all_tracks, tracksum, by="tagid")
all_tracks$datetime<-all_tracks$date

all_tracks$date<-as.Date(all_tracks$date)

#order dataset by date
all_tracks <- all_tracks[order(all_tracks$date),]

all_clean_pre<-all_tracks #or
alltracks<-all_tracks

## Download environmental data ----
## For each environmental extraction we use the Fxn_Extract_EnviroRaw function, which needs the following items for running:
    #input_file: nc file containing environmental data with time, lat/long/depth
    #target_df: the dataframe that has your tracking data and pseudo-absences
    #variable_name: the name mostly given by the source you acquired data from, for example Copernicus has thetao for temperature. 
    #new_var_name: the new variable name you want to give in your df, eg SST, SSH, etc

##Note, if you need to pre-inspect your netcdf file and gather information you can run the following code: 
  #nc<-nc_open(input_file)
  #print(nc)

##For seabottom salinity in 2018, we received a folder with several nc files that represent individual days for the timeframe we requested it. Thus, we applied the Fxn_Extract_EnviroRaw_MultiNC function instead


### Here we split data into 2018, 2023-2024 for extracting dynamic variables


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

all_tracks_a
saveRDS(all_tracks_a, file = paste0(Results, "/SpWH_2023-24_EnviroDynamic_DF_wide.rds")) 

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
all_tracks=all_clean_pre
all_tracks_b<-filter(all_tracks, TagYear %in% c("2018"))
# all_tracks_b<-filter(all_tracks_b, sen_rep %in% c("0","1"))

NC_2018_WD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Environmental_variables/Products_2018"

require(hms) #V2
StartTime_Dyn18=Sys.time() ; StartTime
print(paste0("Start Time Dynamic 2018:", StartTime_Dyn18))

#### potential temperature (sea surface temperature) ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"),
  target_df=all_tracks_b,
  variable_name="thetao",
  new_var_name = "surface_temperature"
)

#### potential temperature (mid temperature) ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_mid.nc"), 
  target_df=all_tracks_b,
  variable_name="thetao",
  new_var_name = "mid_temperature"
)

#### potential temperature (deep temperature) ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_deep.nc"), 
  target_df=all_tracks_b,
  variable_name="thetao",
  new_var_name = "deep_temperature"
)

#### shallow_sea_water_salinity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"), 
  target_df=all_tracks_b,
  variable_name="so",
  new_var_name = "shallow_sea_water_salinity"
)

#### mid_sea_water_salinity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_mid.nc"), 
  target_df=all_tracks_b,
  variable_name="so",
  new_var_name = "mid_sea_water_salinity"
)

#### deep_sea_water_salinity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_deep.nc"), 
  target_df=all_tracks_b,
  variable_name="so",
  new_var_name = "deep_sea_water_salinity"
)

#### shallow_eastward_sea_water_velocity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"),
  target_df=all_tracks_b,
  variable_name="uo",
  new_var_name = "shallow_eastward_sea_water_velocity"
)

#### shallow_northward_sea_water_velocity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"), 
  target_df=all_tracks_b,
  variable_name="vo",
  new_var_name = "shallow_northward_sea_water_velocity"
)

#### mid_eastward_sea_water_velocity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_mid.nc"), 
  target_df=all_tracks_b,
  variable_name="uo",
  new_var_name = "mid_eastward_sea_water_velocity"
)

#### mid_northward_sea_water_velocity ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_mid.nc"), 
  target_df=all_tracks_b,
  variable_name="vo",
  new_var_name = "mid_northward_sea_water_velocity"
)

#### deep_eastward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_deep.nc"), 
  target_df=all_tracks_b,
  variable_name="uo",
  new_var_name = "deep_eastward_sea_water_velocity"
)

#### deep_northward_sea_water_velocity ready####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_deep.nc"), 
  target_df=all_tracks_b,
  variable_name="vo",
  new_var_name = "deep_northward_sea_water_velocity"
)

#### sea_bottom_temperature  ready####
# Run_Extract<-Fxn_Extract_EnviroRaw_3D( #changed to 3D but could use something else
Run_Extract<-Fxn_Extract_EnviroRaw( 
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"), 
  target_df=all_tracks_b,
  variable_name="bottomT",
  new_var_name = "sea_bottom_temperature"
)

#### Sea bottom Salinity ready -----

##Using multiple nc files
all_tracks_dyn_pre <- Fxn_Extract_EnviroRaw_MultiNC(
  dir           = paste0(NC_2018_WD,"/Daily_BotSal2018"),
  variable_name = "sob",
  target_df     = all_tracks_b,
  new_var_name  = "bottomSal"
)
#### sea_surface_height ready####
Run_Extract<-Fxn_Extract_EnviroRaw_3D(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_phy_my_All2018_vars_daily_surf.nc"),
  target_df=all_tracks_b,
  variable_name="zos",
  new_var_name = "sea_surface_height"
)


#### chloro ready ####
Run_Extract<-Fxn_Extract_EnviroRaw(
  input_file = paste0(NC_2018_WD,"/cmems_mod_glo_bgc_my_AllYrs2018_nutrients_daily_deep.nc"),
  target_df=all_tracks_b,
  variable_name="chl",
  new_var_name = "chla"
)

### Save 2018 enviro dataframe dynamic ----

all_tracks_b
saveRDS(all_tracks_b, file = paste0(Results, "/SpWH_AllYr2018_EnviroDynamic_Pred_DF_wide.rds"))

require(hms) #V2
EndTime_Dyn18=Sys.time()
RunTimeAll_Dyn18=EndTime_Dyn18 - StartTime_Dyn18;
TimeString_RunB=paste0("Dynamic 2018 Start Time: ",
                       StartTime_Dyn18,
                       ";  End Time after saving files: ",
                       EndTime_Dyn18,
                       ";  Duration ",
                       as_hms(round(RunTimeAll_Dyn18,2))) ##TimeString
print(TimeString_RunB)

## Merge Dynamic 2018-2024 ----
rm(all_tracks)
all_tracks_a_pre<-readRDS(file = paste0(Results, "/SpWH_AllYr2023-24_EnviroDynamic_Pred_DF_wide.rds"))
all_tracks_a=all_tracks_a_pre

all_tracks_b_pre<-readRDS(file = paste0(Results, "/SpWH_AllYr2018_EnviroDynamic_Pred_DF_wide.rds"))
all_tracks_b=all_tracks_b_pre

##Compare missing variables
missing_in_b <- setdiff(names(all_tracks_a), names(all_tracks_b))
print(missing_in_b)
# all_tracks_b<-subset(all_tracks_b, select=-c(datetime))
all_tracks<-rbind(all_tracks_a, all_tracks_b)
all_tracks_dynamic<-all_tracks

saveRDS(all_tracks, file = paste0(Results, "/SpWH_AllYr_EnviroDynamic_Pred_DF_wide.rds"))

EndTime_Dyn=Sys.time()
RunTimeAll_Dyn=EndTime_Dyn - StartTime;
TimeString_RunDyn=paste0("Dynamic Start Time: ",
                         StartTime,
                         ";  End Time after saving files: ",
                         EndTime_Dyn,
                         ";  Duration ",
                         as_hms(round(RunTimeAll_Dyn,2))) ##TimeString
print(TimeString_RunDyn)

###Merging 2018-2024 data for subsetted data ----
# all_tracks_dynamic<-rbind(all_tracks_a, all_tracks_b) ##only if you go back and need to re-arrange add new variables.

## Static variables  ####
require(hms) #V2
StartTime_Stat=Sys.time() ; StartTime_Stat
print(paste0("Start Time Static:", StartTime_Stat))

library(sp); library(sf); library(raster)
plot(all_tracks$lat~all_tracks$lon)
dist<-all_tracks

coordinates(dist) <- ~lon+lat

proj4string(dist)<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# dist <-sp::spTransform(dist, CRS("+proj=utm +zone=26 +ellps=intl +towgs84=-104,167,-38,0,0,0,0 +units=m +no_defs"))  
dist <-sp::spTransform(dist, CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")) ##MDR updated this to keep same across the border

dist<-cbind(dist, dist@coords)
plot(dist$coords.x2~dist$coords.x1)

dist$lon<-dist$coords.x1
dist$lat<-dist$coords.x2

xypel<-dist[c("lon", "lat")]
## Convert xypel to sf SpatialPointsDataFrame
xypel_sf <- st_as_sf(xypel)

### DISTANCE TO COAST done----
print("Calculate Distance Coast")
coast <- shapefile(paste0(Sergi_HD, "/", "Shapefiles/","Azores_UTM26N.shp")) #or coast_SG
# coast <- st_read(paste0(Sergi_HD, "/", "Shapefiles/","Azores_UTM26N.shp"))
sp::proj4string(coast) <-CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Convert to coast to sf if it is a SpatialPolygonsDataFrame 
coast_sf <- st_as_sf(coast)

## Calculate the distance to the nearest coastline for each point
distcoast <- st_distance(xypel_sf, coast_sf)

## If coast_sf has multiple polygons, find the minimum distance for each point
distcoast_min <- apply(distcoast, 1, min)

## Convert distances to kilometers
distcoast_km <- distcoast_min / 1000

## Combine with the original data
dist_coast <- cbind(xypel_sf, distcoast_km)

## View the result
print(dist_coast)

## Extract dataframe
distcoast_df <- st_drop_geometry(dist_coast) #to get df

##Save files:
ItemName="dist_coast"
write.csv(distcoast_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_coast, output_path, append = FALSE) ##change append=true if you want to add on to the layer, but false is for overwriting the shapefile if it is already there

###DISTANCE TO 200M ISOBATH ready v3.3----
require(hms) #V2
StartTime_Stat=Sys.time() ; StartTime_Stat
print(paste0("Start Time Static:", StartTime_Stat))
print("Calculate Distance 200m")

# iso_200m <- shapefile(paste0(Sergi_HD, "/", "Shapefiles/","bat_UTM26N_200m.shp")) ##v3.0 with polygon
iso_200m <- shapefile(paste0(Sergi_HD, "/", "Shapefiles/Bat_CorrigidoLines/","bat_corrigido_PA.shp"));iso_200m ##new v3.3. new shapefile lines
sp::proj4string(iso_200m) <-CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# #dist_200m<-NULL 

## Convert to item to sf if it is a SpatialPolygonsDataFrame 
iso_200m_sf_pre <- st_as_sf(iso_200m)
iso_200m_sf<-iso_200m_sf_pre[iso_200m_sf_pre$Bat=="-200",]

## Calculate the distance to the nearest coastline for each point
dist200m <- st_distance(xypel_sf, iso_200m_sf)

## If iso_200m_sf has multiple polygons, find the minimum distance for each point
dist200m_min <- apply(dist200m, 1, min)

## Convert distances to kilometers
dist200m_km <- dist200m_min / 1000

## Combine with the original data
dist_200 <- cbind(xypel_sf, dist200m_km)

## View the result
print(dist_200)

## Extract dataframe
dist200m_df <- st_drop_geometry(dist_200) #to get df

##Save files:
ItemName="dist_200m"
write.csv(dist200m_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_200, output_path, append = FALSE)

###DISTANCE TO 500M ISOBATH v3.3----
print("Calculate Distance 500m")

# iso_500m_SG <- shapefile(paste0(Sergi_HD, "/", "Shapefiles/","bat_UTM26N_500m.shp"))
# iso_500m <- st_read(paste0(Sergi_HD, "/", "Shapefiles/","bat_UTM26N_500m.shp")) 
##since this is an sf object we need to skip some things below and convert to the crs we want
# sp::proj4string(iso_500m) <-CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Convert to item to sf if it is a SpatialPolygonsDataFrame 
# iso_500m_sf <- st_as_sf(iso_500m)
iso_500m_sf<-iso_200m_sf_pre[iso_200m_sf_pre$Bat=="-500",]

##Change projection to what we want for the sf object
st_crs(iso_500m_sf) <- "+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Calculate the distance to the nearest coastline for each point
dist500m <- st_distance(xypel_sf, iso_500m_sf)

## If iso_500m_sf has multiple polygons, find the minimum distance for each point
dist500m_min <- apply(dist500m, 1, min)

## Convert distances to kilometers
dist500m_km <- dist500m_min / 1000

## Combine with the original data
dist_500 <- cbind(xypel_sf, dist500m_km)

## View the result
print(dist_500)

## Extract dataframe
dist500m_df <- st_drop_geometry(dist_500) #to get df

##Save files:
ItemName="dist_500m"
write.csv(dist500m_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_500, output_path, append = FALSE)


###DISTANCE TO 1000M ISOBATH v3.3----
print("Calculate Distance 1000m")

## Convert to item to sf if it is a SpatialLinesDataFrame 

iso_1000m_sf<-iso_200m_sf_pre[iso_200m_sf_pre$Bat=="-1000",]

##Change projection to what we want for the sf object
st_crs(iso_1000m_sf) <- "+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Calculate the distance to the nearest coastline for each point
dist1000m <- st_distance(xypel_sf, iso_1000m_sf)

## If iso_1000m_sf has multiple polygons, find the minimum distance for each point
dist1000m_min <- apply(dist1000m, 1, min)

## Convert distances to kilometers
dist1000m_km <- dist1000m_min / 1000

## Combine with the original data
dist_1000 <- cbind(xypel_sf, dist1000m_km)

## View the result
print(dist_1000)

## Extract dataframe
dist1000m_df <- st_drop_geometry(dist_1000) #to get df

##Save files:
ItemName="dist_1000m"
write.csv(dist1000m_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_1000, output_path, append = FALSE)


#### Merge All tracks with Distance Coast to 1000m ready----
##MDR made changes here:  and added the dataframe from each calc
all_tracks<-cbind(all_tracks, distcoast_df, dist200m_df, dist500m_df, dist1000m_df)
saveRDS(all_tracks, file = paste0(Results, "/SpWH_All_EnviroDynamicStatDist_part1_DF.rds"))

#for filtered 
# all_tracks_filt2=all_tracks_filt
# all_tracks_combo<-cbind(all_tracks_filt2, distcoast_df, dist200m_df, dist500m_df, dist1000m_df)


###DISTANCE TO SEAMOUNTS ready----
print("Calculate Distance Seamounts")

# smnt <- shapefile("D:/R_temp/TrackS/SMNTSall.shp")#SMNTSall.shp=SMNTS_tab1all.shp
smnt <- shapefile(paste0(Sergi_HD, "/", "Shapefiles/","SMNTSall.shp"))#; st_crs(smnt)
sp::proj4string(smnt) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
smnts26n <-sp::spTransform(smnt, CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Convert to item to sf if it is a SpatialPolygonsDataFrame 
smnts26n_sf <- st_as_sf(smnts26n)

## Calculate the distance to the nearest coastline for each point
distsmnt <- st_distance(xypel_sf, smnts26n_sf)

## If smnts26n_sf has multiple polygons, find the minimum distance for each point
distsmnt_min <- apply(distsmnt, 1, min)

## Convert distances to kilometers
distsmnt_km <- distsmnt_min / 1000

## Combine with the original data
dist_smnt <- cbind(xypel_sf, distsmnt_km)

## View the result
print(dist_smnt)

## Extract dataframe
distsmnt_df <- st_drop_geometry(dist_smnt) #to get df
# print(distsmnt_df)

##Save files:
ItemName="dist_smnt"
write.csv(distsmnt_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_smnt, output_path, append = FALSE)

###DISTANCE TO SEAMOUNTS (large) ready----
print("Calculate Distance Seamounts Lrg")
## select seamounts within large seamounts
smnts26nLRG<-smnts26n[smnts26n$Category=="Large",] #distance all seamounts(461 smnts)

# dist_smntLRG<-NULL  

## Convert to item to sf if it is a SpatialPolygonsDataFrame 
smnts26nLRG_sf <- st_as_sf(smnts26nLRG)

## Calculate the distance to the nearest coastline for each point
distsmntLRG <- st_distance(xypel_sf, smnts26nLRG_sf)

## If smnts26nLRG_sf has multiple polygons, find the minimum distance for each point
distsmntLRG_min <- apply(distsmntLRG, 1, min)

## Convert distances to kilometers
distsmntLRG_km <- distsmntLRG_min / 1000

## Combine with the original data
dist_smntLRG <- cbind(xypel_sf, distsmntLRG_km)

## View the result
print(dist_smntLRG)

## Extract dataframe
distsmntLRG_df <- st_drop_geometry(dist_smntLRG) #to get df

# print(distsmntLRG_df)

##Save files:
ItemName="dist_smntLRG"
write.csv(distsmntLRG_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_smntLRG, output_path, append = FALSE)


###DISTANCE TO SEAMOUNTS (Small) ready----
print("Calculate Distance Seamounts Small")

## select seamounts within small seamounts
smnts26nSML<-smnts26n[smnts26n$Category=="Small",]

dist_smntSML<-NULL  #distance all seamounts(461 smnts)

## Convert to item to sf if it is a SpatialPolygonsDataFrame 
smnts26nSML_sf <- st_as_sf(smnts26nSML)

## Calculate the distance to the nearest coastline for each point
distsmntSML <- st_distance(xypel_sf, smnts26nSML_sf)

## If smnts26nSML_sf has multiple polygons, find the minimum distance for each point
distsmntSML_min <- apply(distsmntSML, 1, min)

## Convert distances to kilometers
distsmntSML_km <- distsmntSML_min / 1000

## Combine with the original data
dist_smntSML <- cbind(xypel_sf, distsmntSML_km)

## View the result
print(dist_smntSML)

## Extract dataframe
distsmntSML_df <- st_drop_geometry(dist_smntSML) #to get df

# print(distsmntSML_df)

##Save files:
ItemName="dist_smntSML"
write.csv(distsmntSML_df, file = paste0(Results,"/",ItemName,"_df.csv"))
output_path <- paste0(Results, "/", "StaticShapefiles", "/",ItemName, ".shp", sep="")
st_write(dist_smntSML, output_path, append = FALSE)

### Slope ready but require terra updates IP ----
print("Calculate Distance Slope and Depth")
# depth<- raster("C:/R/Azores/Batimetria/gebco_08_26N.img")#26N
depth<- raster(paste0(Sergi_HD, "/", "Shapefiles/","gebco_08_26N.img"))#26N

##
slope28620d<-terrain(depth, opt='slope', unit='degrees', neighbors=4, filename='slope28620d',overwrite=TRUE)


writeRaster(slope28620d,overwrite=TRUE, 'slope28620d.tif')
slope<-extract(slope28620d, xypel, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
               fun=NULL, na.rm=TRUE, layer, nl, df=FALSE, factors=FALSE)

#5) depth
#
depth<-extract(depth, xypel, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
               fun=NULL, na.rm=TRUE, layer, nl, df=FALSE, factors=FALSE)


all_tracks<-cbind(all_tracks, depth, slope, distsmnt_df,distsmntLRG_df,distsmntSML_df)
# all_tracks_combo<-cbind(all_tracks_combo, depth, slope, distsmnt_df,distsmntLRG_df,distsmntSML_df)

EndTimeStat=Sys.time()
print(paste0("End Time Static:", EndTimeStat))

RunTimeStat=EndTimeStat - StartTime_Stat

TimeString_Stat=paste0("Download Stat Data Start Time: ",
                       StartTime_Stat,
                       ";  End Time after saving files: ",
                       EndTimeStat,
                       ";  Duration ",
                       as_hms(round(RunTimeStat,2))) ##TimeString
print(TimeString_Stat)
##Save Final Output ----

# all_tracks<-cbind(all_tracks, dist200m_df,dist500m_df,dist1000m_df)
# all_tracks<- all_tracks[,-c(30,31)]
write.csv(all_tracks, file = paste0(Results,"/all_tracks_env_data.csv"))
saveRDS(all_tracks, file = paste0(Results,"/all_tracks_env_data.rds"))

###Create shapefile 
dist_all<-all_tracks
coordinates(dist_all) <- ~lon+lat
proj4string(dist_all)<- CRS("+proj=utm +zone=26 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

library(sf)
dist_all_sf <- st_as_sf(dist_all)

st_write(dist_all_sf,  paste0(Results, "/", "SpWh_alltracks_enviro.shp"))

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


##Save final output with Dynamic additions ----
all_tracks_combo<-cbind(all_tracks_dynamic, all_tracks_stat)
saveRDS(all_tracks_combo, file = paste0(Results,"/all_tracks_env_data_withLagChloro.rds"))

## Save final output with Bottom salinity ----
all_tracks_combo<-cbind(all_tracks_dynamic, all_tracks_stat)
# saveRDS(all_tracks_combo, file = paste0(Results,"/all_tracks_env_data_wBotSal_DataYrs2324.rds"))
saveRDS(all_tracks_combo, file = paste0(Results,"/all_tracks_env_data_wBotSal.rds"))

