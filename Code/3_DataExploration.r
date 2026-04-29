###                                 Data exploration ###
##Script 3: Data exploration for assessing collinereaty, transform variables and then select variables to be used for modeling approaches.

## Load Files ----
#library
library(pastecs)
library(lattice)
library(psych)
library(fitdistrplus)
library(logspline)
library(lubridate)

#DATABASE#uw
#database
#Attention#mercator files #CHLA and NPP have NA value# these were clean before to include in our database!!! 
##MDR added a save code of PA variables 
##Setting working directories 
SDM_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD) ## change this in work computer
# SDM_HD<-"C:/Mareike/PhD Azores/Data/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD) ## change this in work computer

Results<-"~/Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/EnviroResults"  ## change this in work computer
# Results<-"C:/Mareike/PhD Azores/Data/Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/EnviroResults"  ## change this in work computer 

Sergi_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi" ##MDR added and change this on the work computer
# Sergi_HD<-"C:/Mareike/PhD Azores/Data/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi" ##MDR added and change this on the work computer

# FFDB<- read.csv(file = "all_tracks_env_data.csv",sep = ",",dec = ".", header = TRUE)
# FFDB_pre<- readRDS(paste0(Results,"/all_tracks_env_data.rds"))
FFDB_pre<- readRDS(paste0(Results,"/all_tracks_env_data.rds"))
FFDB<-FFDB_pre
str(FFDB)

# Check the column names to confirm duplicates are removed
colnames(FFDB)
boxplot(FFDB$depth)

## Data cleaning Prep ----
#cleaning NA
#FFDB$SLPd<-FFDB$Slope26N
FFDB<-FFDB[FFDB$depth<=0,]## 
FFDB<-FFDB[complete.cases(FFDB$depth),]
FFDB<-FFDB[FFDB$depth!=0,] 
FFDB$depth<-FFDB$depth*-1
FFDB<-FFDB[FFDB$distcoast_km!=0,] ##updated the variable name 


FFDB$date<-as.Date(FFDB$date)

FFDB$MONTH<-month(as.POSIXlt(FFDB$date, format="%d/%m/%Y"))

FFDB$fMONTH <- factor(FFDB$MONTH, levels = c(2,3,4,5, 6, 7,8,9,10),
                      labels = c("February","March","April","May", "June","July","August", "September","October")) ## change this to basic as date=month



FFDB$shallow_EKE<-0.5* ((FFDB$shallow_eastward_sea_water_velocity^2)+(FFDB$shallow_northward_sea_water_velocity^2))

FFDB$mid_EKE<-0.5* ((FFDB$mid_eastward_sea_water_velocity^2)+(FFDB$mid_northward_sea_water_velocity^2))

FFDB$deep_EKE<-0.5* ((FFDB$deep_eastward_sea_water_velocity^2)+(FFDB$deep_northward_sea_water_velocity^2))
##variables
tagid<-FFDB$tagid
tagid2<-as.character(FFDB$tagid)
Lat<-FFDB$lat 
Long<-FFDB$lon
# date<-FFDB$date
date<-as.Date(FFDB$date)
year<-year(as.POSIXlt(FFDB$date, format="%d/%m/%Y"))
month<-FFDB$MONTH
sen_rep<-FFDB$sen_rep
ID<-FFDB$ID ##this means itiration from sim
type<-FFDB$type
surface_temp<-FFDB$surface_temperature ##updated from sst to surface temperature 
mid_temp<-FFDB$mid_temperature
deep_temp<-FFDB$deep_temperature ##updated from sst to surface temperature 

shallow_sea_water_salinity<-FFDB$shallow_sea_water_salinity
mid_sea_water_salinity<-FFDB$mid_sea_water_salinity
deep_sea_water_salinity<-FFDB$deep_sea_water_salinity

shallow_eastward_sea_water_velocity<-FFDB$shallow_eastward_sea_water_velocity
shallow_northward_sea_water_velocity<-FFDB$shallow_northward_sea_water_velocity
shallow_EKE<-FFDB$shallow_EKE

mid_eastward_sea_water_velocity<-FFDB$mid_eastward_sea_water_velocity
mid_northward_sea_water_velocity<-FFDB$mid_northward_sea_water_velocity
mid_EKE<-FFDB$mid_EKE

deep_eastward_sea_water_velocity<-FFDB$deep_eastward_sea_water_velocity
deep_northward_sea_water_velocity<-FFDB$deep_northward_sea_water_velocity
deep_EKE<-FFDB$deep_EKE

sea_bottom_temperature<-FFDB$sea_bottom_temperature
sea_surface_height<-FFDB$sea_surface_height
dist_coast<-FFDB$distcoast_km
dist_200m<-FFDB$dist200m_km
dist_500m<-FFDB$dist500m_km
dist_1000m<-FFDB$dist1000m_km
depth<-FFDB$depth
slope<-FFDB$slope
dist_smnt<-FFDB$distsmnt_km
dist_smntLRG<-FFDB$distsmntLRG_km
dist_smntSML<-FFDB$distsmntSML_km
chloro<-FFDB$chla
## Data visualization ----

library(ggplot2)

cols <- c("#F76D5E", "#FFFFBF", "#72D8FF")
# cols <- c("#72D8FF", "red")
# Basic density plot in ggplot2
### Surface temperature ----
##updated from sst
##Note: we pulled out data from two stations, and thus we have a bimodal distribution
ItemVar="SurfaceTemp"
hist(surface_temp)
hist(sqrt(surface_temp))
hist(log10(surface_temp))

dens_temp<-ggplot(FFDB, aes(x = surface_temperature, fill = type)) + ##changed sst variable to surface temp
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_temp

freq_temp<-ggplot(FFDB, aes(x=surface_temperature, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_temp

ggsave(dens_temp, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
# ggsave(freq_temp, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)


### Mid temperature ----

ItemVar="MidTemp"
hist(mid_temp)
hist(sqrt(mid_temp))
hist(log10(mid_temp))

dens_tempMid<-ggplot(FFDB, aes(x = mid_temp, fill = type)) + ##changed sst variable to surface temp
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_tempMid

freq_tempMid<-ggplot(FFDB, aes(x=mid_temp, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_tempMid

ggsave(dens_tempMid, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
# ggsave(freq_tempMid, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)



###Deep temperature----
ItemVar="DeepTemp"

hist(deep_temp) # we may not need to transform since it is normal
hist(sqrt(deep_temp))
hist(log10(deep_temp))

dens_dtemp<-ggplot(FFDB, aes(x = deep_temp, fill = type)) + ##changed variable to deep temp
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_dtemp

freq_dtemp<-ggplot(FFDB, aes(x=deep_temp, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_dtemp

ggsave(dens_dtemp, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
ggsave(freq_dtemp, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)

###shallow_sea_water_salinity ----
ItemVar="ShallowSal"
hist(shallow_sea_water_salinity)
hist(sqrt(shallow_sea_water_salinity))
hist(log10(shallow_sea_water_salinity))

dens_plo<-ggplot(FFDB, aes(x = shallow_sea_water_salinity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=shallow_sea_water_salinity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)


FFDB$surface_sea_water_salinity<-FFDB$shallow_sea_water_salinity ##added April 2026
dens_plo<-ggplot(FFDB, aes(x = surface_sea_water_salinity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo
ItemVar="SurfaceSal"
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")


###mid_sea_water_salinity ----
ItemVar="MidSal"
hist(mid_sea_water_salinity)
hist(sqrt(mid_sea_water_salinity))
hist(log10(mid_sea_water_salinity))

dens_plo<-ggplot(FFDB, aes(x = mid_sea_water_salinity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=mid_sea_water_salinity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)

###deep_sea_water_salinity----
ItemVar="DeepSal"

hist(deep_sea_water_salinity)
hist(sqrt(deep_sea_water_salinity))
hist(log10(deep_sea_water_salinity))

dens_sal<-ggplot(FFDB, aes(x = deep_sea_water_salinity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_sal

freq_sal<-ggplot(FFDB, aes(x=deep_sea_water_salinity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_sal

ggsave(dens_sal, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
ggsave(freq_sal, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)


### shallow_eastward_sea_water_velocity----
ItemVar="ShallowEastW_SeaWatVel"
hist(shallow_eastward_sea_water_velocity)
hist(sqrt(shallow_eastward_sea_water_velocity))
hist(log10(shallow_eastward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = shallow_eastward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=shallow_eastward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

FFDB$surface_eastward_sea_water_velocity<-FFDB$shallow_eastward_sea_water_velocity ##added April 2026
dens_plo<-ggplot(FFDB, aes(x = surface_eastward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo
ItemVar="Surface_EW_vel"
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###shallow_northward_sea_water_velocity----
ItemVar="ShallowNorthW_SeaWatVel"
hist(shallow_northward_sea_water_velocity)
hist(sqrt(shallow_northward_sea_water_velocity))
hist(log10(shallow_northward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = shallow_northward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=shallow_northward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

FFDB$surface_northward_sea_water_velocity<-FFDB$shallow_northward_sea_water_velocity ##added April 2026
dens_plo<-ggplot(FFDB, aes(x = surface_northward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo
ItemVar="Surface_NW_vel"
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

### Shallow EKE (transform to log)----
ItemVar="ShallowEKE"

hist(shallow_EKE)
hist(sqrt(shallow_EKE))
hist(log10(shallow_EKE))

dens_plo<-ggplot(FFDB, aes(x = shallow_EKE, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=shallow_EKE, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

FFDB$surface_EKE<-FFDB$shallow_EKE ##added April 2026
dens_plo<-ggplot(FFDB, aes(x = surface_EKE, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo
ItemVar="Surface_EKE"
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

dens_plo<-ggplot(FFDB, aes(x = surface_EKE, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(fill = "Presence Type")
dens_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot_bottomLegend.tif"), width=7.0, height = 5, bg="white", compression="lzw")


### mid_eastward_sea_water_velocity----
ItemVar="MidEastW_SeaWatVel"
hist(mid_eastward_sea_water_velocity)
hist(sqrt(mid_eastward_sea_water_velocity))
hist(log10(mid_eastward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = mid_eastward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=mid_eastward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
###mid_northward_sea_water_velocity----
ItemVar="MidNorthW_SeaWatVel"
hist(mid_northward_sea_water_velocity)
hist(sqrt(mid_northward_sea_water_velocity))
hist(log10(mid_northward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = mid_northward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=mid_northward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

### Mid EKE (transform log)----
ItemVar="MidEKE"

hist(shallow_EKE)
hist(sqrt(shallow_EKE))
hist(log10(shallow_EKE))

dens_plo<-ggplot(FFDB, aes(x = mid_EKE, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=mid_EKE, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")






###deep_eastward_sea_water_velocity----
ItemVar="DeepEastW_SeaWatVel"
hist(deep_eastward_sea_water_velocity) ##more normal
hist(sqrt(deep_eastward_sea_water_velocity))
hist(log10(deep_eastward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = deep_eastward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=deep_eastward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw",)

###deep_northward_sea_water_velocity----
ItemVar="DeepNorthW_SeaWatVel"
hist(deep_northward_sea_water_velocity) ##more normal
hist(sqrt(deep_northward_sea_water_velocity))
hist(log10(deep_northward_sea_water_velocity))

dens_plo<-ggplot(FFDB, aes(x = deep_northward_sea_water_velocity, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=deep_northward_sea_water_velocity, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###deep_EKE (needs transformation Log)----
##left skewed ori
ItemVar="DeepEKE"
hist(deep_EKE)
hist(sqrt(deep_EKE))
hist(log10(deep_EKE))

dens_plo<-ggplot(FFDB, aes(x = deep_EKE, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=deep_EKE, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")



###sea_bottom_temperature (transform to log)----
ItemVar="BottomTemp"
hist(sea_bottom_temperature) 
hist(sqrt(sea_bottom_temperature))
hist(log10(sea_bottom_temperature)) 

dens_plo<-ggplot(FFDB, aes(x = sea_bottom_temperature, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=sea_bottom_temperature, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

##need transformation to log

### sea_surface_height ----
ItemVar="SeaSurfaceHeight"
hist(sea_surface_height)
hist(sqrt(sea_surface_height))
hist(log10(sea_surface_height))

dens_plo<-ggplot(FFDB, aes(x = sea_surface_height, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

freq_plo<-ggplot(FFDB, aes(x=sea_surface_height, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###chlorophyll (log transform)-----
ItemVar="Chloro"
hist(chloro)
hist(sqrt(chloro))
hist(log10(chloro))

dens_plo<-ggplot(FFDB, aes(x = chla, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw() #or theme_classic
dens_plo

freq_plo<-ggplot(FFDB, aes(x=chla, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")



###dist_coast (transform to  sqrt)----
ItemVar="DistCoast"
hist(dist_coast) ##left skewed
hist(sqrt(dist_coast))
hist(log10(dist_coast))

dens_plo<-ggplot(FFDB, aes(x = dist_coast, fill = type)) +
  geom_density(alpha = 0.7) +  
  scale_fill_manual(values = cols)+
  theme_bw()+
  ylab("")+
  theme(legend.position = "none")
dens_plo

dens_plo_byIndi<-ggplot(FFDB, aes(x = dist_coast, fill = type)) +
  geom_density(alpha = 0.7) +  
  scale_fill_manual(values = cols)+
  facet_wrap(~tagid)+
  theme_bw()
dens_plo_byIndi

freq_plo<-ggplot(FFDB, aes(x=dist_coast, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

ggsave(dens_plo_byIndi, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot_byIndi.tif"),  width=7.0, height = 6.5, bg="white", compression="lzw") #width=7.0, height = 5, or 7.17 x 6.43 in image

ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

##need transformation sqrt best

###dist_200m ( transform to sqrt)----
ItemVar="Dist200m"
hist(dist_200m) ##left skewed
hist(sqrt(dist_200m))
hist(log10(dist_200m))

dens_plo<-ggplot(FFDB, aes(x = dist_200m, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_200m, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

##need transformation
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###dist_500m (possibly transform to sqrt)----
ItemVar="Dist500m"
hist(dist_500m) ##left skewed
hist(sqrt(dist_500m))
hist(log10(dist_500m))

dens_plo<-ggplot(FFDB, aes(x = dist_500m, fill = type)) +
  # ggplot(FFDB, aes(x = sqrt(dist_500m), fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_500m, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo
##need transformation
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###dist_1000m ----
ItemVar="Dist1000m"
hist(dist_1000m)
hist(sqrt(dist_1000m))
hist(log10(dist_1000m))

dens_plo<-ggplot(FFDB, aes(x = dist_1000m, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_1000m, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###depth ----
ItemVar="Depth"
hist(depth) ##semi normal dist
hist(sqrt(depth))
hist(log10(depth))

dens_plo<-ggplot(FFDB, aes(x = depth, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=depth, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###slope (need to transform, log)----
ItemVar="Slope"
hist(slope) ##left skewed
hist(sqrt(slope))
hist(log10(slope))

dens_plo<-ggplot(FFDB, aes(x = slope, fill = type)) +
  # ggplot(FFDB, aes(x = log10(slope), fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=slope, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo
##need transformation
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###dist_smnt (need to transform to sqrt poss)----
ItemVar="DistSeaMnt"
hist(dist_smnt) ##left skewed
hist(sqrt(dist_smnt))
hist(log10(dist_smnt))

dens_plo<-ggplot(FFDB, aes(x = dist_smnt, fill = type)) +
  geom_density(alpha = 0.7) +  
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_smnt, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo

##need transformation
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###dist_smntLRG (need to transform, sqrt best)----
ItemVar="DistSeaMnt_LRG"
hist(dist_smntLRG) ##left skewed
hist(sqrt(dist_smntLRG))
hist(log10(dist_smntLRG))

dens_plo<-ggplot(FFDB, aes(x = dist_smntLRG, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_smntLRG, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo
##need transformation
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")

###dist_smntSML(need to transform sqrt)-----
ItemVar="DistSeaMnt_SML"
hist(dist_smntSML)
hist(sqrt(dist_smntSML))
hist(log10(dist_smntSML))

dens_plo<-ggplot(FFDB, aes(x = dist_smntSML, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)+
  theme_bw()
dens_plo

freq_plo<-ggplot(FFDB, aes(x=dist_smntSML, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = cols) +
  theme_bw()
freq_plo
ggsave(dens_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_DensityPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")
ggsave(freq_plo, filename=paste0(Results, "/", "Enviro_PA_SumResults/ForPub", "/", ItemVar, "_FreqPlot.tif"), width=7.0, height = 5, bg="white", compression="lzw")


## Transformations----
### For our data Sperm whales (MDR)----
##Need to transform: Surface temperature,  deep eke, shallow eke, sea bottom temp, dist coast, dist 200, dist 500, slope, smnt, smnt Lrg, smnt small.
###Review newest Dec2024 data and transform new values. Don't forget to include the vars for all depth ranges 
# log_sst<-log10(FFDB$sst)
log_sst<-log10(FFDB$surface_temperature)

log_shallow_EKE<-log10(FFDB$shallow_EKE)
log_mid_EKE<-log10(FFDB$mid_EKE)
log_deep_EKE<-log10(FFDB$deep_EKE)

log_sea_bottom_temperature<-log10(FFDB$sea_bottom_temperature)
FFDB$slope<-(FFDB$slope)+0.0001
log_slope<-log10(FFDB$slope)
log_chloro<-log10(FFDB$chla) ##we had forgotten about this.

sqrt_dist_coast<-sqrt(FFDB$distcoast_km)
sqrt_dist_200m<-sqrt(FFDB$dist200m_km)
sqrt_dist_500m<-sqrt(FFDB$dist500m_km)
sqrt_dist_1000m<-sqrt(FFDB$dist1000m_km) ##added mdr
sqrt_dist_smnt<-sqrt(FFDB$distsmnt_km)
sqrt_dist_smntLRG<-sqrt(FFDB$distsmntLRG_km)
sqrt_dist_smntSML<-sqrt(FFDB$distsmntLRG_km)

##Covariate analysis----
#1)desc stat
#install.packages("pastecs")
library(pastecs)
###1_All_StatDyn ----
##need to add ID and presence
##remove upward sea water velocity; add deep temp (mdr)
scores_sel<-as.data.frame(cbind(tagid, tagid2,Lat, Long, date,year,month,sen_rep,type,
                                log_sst,
                                deep_temp,
                                mid_temp,
                                shallow_sea_water_salinity,
                                mid_sea_water_salinity,
                                deep_sea_water_salinity,
                                shallow_eastward_sea_water_velocity,
                                shallow_northward_sea_water_velocity,
                                mid_eastward_sea_water_velocity, 
                                mid_northward_sea_water_velocity,
                                deep_eastward_sea_water_velocity,
                                deep_northward_sea_water_velocity,
                                log_shallow_EKE,
                                log_mid_EKE,
                                log_deep_EKE,
                                log_sea_bottom_temperature,
                                sea_surface_height, 
                                sqrt_dist_coast, 
                                sqrt_dist_200m,
                                sqrt_dist_500m, 
                                sqrt_dist_1000m, 
                                depth,
                                log_slope,
                                sqrt_dist_smnt, 
                                sqrt_dist_smntLRG, 
                                sqrt_dist_smntSML, 
                                log_chloro))

# write.csv(scores_sel,file="scores_sel.csv")
# write.csv(scores_sel,file=paste0(Results, "/Enviro_scores_sel.csv"))

###environmental and SEAPODYM variables
##Convert some columns to numeric (mdr)
##For correlation and output: a) all variables b) no deep values (ie. deep temp, deep salinity, deep eke)
library(dplyr)
scores_sel <- scores_sel %>% mutate_at(c("log_sst","mid_temp", "deep_temp","deep_sea_water_salinity", "mid_sea_water_salinity", "shallow_sea_water_salinity", "mid_eastward_sea_water_velocity", "mid_northward_sea_water_velocity",  "deep_eastward_sea_water_velocity","deep_northward_sea_water_velocity","shallow_eastward_sea_water_velocity", "shallow_northward_sea_water_velocity", "log_deep_EKE", "log_mid_EKE","log_shallow_EKE", "log_sea_bottom_temperature", "sea_surface_height","sqrt_dist_coast", "sqrt_dist_200m", "sqrt_dist_500m", "sqrt_dist_1000m","depth", "log_slope", "sqrt_dist_smnt", "sqrt_dist_smntLRG", "sqrt_dist_smntSML", "log_chloro"), as.numeric)
detach("package:dplyr", unload = TRUE)

corPearson_scores_sel<-cor(scores_sel[,c("log_sst", 
                                         "mid_temp",
                                         "deep_temp",
                                         "shallow_sea_water_salinity", 
                                         "mid_sea_water_salinity",     
                                         "deep_sea_water_salinity", 
                                         "shallow_eastward_sea_water_velocity", 
                                         "shallow_northward_sea_water_velocity",
                                         "mid_eastward_sea_water_velocity",
                                         "mid_northward_sea_water_velocity",
                                         "deep_eastward_sea_water_velocity",
                                         "deep_northward_sea_water_velocity",
                                         "log_shallow_EKE", 
                                         "log_mid_EKE",
                                         "log_deep_EKE",
                                         "log_sea_bottom_temperature", 
                                         "sea_surface_height",
                                         "sqrt_dist_coast", 
                                         "sqrt_dist_200m", 
                                         "sqrt_dist_500m", 
                                         "sqrt_dist_1000m",
                                         "depth", 
                                         "log_slope",
                                         "sqrt_dist_smnt", 
                                         "sqrt_dist_smntLRG", 
                                         "sqrt_dist_smntSML", 
                                         "log_chloro")],
                           # use="complete"
                           # use="complete.obs" #second attempt
                           use = "pairwise.complete.obs"
)

corPearson_scores_sel
# write.csv(corPearson_scores_sel,file="corPearson_scores.csv")
# write.csv(corPearson_scores_sel,file=paste0(Results, "/corPearson_scores.csv"))
# write.csv(corPearson_scores_sel,file=paste0(Results, "/corPearson_scores_v2.csv"))

##Aditional step to change NAs to Zero and run correlation matrix. (from chatgpt)
# Replace NA correlations with 0
corPearson_scores_sel[is.na(corPearson_scores_sel)] <- 0


##
signif(corPearson_scores_sel, digits=4)
c<-corPearson_scores_sel


###2_Dynamic values ----
scores_sel_dyn<-as.data.frame(cbind(tagid, tagid2,Lat, Long, date,year,month,sen_rep,type,
                                    log_sst,
                                    deep_temp,
                                    deep_sea_water_salinity,
                                    shallow_sea_water_salinity,
                                    deep_eastward_sea_water_velocity,
                                    deep_northward_sea_water_velocity,
                                    shallow_eastward_sea_water_velocity,
                                    shallow_northward_sea_water_velocity,
                                    log_deep_EKE,
                                    log_shallow_EKE,
                                    log_sea_bottom_temperature,
                                    sea_surface_height, 
                                    # sqrt_dist_coast, 
                                    # sqrt_dist_200m,
                                    # sqrt_dist_500m, 
                                    # sqrt_dist_1000m, 
                                    # depth,
                                    # log_slope,
                                    # sqrt_dist_smnt, 
                                    # sqrt_dist_smntLRG, 
                                    # sqrt_dist_smntSML, 
                                    log_chloro))

# write.csv(scores_sel_dyn,file="scores_sel_dyn.csv")
# write.csv(scores_sel_dyn,file=paste0(Results, "/Enviro_scores_sel_dyn.csv"))

library(dplyr)
scores_sel_dyn <- scores_sel_dyn %>% mutate_at(c("log_sst", "deep_temp","deep_sea_water_salinity", "shallow_sea_water_salinity", "deep_eastward_sea_water_velocity","deep_northward_sea_water_velocity","shallow_eastward_sea_water_velocity", "shallow_northward_sea_water_velocity", "log_deep_EKE","log_shallow_EKE", "log_sea_bottom_temperature", "sea_surface_height",
                                                 # "sqrt_dist_coast", 
                                                 # "sqrt_dist_200m", "sqrt_dist_500m", "sqrt_dist_1000m","depth", "log_slope", "sqrt_dist_smnt", "sqrt_dist_smntLRG", "sqrt_dist_smntSML", 
                                                 "log_chloro"), as.numeric)
detach("package:dplyr", unload = TRUE)

corPearson_scores_sel_dyn<-cor(scores_sel_dyn[,c("log_sst", 
                                                 "deep_temp",
                                                 "deep_sea_water_salinity", 
                                                 "shallow_sea_water_salinity", 
                                                 "deep_eastward_sea_water_velocity",
                                                 "deep_northward_sea_water_velocity",
                                                 "shallow_eastward_sea_water_velocity", 
                                                 "shallow_northward_sea_water_velocity",
                                                 "log_deep_EKE",
                                                 "log_shallow_EKE", 
                                                 "log_sea_bottom_temperature", 
                                                 "sea_surface_height",
                                                 # "sqrt_dist_coast", 
                                                 # "sqrt_dist_200m", 
                                                 # "sqrt_dist_500m", 
                                                 # "sqrt_dist_1000m",
                                                 # "depth", 
                                                 # "log_slope",
                                                 # "sqrt_dist_smnt", 
                                                 # "sqrt_dist_smntLRG", 
                                                 # "sqrt_dist_smntSML", 
                                                 "log_chloro")],use="complete")

corPearson_scores_sel_dyn
# write.csv(corPearson_scores_sel_dyn,file="corPearson_scores_Dyn.csv")
# write.csv(corPearson_scores_sel_dyn,file=paste0(Results, "/corPearson_scores_Dyn.csv"))

signif(corPearson_scores_sel_dyn, digits=4)
d<-corPearson_scores_sel_dyn

###3_Static values ----
scores_sel_stat<-as.data.frame(cbind(tagid, tagid2,Lat, Long, date,year,month,sen_rep,type,
                                     # log_sst,
                                     # deep_temp,
                                     # deep_sea_water_salinity,
                                     # shallow_sea_water_salinity,
                                     # deep_eastward_sea_water_velocity,
                                     # deep_northward_sea_water_velocity,
                                     # shallow_eastward_sea_water_velocity,
                                     # shallow_northward_sea_water_velocity,
                                     # log_deep_EKE,
                                     # log_shallow_EKE,
                                     # log_sea_bottom_temperature,
                                     # sea_surface_height, 
                                     sqrt_dist_coast,
                                     sqrt_dist_200m,
                                     sqrt_dist_500m,
                                     sqrt_dist_1000m,
                                     depth,
                                     log_slope,
                                     sqrt_dist_smnt,
                                     sqrt_dist_smntLRG,
                                     sqrt_dist_smntSML#,
                                     # log_chloro
))

# #write.csv(scores_sel_stat,file="scores_sel_dyn.csv")
# write.csv(scores_sel_stat,file=paste0(Results, "/Enviro_scores_sel_stat.csv"))

library(dplyr)
scores_sel_stat <- scores_sel_stat %>% mutate_at(c(
  # "log_sst", "deep_temp","deep_sea_water_salinity", "shallow_sea_water_salinity", "deep_eastward_sea_water_velocity","deep_northward_sea_water_velocity","shallow_eastward_sea_water_velocity", "shallow_northward_sea_water_velocity", "log_deep_EKE","log_shallow_EKE", "log_sea_bottom_temperature", "sea_surface_height",
  "sqrt_dist_coast",
  "sqrt_dist_200m", "sqrt_dist_500m", "sqrt_dist_1000m","depth", "log_slope", "sqrt_dist_smnt", "sqrt_dist_smntLRG", "sqrt_dist_smntSML"#,
  # "log_chloro"
), as.numeric)
detach("package:dplyr", unload = TRUE)

corPearson_scores_sel_stat<-cor(scores_sel_stat[,c(
  # "log_sst", 
  # "deep_temp",
  # "deep_sea_water_salinity", 
  # "shallow_sea_water_salinity", 
  # "deep_eastward_sea_water_velocity",
  # "deep_northward_sea_water_velocity",
  # "shallow_eastward_sea_water_velocity", 
  # "shallow_northward_sea_water_velocity",
  # "log_deep_EKE",
  # "log_shallow_EKE", 
  # "log_sea_bottom_temperature", 
  # "sea_surface_height",
  "sqrt_dist_coast",
  "sqrt_dist_200m",
  "sqrt_dist_500m",
  "sqrt_dist_1000m",
  "depth",
  "log_slope",
  "sqrt_dist_smnt",
  "sqrt_dist_smntLRG",
  "sqrt_dist_smntSML"#,
  # "log_chloro"
)],use="complete")

corPearson_scores_sel_stat
# #write.csv(corPearson_scores_sel_dyn,file="corPearson_scores_Dyn.csv")
# write.csv(corPearson_scores_sel_stat,file=paste0(Results, "/corPearson_scores_Stat.csv"))

signif(corPearson_scores_sel_stat, digits=4)
s<-corPearson_scores_sel_stat


### Assessing correlations ----
## From Maria Ines 
library(corrplot)
library(RColorBrewer)
library(wesanderson)
#calculate correlation p-values !mat is a matrix of data!
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
df<-scores_sel
df.m<-as.matrix(df[,c("log_sst", 
                      "mid_temp",
                      "deep_temp",
                      "shallow_sea_water_salinity", 
                      "mid_sea_water_salinity",     
                      "deep_sea_water_salinity", 
                      "shallow_eastward_sea_water_velocity", 
                      "shallow_northward_sea_water_velocity",
                      "mid_eastward_sea_water_velocity",
                      "mid_northward_sea_water_velocity",
                      "deep_eastward_sea_water_velocity",
                      "deep_northward_sea_water_velocity",
                      "log_shallow_EKE", 
                      "log_mid_EKE",
                      "log_deep_EKE",
                      "log_sea_bottom_temperature", 
                      "sea_surface_height",
                      "sqrt_dist_coast", 
                      "sqrt_dist_200m", 
                      "sqrt_dist_500m", 
                      "sqrt_dist_1000m",
                      "depth", 
                      "log_slope",
                      "sqrt_dist_smnt", 
                      "sqrt_dist_smntLRG", 
                      "sqrt_dist_smntSML", 
                      "log_chloro")])
##Extra code from chatgpt v2 Jan2024
# Convert to numeric and handle NAs/Inf
# df.m <- as.data.frame(lapply(df.m, function(x) ifelse(is.finite(x), x, NA)))

# Ensure all rows with NA are removed consistently
# df.m <- na.omit(df.m)
# c <- cor(df.m, use = "pairwise.complete.obs")

##
p.cor<-cor.mtest(df.m)
pal <- wes_palette("Zissou1", 20, type = "continuous")
##Version Ines
# corrplot(c, method="color", type="lower", order="hclust",
#         col=pal,p.mat=p.cor,sig.level = 0.05,insig = "blank",
#         diag=F,tl.col = 'black',tl.srt = 90,tl.cex = .65,cl.pos="r",cl.cex = .7,cl.ratio = .1) ##tl.cex was0.5 

##c above is whatever stat variance test your run 

##Other code from chgpt
plot.new(); dev.off()  # Reset plotting device
# par(mar = c(5,5,5,5))  # Increase margins
par(mar = c(5, 5, 2, 10))  # Increase right margin for legend spacing


corrplot(c, method="color", type="lower", order="hclust",
         col=pal, p.mat=p.cor, sig.level=0.05, insig="blank",
         diag=FALSE, tl.col='black', tl.srt=45, tl.cex=x0.5,  # Reduce text size
         cl.pos="r", cl.cex=0.6, cl.ratio=0.1)


##v2 
corrplot(c, method="color", type="lower", order="hclust",
         col=pal, p.mat=p.cor, sig.level=0.05, insig="blank",
         diag=FALSE, tl.col='black', tl.srt=45, tl.cex=0.5,  # Adjust text size
         cl.pos="r", cl.cex=0.7, cl.ratio=0.4, cl.align="r")  # Move legend right

##Variable correlated (From Oct2024 deadline)

#sea_bottom_temperature correlated to dist_200m, dist_500m, depth (>0.75)
#dist_coast correlated to dist_200m, dist_500m (>0.75)
#dist_200m correlated to depth, dist_500m (>0.75)
#dist_500m correlated to depth (>0.75)
#dist_smnt correlated to dist_smntSML (>0.75)

## Variable to eliminate
#dist_200m
#dist_500m
#dist_smntSML

###Dynamic -----
df<-scores_sel_dyn
df.m<-as.matrix(df[,c("log_sst", 
                      "deep_temp",
                      "deep_sea_water_salinity", 
                      "shallow_sea_water_salinity", 
                      "deep_eastward_sea_water_velocity",
                      "deep_northward_sea_water_velocity",
                      "shallow_eastward_sea_water_velocity", 
                      "shallow_northward_sea_water_velocity",
                      "log_deep_EKE",
                      "log_shallow_EKE", 
                      "log_sea_bottom_temperature", 
                      "sea_surface_height",
                      # "sqrt_dist_coast", 
                      # "sqrt_dist_200m", 
                      # "sqrt_dist_500m", 
                      # "sqrt_dist_1000m",
                      # "depth", 
                      # "log_slope",
                      # "sqrt_dist_smnt", 
                      # "sqrt_dist_smntLRG", 
                      # "sqrt_dist_smntSML", 
                      "log_chloro")])

any(is.na(df.m))  # Check for NA values
any(is.infinite(df.m))  # Check for Inf values

df.m[is.na(df.m)] <- 0
df.m[is.infinite(df.m)] <- 0

p.cor<-cor.mtest(df.m)
pal <- wes_palette("Zissou1", 20, type = "continuous")
corrplot(c, method="color", type="lower", order="hclust",
         col=pal,p.mat=p.cor,sig.level = 0.05,insig = "blank",
         diag=F,tl.col = 'black',tl.srt = 90,tl.cex = .65,cl.pos="r",cl.cex = .7,cl.ratio = .1) ##tl.cex was0.5 


###Static -----
df<-scores_sel_stat
df.m<-as.matrix(df[,c("sqrt_dist_coast",
                      "sqrt_dist_200m",
                      "sqrt_dist_500m",
                      "sqrt_dist_1000m",
                      "depth",
                      "log_slope",
                      "sqrt_dist_smnt",
                      "sqrt_dist_smntLRG",
                      "sqrt_dist_smntSML"#,
                      # "log_chloro"
)])

any(is.na(df.m))  # Check for NA values
any(is.infinite(df.m))  # Check for Inf values

# df.m[is.na(df.m)] <- 0
# df.m[is.infinite(df.m)] <- 0

p.cor<-cor.mtest(df.m)
pal <- wes_palette("Zissou1", 20, type = "continuous")
corrplot(c, method="color", type="lower", order="hclust",
         col=pal,p.mat=p.cor,sig.level = 0.05,insig = "blank",
         diag=F,tl.col = 'black',tl.srt = 90,tl.cex = .65,cl.pos="r",cl.cex = .7,cl.ratio = .1) ##tl.cex was0.5 



## FINAL TABLE FOR MODELLING-----

#add ID and presence
data_PMA_satellite<-as.data.frame(scores_sel[,c("tagid","tagid2","Lat","Long","date","year","month","sen_rep","type",
                                                "log_sst", 
                                                "mid_temp",
                                                "deep_temp", 
                                                "shallow_sea_water_salinity",
                                                "mid_sea_water_salinity",
                                                "deep_sea_water_salinity", 
                                                "shallow_eastward_sea_water_velocity",
                                                "shallow_northward_sea_water_velocity", 
                                                "mid_eastward_sea_water_velocity",
                                                "mid_northward_sea_water_velocity", 
                                                "deep_eastward_sea_water_velocity", 
                                                "deep_northward_sea_water_velocity",
                                                "log_shallow_EKE", 
                                                "log_mid_EKE", 
                                                "log_deep_EKE",
                                                "log_sea_bottom_temperature", 
                                                "sea_surface_height", 
                                                "sqrt_dist_coast",
                                                "sqrt_dist_200m", 
                                                "sqrt_dist_500m", 
                                                "sqrt_dist_1000m", 
                                                "depth", 
                                                "log_slope",
                                                "sqrt_dist_smnt", 
                                                "sqrt_dist_smntSML",
                                                "sqrt_dist_smntLRG",
                                                "log_chloro"#,
                                                # "surface_temp", 
                                                # "shallow_EKE", 
                                                # "mid_EKE",
                                                # "deep_EKE"
)])
data_PMA_satellite$date<-as.Date(data_PMA_satellite$date)
data_PMA_satellite$date <- as.Date(as.numeric(data_PMA_satellite$date), origin = "1970-01-01")



# #write.csv(data_PMA_satellite,file="data_PMA_satellite.csv")
# write.csv(data_PMA_satellite,file=paste0(Results,"/data_PMA_satellite.csv"))
# saveRDS(data_PMA_satellite,file=paste0(Results,"/data_PMA_satellite.rds"))

##
#choose one value that is correlated and run it on in the GAMs

###############################################################
## Variables with Means & SD ----
FFDB_pre<- readRDS(paste0(Results,"/all_tracks_env_data_wMeanSD.rds"))
FFDB<-FFDB_pre
str(FFDB)

## Data cleaning Prep ----
#cleaning NA
#FFDB$SLPd<-FFDB$Slope26N
FFDB<-FFDB[FFDB$depth<=0,]## 
FFDB<-FFDB[complete.cases(FFDB$depth),]
FFDB<-FFDB[FFDB$depth!=0,] 
FFDB$depth<-FFDB$depth*-1
FFDB<-FFDB[FFDB$distcoast_km!=0,] ##updated the variable name 


FFDB$date<-as.Date(FFDB$date)

FFDB$MONTH<-month(as.POSIXlt(FFDB$date, format="%d/%m/%Y"))

FFDB$fMONTH <- factor(FFDB$MONTH, levels = c(2,3,4,5, 6, 7,8,9,10),
                      labels = c("February","March","April","May", "June","July","August", "September","October")) ## change this to basic as date=month

FFDB$shallow_EKE_Mean<-0.5* ((FFDB$shallow_eastward_sea_water_velocity_Mean^2)+(FFDB$shallow_northward_sea_water_velocity_Mean^2))

FFDB$mid_EKE_Mean<-0.5* ((FFDB$mid_eastward_sea_water_velocity_Mean^2)+(FFDB$mid_northward_sea_water_velocity_Mean^2))

FFDB$deep_EKE_Mean<-0.5* ((FFDB$deep_eastward_sea_water_velocity^2)+(FFDB$deep_northward_sea_water_velocity^2))
##variables
tagid<-FFDB$tagid
tagid2<-as.character(FFDB$tagid)
Lat<-FFDB$lat 
Long<-FFDB$lon
# date<-FFDB$date
date<-as.Date(FFDB$date)
year<-year(as.POSIXlt(FFDB$date, format="%d/%m/%Y"))
month<-FFDB$MONTH
sen_rep<-FFDB$sen_rep
ID<-FFDB$ID ##this means itiration from sim
type<-FFDB$type

surface_temp_Mean<-FFDB$surface_temperature_Mean 
mid_temp_Mean<-FFDB$mid_temperature_Mean
deep_temp_Mean<-FFDB$deep_temperature_Mean ##updated from sst to surface temperature 

shallow_sea_water_salinity<-FFDB$shallow_sea_water_salinity
mid_sea_water_salinity<-FFDB$mid_sea_water_salinity
deep_sea_water_salinity<-FFDB$deep_sea_water_salinity

shallow_eastward_sea_water_velocity<-FFDB$shallow_eastward_sea_water_velocity
shallow_northward_sea_water_velocity<-FFDB$shallow_northward_sea_water_velocity
shallow_EKE<-FFDB$shallow_EKE

mid_eastward_sea_water_velocity<-FFDB$mid_eastward_sea_water_velocity
mid_northward_sea_water_velocity<-FFDB$mid_northward_sea_water_velocity
mid_EKE<-FFDB$mid_EKE

deep_eastward_sea_water_velocity<-FFDB$deep_eastward_sea_water_velocity
deep_northward_sea_water_velocity<-FFDB$deep_northward_sea_water_velocity
deep_EKE<-FFDB$deep_EKE

sea_bottom_temperature<-FFDB$sea_bottom_temperature
sea_surface_height<-FFDB$sea_surface_height
chloro<-FFDB$chla