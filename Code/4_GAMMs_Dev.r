                    ## 4. Generalized Additive Mixed Models (GAMMs)
##Script 4. Development of SDMs by applying GAMMs via similar approaches by Hazen et al (2017) and then updated by Pérez-Jorge et al. (2020). 
##Mareike  developed Functionss to help with the model evaluation, variable contribution and create prediction maps. The different functions were based on previous code by Elliott Hazen, Sergi Pérez-Jorge and one mapping code from Maria Inês Pinheiro da Silva.
##Authors: Mareike D. Duffing Romero, Elliott Hazen, Sergi Pérez-Jorge
##Version controls----
##v6 Models for updated data for Nov 2024 and MS 
##v7 some January 2026 updates with a new updated prediction map that includes more months of data and bottomSal, plus new function to create a new map for pred and se at a 3year 9 month facet. 
##Load files ----
# library(MuMIn) ##i dont have this
library(mgcv)
library (MASS)
library(Hmisc) 
library(dismo)
library(dsm)
# library(rgdal)
library(nlme)
# library(gamm4)
# library(biomod2) ##load this package after running the gams
library(mgcViz)
library(dplyr)

rm(list = ls())

##Setting working directories 
SDM_HD<-"~/Documents/Documents - Mareike’s MacBook Pro/PhD Azores/Data/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD) 
InResults<-"~/Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/EnviroResults" 
Results<-"~//Projects/SpermWhale/3_SpeciesDistMod/SDM_Workflow_Code/GAM_Results"  ## change this in work computer
Sergi_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi" ##MDR added and change this on the work computer
library(sf)
Az<-sf::st_read(paste0(Sergi_HD, "/Shapefiles/","Azores.shp"))

## source input functions for evlations, maps, etc
Version="v.1.8.R"# Version
source(paste0(SDM_HD, "/Function_GAM_EvalMaps_", Version))

Ori_tracks_pre <- readRDS(file = paste0(SDM_HD, "/Data/SpermWhale_AllYrsTag_Predicted_6hrSSM_13Wh_v2.rds")) 
Ori_tracks<-filter(Ori_tracks_pre, !id=="244386")
###Vars----
var.dyn <- c("surface_temperature",
             "mid_temp", #"mid_temperature",
             "deep_temp", #"deep_temperature", 
             "deep_sea_water_salinity", "mid_sea_water_salinity", "shallow_sea_water_salinity" , "deep_eastward_sea_water_velocity", "deep_northward_sea_water_velocity",  "shallow_eastward_sea_water_velocity",  "shallow_northward_sea_water_velocity", "mid_eastward_sea_water_velocity", "mid_northward_sea_water_velocity","sea_bottom_temperature",  "sea_surface_height", "chloro", "deep_EKE", "shallow_EKE", "mid_EKE")
var.dyn_Mean<-paste0(var.dyn, "_Mean")
var.dyn_SD<-paste0(var.dyn, "_SD")
var_stat<-c("dist_coast", "dist_200m", "dist_500m", "dist_1000m", "depth", "slope", "dist_smnt", "dist_smntLRG", "dist_smntSML")
var.d.trans<-c("log_sst", 
               "log_shallow_EKE", 
               "log_mid_EKE", 
               "log_deep_EKE", 
               "log_sea_bottom_temperature",
               "log_chloro") 
var.d.trans.Mean<-paste0(var.d.trans, "_Mean")
var.d.trans.SD<-paste0(var.d.trans, "_SD")
var.st.trans<-c("log_slope",
                "sqrt_dist_coast",
                "sqrt_dist_200m",
                "sqrt_dist_500m",
                "sqrt_dist_1000m",
                "sqrt_dist_smnt",
                "sqrt_dist_smntLRG",
                "sqrt_dist_smntSML"
)
## Data Prep ----
library(dplyr); library(lubridate)
##load data
data_pre<- readRDS(file=paste0(InResults,"/data_PMA_satellite.rds"))

data_All<-data_pre
# data_All <- data_All %>%
#   tibble::rowid_to_column("RowID")
# data_All
data_All<-data_All%>% rename("mid_temp"="mid_temperature",
                             "mid_temp_Mean"="mid_temperature_Mean",
                             "mid_temp_SD"="mid_temperature_SD",
                             "deep_temp"="deep_temperature",
                             "deep_temp_Mean"="deep_temperature_Mean",
                             "deep_temp_SD"="deep_temperature_SD"
) #moved this above

data<-data_All %>%
  dplyr::select(-all_of(c(var.dyn_Mean, var.dyn_SD, var.d.trans.Mean, var.d.trans.SD))) ##had like this data_NoSpatMeanSD

str(data)
data$year<-format(as.Date(data$datetime, tz="GMT"),"%Y")
data$tagid<-as.factor(data$tagid)
data$tagid2<-as.factor(data$tagid2)

data$ftagid<-factor(data$tagid) ##useless if by number 1-12
# data$ftagid<-factor(data$tagid2)
# data$tagid<-factor(data$tagid2)

data$Month<-as.factor(data$month)
data$dum<-1

tags<-data%>%
  group_by(tagid) %>% #, tagid2
  dplyr::summarise(tot_detections=n())
tags

tags_type<-data%>%
  group_by(type, tagid) %>% #, tagid2
  dplyr::summarise(tot_detections=n())
tags_type
tags_iteration<-data%>%
  group_by(type, iteration,tagid) %>% #, tagid2
  dplyr::summarise(tot_detections=n())
tags_iteration

##select only pseudo absence CRW

pseudo_data <- data[which(data$type =="CRW"),]

pseudo_data$presence<-0

pseudo_data_order <- pseudo_data[order(pseudo_data$tagid, pseudo_data$sen_rep),]

### select 1 simulation by tagid
original_data <- data[which(data$type =="original"),]
original_data$presence<-1
names(original_data)

# write.csv(data, file = paste0(InResults,"/data_PMA_sat_enviro.csv"))

###Data assessment for NA's -----
##Here we assess number of NAs, espcially for some columns, such as deep variables that may contain NAs
library(dplyr)
var.dyn <- c("surface_temperature","mid_temp","deep_temp", "deep_sea_water_salinity", "mid_sea_water_salinity", "shallow_sea_water_salinity" , "deep_eastward_sea_water_velocity", "deep_northward_sea_water_velocity",  "shallow_eastward_sea_water_velocity",  "shallow_northward_sea_water_velocity", "mid_eastward_sea_water_velocity", "mid_northward_sea_water_velocity","sea_bottom_temperature",  "sea_surface_height", "chloro",
             # "chloro_lag1Mo", "chloro_lag2Mo",
             "deep_EKE", "shallow_EKE", "mid_EKE")

var.stat<-c("dist_coast", "dist_200m", "dist_500m", "dist_1000m", "depth", "slope", "dist_smnt", "dist_smntLRG", "dist_smntSML")

var.st.interst<-c("depth",
                  "log_slope",
                  "sqrt_dist_coast",
                  "sqrt_dist_200m",
                  "sqrt_dist_500m",
                  "sqrt_dist_1000m",
                  "sqrt_dist_smnt",
                  "sqrt_dist_smntLRG",
                  "sqrt_dist_smntSML"
)
vars<-c(var.dyn, var.st.interst)

df_summary2 <- data %>%
  group_by(type, sen_rep) %>%
  summarise(
    total_na = sum(across(all_of(vars), ~ sum(is.na(.)))),
    total_inf = sum(across(all_of(vars), ~ sum(is.infinite(.)))),
    .groups = "drop"
  )

n_vars <- length(vars)

df_summary3 <- data %>%
  group_by(type, sen_rep) %>%
  summarise(
    n_rows = n(),  # total rows in this group
    total_na = sum(across(all_of(vars), ~ sum(is.na(.)))),
    # total_inf = sum(across(all_of(vars), ~ sum(is.infinite(.)))),
    .groups = "drop"
  ) %>%
  mutate(
    total_cells = n_rows * n_vars,
    pct_na = (total_na / total_cells) * 100,
    # pct_inf = (total_inf / total_cells) * 100
  ) %>%
  dplyr::select(type, sen_rep, total_na, pct_na#, 
                # total_inf, pct_inf
  )

df_summary_by_type <- data %>%
  group_by(type) %>%
  summarise(
    n_rows = n(),
    total_na = sum(across(all_of(vars), ~ sum(is.na(.)))),
    .groups = "drop"
  ) %>%
  mutate(
    total_cells = n_rows * n_vars,
    pct_na = (total_na / total_cells) * 100
  ) %>%
  dplyr::select(type, total_na, pct_na); print(df_summary_by_type)
# type     total_na pct_na
# CRW         85485   1.71
# original     1834   2.97

##By Type and variable
vars2<-c(var.dyn)
df_na_by_var <- data %>%
  group_by(type) %>%
  summarise(
    across(
      all_of(vars2),
      list(
        total_na = ~ sum(is.na(.)),
        pct_na = ~ mean(is.na(.)) * 100
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

df_na_long <- data %>%
  pivot_longer(cols = all_of(vars2), names_to = "variable", values_to = "value") %>%
  group_by(type, variable) %>%
  summarise(
    total_na = sum(is.na(value)),
    pct_na = mean(is.na(value)) * 100,
    .groups = "drop"
  ) %>%
  arrange(type, variable);df_na_long

df_na_long2 <- data %>%
  pivot_longer(cols = all_of(vars2), names_to = "variable", values_to = "value") %>%
  group_by(type, sen_rep,variable) %>%
  summarise(
    total_na = sum(is.na(value)),
    pct_na = mean(is.na(value)) * 100,
    .groups = "drop"
  ) %>%
  arrange(type, variable);df_na_long2


totalrows<-nrow(data)
df_na_long$Percentage_NA_byCol<-(df_na_long$total_na/totalrows)*100 ##perhaps not but lets see

df_na_long_vars <- data %>%
  pivot_longer(cols = all_of(vars2), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>% #type, 
  summarise(
    total_na = sum(is.na(value)),
    pct_na = mean(is.na(value)) * 100,
    .groups = "drop"
  ) %>%
  arrange(variable) #type, 

totalrows<-nrow(data)
df_na_long$Percentage_NA_byCol<-(df_na_long$total_na/totalrows)*100 ##perhaps not but lets see


df_sum_crw_dyn <- data %>%
  group_by(type, sen_rep) %>%
  summarise(
    n_rows = n(),  # total rows in this group
    total_na = sum(across(all_of(vars2), ~ sum(is.na(.)))),
    # total_inf = sum(across(all_of(vars), ~ sum(is.infinite(.)))),
    .groups = "drop"
  ) %>%
  mutate(
    total_cells = n_rows * n_vars,
    pct_na = (total_na / total_cells) * 100,
    # pct_inf = (total_inf / total_cells) * 100
  ) %>%
  dplyr::select(type, sen_rep, total_na, pct_na#, 
                # total_inf, pct_inf
  ); df_sum_crw_dyn

# library(openxlsx)
# write.xlsx(df_na_long, file = paste0(Results, "/", "SummaryTable_TotalNas_byDynamicVar_Type.xlsx"))
# write.xlsx(df_na_long_vars, file = paste0(Results, "/", "SummaryTable_TotalNas_byDynamicVar.xlsx"))
# write.xlsx(df_sum_crw_dyn, file = paste0(Results, "/", "SummaryTable_TotalNas_SenRep_forDynamic.xlsx"))

##Load Prediction Data ----
library(dplyr)
PRED_DATA_pre<- readRDS(file = paste0(SDM_HD, "/SDM_Workflow_Code/PredictionData/pred_data_env_data.rds")) ##New Pred Code July 2025

# PRED_DATA_pre<- readRDS(file = paste0(SDM_HD, "/SDM_Workflow_Code/PredictionData/pred_data_env_data_AllYrs_wide_Ensemble.rds")) ##new data Jan 2026 for more months -->apply this for Fxn_AllMaps_pred

PRED_DATA<-PRED_DATA_pre
PRED_DATA <- PRED_DATA %>%
  filter(lon >= -32 & lon <= -23.5)

dplyr::glimpse(PRED_DATA)
# PRED_DATA <- PRED_DATA %>% mutate(Grid_ID= row_number())
PRED_DATA <- PRED_DATA %>% mutate(ID= row_number())

PRED_DATA$date_final<-as.Date(PRED_DATA$date_final)
PRED_DATA$month<-month(PRED_DATA$date_final)
PRED_DATA$year<-year(PRED_DATA$date_final)

# PRED_DATA$deep_temp<-PRED_DATA$deep_temperature 
PRED_DATA$ftagid<-1 ###HUH? SHOULDN'T WE HAVE ALL THE WHALES?
PRED_DATA$dum<-0 ##changed to 1 but was 0; 0 is turned off
PRED_DATA$Month<-as.factor(PRED_DATA$month)

PRED_DATA<-PRED_DATA%>% rename("chloro"="chla",
                               # "chloro_Mean"="chla_Mean",
                               # "chloro_SD"="chla_SD", 
                               "mid_temp"="mid_temperature",
                               # "mid_temp_Mean"="mid_temperature_Mean",
                               # "mid_temp_SD"="mid_temperature_SD",
                               "deep_temp"="deep_temperature",
                               # "deep_temp_Mean"="deep_temperature_Mean",
                               # "deep_temp_SD"="deep_temperature_SD"
)
colnames(PRED_DATA)
colSums(is.na(PRED_DATA))

Pred2018<-filter(PRED_DATA, year=="2018")
colSums(is.na(Pred2018))
sum(is.na(Pred2018))

Pred2023<-filter(PRED_DATA, year=="2023")
colSums(is.na(Pred2023))
sum(is.na(Pred2023))

##Summarize Nas
# sum_nas<-PRED_DATA%>%
#         group_by(year, Month)%>%
#         dplyr::summarise(
#           total_detections=n(),
#           total_na = sum(across(all_of(vars), ~ sum(is.na(.)))),
#           .groups = "drop")%>%
#         mutate(PercNa=total_na/total_detections*100)
# 
# sum_na_by_var <- PRED_DATA %>%
#   group_by(year) %>%
#   summarise(across(all_of(vars2),
#       list(total_na = ~ sum(is.na(.)),
#             pct_na = ~ mean(is.na(.)) * 100),
#       .names = "{.col}_{.fn}"),
#       .groups = "drop")
# 
# sum_na_long_vars <- PRED_DATA %>%
#   pivot_longer(cols = all_of(vars2), names_to = "variable", values_to = "value") %>%
#   group_by(year, variable) %>% #type,
#   summarise(
#     total_na = sum(is.na(value)),
#     pct_na = mean(is.na(value)) * 100,
#     .groups = "drop"
#   ) %>%
#   arrange(year, variable) #type,
# sum_na_long_vars

# PRED_DATA2<-filter(PRED_DATA, !date_final=="2018-05-01") ##previously
# PRED_DATA2_cl<-PRED_DATA2%>%
#   filter(!is.na(deep_temp))

PRED_DATA_cl<-PRED_DATA%>%
  filter(!is.na(deep_temp))
##Test this: 
# PRED_DATA_cl<-filter(PRED_DATA_cl, !date_final=="2023-01-01")


##Review some issues with our data
bad_summary <- PRED_DATA %>%
  dplyr::summarise(
    dplyr::across(
      everything(),
      ~ sum(!is.finite(.)),
      .names = "bad_{.col}"
    )
  ) %>%
  tidyr::pivot_longer(everything())

bad_summary %>% dplyr::filter(value > 0)

bad_summary2 <- PRED_DATA_cl %>%
  dplyr::summarise(
    dplyr::across(
      everything(),
      ~ sum(!is.finite(.)),
      .names = "bad_{.col}"
    )
  ) %>%
  tidyr::pivot_longer(everything())

bad_summary2 %>% dplyr::filter(value > 0)


## Applying GAM Models ----
#Deep Models without Nas in CRW-----
##Clean NAs from Pseudo-abs -----
# cleaned_data <- data %>%
#   filter(!is.na(surface_temperature))
nrow(pseudo_data)
sum(is.na(pseudo_data$deep_temp)) #12920
sum(is.na(pseudo_data$surface_temperature)) #1235

pseudo_data_cl<-pseudo_data%>%
  filter(!is.na(deep_temp))
sum(is.na(pseudo_data_cl$deep_temp)) 
sum(is.na(pseudo_data_cl$surface_temperature)) 

data2_cl<-rbind(pseudo_data_cl, original_data)

data2_sum<-data2_cl%>%
  group_by(year, Month)%>%
  dplyr::summarise(totaldets=n())
## Data selection by Simulation ----
##select best model
x<-1:40

model_ID<-sample(x,size=1)  ##For randomly selecting a pseudo-abs-presence dataframe
##If you are testing various combinations of variables, one should select the same set of pseudo-absence/presence replicate

model_ID="30" #or 40
print(model_ID)
data_1b<-subset(pseudo_data_cl, pseudo_data_cl$sen_rep==model_ID)

data_1b<-rbind(data_1b,original_data)
data_1b_backup=data_1b


##Modeling -----
###Model 1 Previous best 4d plus deep vars----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1b <- mgcv::gam(presence ~ s(log_sst)+ 
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1b) #4655.337 #no change
par(mfrow= c (1,1))
plot(Models_1)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1b)

# coef(Models_1)

###
ItemName="Models_1b"
ModelName="BestPrevMod_PlusDeep"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1b_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="BestPrevMod_PlusDeep",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1b",
                          ModelName="BestPrevMod_PlusDeep",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  grid = PRED_DATA_cl,
  # grid = PRED_DATA2_cl, ##had this open
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsAll<-Fxn_AllMaps_Pred(df=data2_cl,
                          grid = PRED_DATA_cl,
                          Ori_tracks = Ori_tracks,
                          FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath) 
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 1b Previous best 4d plus deep vars no seamount----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1b <- mgcv::gam(presence ~ s(log_sst)+ 
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1b) #4655.337 #no change
# par(mfrow= c (1,1))
# plot(Models_1b)
# # plot(Models_0, residuals=TRUE, pch=1)
# 
# par(mfrow= c (2,2)) 
# gam.check(Models_1b)

# coef(Models_1)

###
ItemName="Models_1b_r"
ModelName="BestPrevMod_PlusDeep_NoSmnts"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1b_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="BestPrevMod_PlusDeep_NoSmnts",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data) 

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1b_r",
                          ModelName="BestPrevMod_PlusDeep_NoSmnts",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsAll<-Fxn_AllMaps_Pred(df=data2_cl,
                          grid = PRED_DATA_cl,
                          Ori_tracks = Ori_tracks,
                          FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath) 

Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 1c Previous best 4d plus deep vars no seamount, swapslope----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            # 'sqrt_dist_coast',
            # 'depth',
            'log_slope',
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               # 'sqrt_dist_coast',
               # 'depth',
               'log_slope',
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1c <- mgcv::gam(presence ~ s(log_sst)+ 
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         # s(sqrt_dist_coast)+
                         # s(depth)+
                         s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1c)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1c) #4655.337 #no change
par(mfrow= c (1,1))
plot(Models_1b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1c)

# coef(Models_1)

###
ItemName="Models_1c"
ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapSlope"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1c
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
# Mod1b_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapSlope",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data) 

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1c",
                          ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapSlope",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath) 

Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 1d Previous best 4d plus deep vars no seamount, swapsIso1k----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            # 'sqrt_dist_coast',
            # 'depth',
            # 'log_slope',
            'sqrt_dist_1000m',
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               # 'sqrt_dist_coast',
               # 'depth',
               # 'log_slope',
               'sqrt_dist_1000m',
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1d <- mgcv::gam(presence ~ s(log_sst)+ 
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         # s(sqrt_dist_coast)+
                         s(sqrt_dist_1000m)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1d)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1d) #4655.337 #no change
par(mfrow= c (1,1))
plot(Models_1d)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1b)

# coef(Models_1)

###
ItemName="Models_1d"
ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapIso1k"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1d
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1b_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapIso1k",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data) 

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1d",
                          ModelName="BestPrevMod_PlusDeep_NoSmnts_SwapIso1k",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath) 

Gamplots<-Fxn_GAM_plot(FilePath = FilePath)




###Model 1e Previous best 4d plus deep vars no seamount, Iso1k and coast----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope',
            'sqrt_dist_1000m',
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope',
               'sqrt_dist_1000m',
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1e <- mgcv::gam(presence ~ s(log_sst)+ 
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         s(sqrt_dist_1000m)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1e)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1e) #4655.337 #no change
par(mfrow= c (1,1))
# plot(Models_1e)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1e)

# coef(Models_1)

###
ItemName="Models_1e"
ModelName="BestPrevMod_PlusDeep_NoSmnts_Coast_wIso1k"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1e
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1b_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="BestPrevMod_PlusDeep_NoSmnts_Coast_wIso1k",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data) 

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1e",
                          ModelName="BestPrevMod_PlusDeep_NoSmnts_Coast_wIso1k",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath) 

Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 1.1 Previous best 4d plus deep vars all EKES----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1.1b <- mgcv::gam(presence ~ s(log_sst)+ 
                           s(mid_temp)+
                           s(deep_temp)+
                           s(shallow_sea_water_salinity)+
                           s(mid_sea_water_salinity)+
                           s(deep_sea_water_salinity)+
                           s(shallow_eastward_sea_water_velocity)+
                           s(shallow_northward_sea_water_velocity)+
                           s(log_shallow_EKE)+
                           # s(mid_eastward_sea_water_velocity)+
                           s(mid_northward_sea_water_velocity)+
                           s(log_mid_EKE)+
                           s(deep_eastward_sea_water_velocity)+ 
                           s(deep_northward_sea_water_velocity)+
                           s(log_deep_EKE)+  
                           s(log_sea_bottom_temperature)+
                           s(sea_surface_height)+
                           # s(log_chloro)+
                           s(sqrt_dist_coast)+
                           # s(depth)+
                           # s(log_slope)+
                           # s(sqrt_dist_smnt)+
                           s(sqrt_dist_smntSML)+
                           # s(sqrt_dist_smntLRG)+
                           s(ftagid,bs = "re", by=dum)
                         ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1.1b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1.1b)  #4647.25
par(mfrow= c (1,1))
plot(Models_1.1b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1.1b)

# coef(Models_1)

###
ItemName="Models_1.1b"
ModelName="BestPrevMod_PlusDeep_AllEkes"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1.1b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1.1b_Eval=EvalGAM(SDM_HD=SDM_HD,
                     FilePath=FilePath,
                     ModelName="BestPrevMod_PlusDeep_AllEkes",
                     pseudo_data=pseudo_data_cl, ##change this to clean
                     ori_data = original_data) 




###Model 1.2 Previous best 4d plus deep vars no mid EKE----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1.2b <- mgcv::gam(presence ~ s(log_sst)+ 
                           s(mid_temp)+
                           s(deep_temp)+
                           s(shallow_sea_water_salinity)+
                           s(mid_sea_water_salinity)+
                           s(deep_sea_water_salinity)+
                           s(shallow_eastward_sea_water_velocity)+
                           s(shallow_northward_sea_water_velocity)+
                           s(log_shallow_EKE)+
                           # s(mid_eastward_sea_water_velocity)+
                           s(mid_northward_sea_water_velocity)+
                           # s(log_mid_EKE)+
                           s(deep_eastward_sea_water_velocity)+ 
                           s(deep_northward_sea_water_velocity)+
                           s(log_deep_EKE)+  
                           s(log_sea_bottom_temperature)+
                           s(sea_surface_height)+
                           # s(log_chloro)+
                           s(sqrt_dist_coast)+
                           # s(depth)+
                           # s(log_slope)+
                           # s(sqrt_dist_smnt)+
                           s(sqrt_dist_smntSML)+
                           # s(sqrt_dist_smntLRG)+
                           s(ftagid,bs = "re", by=dum)
                         ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_1.2b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1.2b)  #44648.654
par(mfrow= c (1,1))
plot(Models_1)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1)

# coef(Models_1)

###
ItemName="Models_1.2b" #Models_1.1b
ModelName="BestPrevMod_PlusDeep_NoMidEke" #BestPrevMod_PlusDeep_AllEkes
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1.1b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
# Mod1.1.2b_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestPrevMod_PlusDeep_NoMidEke",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data) 


EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1.2b",
                          ModelName="BestPrevMod_PlusDeep_NoMidEke",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

PRED_DATA_cl$ftagid<-as.factor(PRED_DATA_cl$ftagid)
MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)

MapsAll<-Fxn_AllMaps_Pred(df=data2_cl,
                          grid = PRED_DATA_cl,
                          Ori_tracks = Ori_tracks,
                          FilePath=FilePath)

Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 1.2a Previous best 4d plus deep vars no mid EKE, no smnts----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1.2a <- mgcv::gam(presence ~ s(log_sst)+ 
                           s(mid_temp)+
                           s(deep_temp)+
                           s(shallow_sea_water_salinity)+
                           s(mid_sea_water_salinity)+
                           s(deep_sea_water_salinity)+
                           s(shallow_eastward_sea_water_velocity)+
                           s(shallow_northward_sea_water_velocity)+
                           s(log_shallow_EKE)+
                           # s(mid_eastward_sea_water_velocity)+
                           s(mid_northward_sea_water_velocity)+
                           # s(log_mid_EKE)+
                           s(deep_eastward_sea_water_velocity)+ 
                           s(deep_northward_sea_water_velocity)+
                           s(log_deep_EKE)+  
                           s(log_sea_bottom_temperature)+
                           s(sea_surface_height)+
                           # s(log_chloro)+
                           s(sqrt_dist_coast)+
                           # s(depth)+
                           # s(log_slope)+
                           # s(sqrt_dist_smnt)+
                           # s(sqrt_dist_smntSML)+
                           # s(sqrt_dist_smntLRG)+
                           s(ftagid,bs = "re", by=dum)
                         ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_1.2a)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1.2a)  #44648.654
par(mfrow= c (1,1))
plot(Models_1a)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_1a)

# coef(Models_1)

###
ItemName="Models_1.2a" #Models_1.1b
ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts" #BestPrevMod_PlusDeep_AllEkes
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1.2a
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
Mod1.1.2b_Eval=EvalGAM(SDM_HD=SDM_HD,
                       FilePath=FilePath,
                       ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts",
                       pseudo_data=pseudo_data_cl, ##change this to clean
                       ori_data = original_data)


EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1.2a",
                          ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

###Model 1.2b_r Previous best 4d plus deep vars no mid EKE, no smnts, swap dist coast w depth----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            # 'sqrt_dist_coast',
            'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               # 'sqrt_dist_coast',
               'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1.2b <- mgcv::gam(presence ~ s(log_sst)+ 
                           s(mid_temp)+
                           s(deep_temp)+
                           s(shallow_sea_water_salinity)+
                           s(mid_sea_water_salinity)+
                           s(deep_sea_water_salinity)+
                           s(shallow_eastward_sea_water_velocity)+
                           s(shallow_northward_sea_water_velocity)+
                           s(log_shallow_EKE)+
                           # s(mid_eastward_sea_water_velocity)+
                           s(mid_northward_sea_water_velocity)+
                           # s(log_mid_EKE)+
                           s(deep_eastward_sea_water_velocity)+ 
                           s(deep_northward_sea_water_velocity)+
                           s(log_deep_EKE)+  
                           s(log_sea_bottom_temperature)+
                           s(sea_surface_height)+
                           # s(log_chloro)+
                           # s(sqrt_dist_coast)+
                           s(depth)+
                           # s(log_slope)+
                           # s(sqrt_dist_smnt)+
                           # s(sqrt_dist_smntSML)+
                           # s(sqrt_dist_smntLRG)+
                           s(ftagid,bs = "re", by=dum)
                         ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_1.2b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1.2b)  #44648.654
par(mfrow= c (1,1))
# plot(Models_1.2b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
# gam.check(Models_1a)

# coef(Models_1)

###
ItemName="Models_1.2b_r" #Models_1.1b
ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast" #BestPrevMod_PlusDeep_AllEkes
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1.2b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
# Mod1.1.2b_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data)


EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1.2b_r",
                          ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

###Model 1.2c Previous best 4d plus deep vars no mid EKE, no smnts, swap dist coast w depth, no seabottom temp----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            # 'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            # 'sqrt_dist_coast',
            'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            # 'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               # 'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               # 'sqrt_dist_coast',
               'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               # 'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_1.2c <- mgcv::gam(presence ~ s(log_sst)+ 
                           s(mid_temp)+
                           s(deep_temp)+
                           s(shallow_sea_water_salinity)+
                           s(mid_sea_water_salinity)+
                           s(deep_sea_water_salinity)+
                           s(shallow_eastward_sea_water_velocity)+
                           s(shallow_northward_sea_water_velocity)+
                           s(log_shallow_EKE)+
                           # s(mid_eastward_sea_water_velocity)+
                           s(mid_northward_sea_water_velocity)+
                           # s(log_mid_EKE)+
                           s(deep_eastward_sea_water_velocity)+ 
                           s(deep_northward_sea_water_velocity)+
                           s(log_deep_EKE)+  
                           # s(log_sea_bottom_temperature)+
                           s(sea_surface_height)+
                           # s(log_chloro)+
                           # s(sqrt_dist_coast)+
                           s(depth)+
                           # s(log_slope)+
                           # s(sqrt_dist_smnt)+
                           # s(sqrt_dist_smntSML)+
                           # s(sqrt_dist_smntLRG)+
                           s(ftagid,bs = "re", by=dum)
                         ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_1.2c)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_1.2c)  #44648.654
# par(mfrow= c (1,1))
# plot(Models_1.2c)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
# gam.check(Models_1c)

# coef(Models_1)

###
ItemName="Models_1.2c" #Models_1.1b
ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast_NoSbTemp" #BestPrevMod_PlusDeep_AllEkes
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_1.2c
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))
# Mod1.1.2b_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast_NoSbTemp",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data)


EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_1.2c",
                          ModelName="BestPrevMod_PlusDeep_NoMidEke_NoSmnts_SwapDepthwCoast_NoSbTemp",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

###Model 2 Previous best 4d plus deep vars no mid EKE and Interactions ----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_2 <- mgcv::gam(presence ~ s(log_sst, shallow_sea_water_salinity)+ 
                        s(mid_temp, mid_sea_water_salinity)+
                        s(deep_temp, deep_sea_water_salinity)+
                        # s(shallow_sea_water_salinity)+
                        # s(mid_sea_water_salinity)+
                        # s(deep_sea_water_salinity)+
                        s(shallow_eastward_sea_water_velocity)+
                        s(shallow_northward_sea_water_velocity)+
                        s(log_shallow_EKE)+
                        # s(mid_eastward_sea_water_velocity)+
                        s(mid_northward_sea_water_velocity)+
                        # s(log_mid_EKE)+
                        s(deep_eastward_sea_water_velocity)+ 
                        s(deep_northward_sea_water_velocity)+
                        s(log_deep_EKE)+  
                        s(log_sea_bottom_temperature)+
                        s(sea_surface_height)+
                        # s(log_chloro)+
                        s(sqrt_dist_coast)+
                        # s(depth)+
                        # s(log_slope)+
                        # s(sqrt_dist_smnt)+
                        s(sqrt_dist_smntSML)+
                        # s(sqrt_dist_smntLRG)+
                        s(ftagid,bs = "re", by=dum)
                      ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_2)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_2)  #4426.925
par(mfrow= c (1,1))
plot(Models_2)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_2)

# coef(Models_1)

###
ItemName="Models_2"
ModelName="BestMod_Inter_NoMidEKE"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_2
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod3_Eval=EvalGAM(SDM_HD=SDM_HD,
#                        FilePath=FilePath,
#                        ModelName="BestMod_Inter_NoMidEKE",
#                        pseudo_data=pseudo_data_cl, ##change this to clean
#                        ori_data = original_data) 

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_2",
                          ModelName="BestMod_Inter_NoMidEKE",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)
c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

###Model 2a Previous best 4d plus deep vars no mid EKE and Surf Interactions ----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_2a <- mgcv::gam(presence ~ s(log_sst, shallow_sea_water_salinity)+ 
                         s(mid_temp)+
                         # s(mid_temp, mid_sea_water_salinity)+
                         # s(deep_temp)+
                         # s(deep_temp, deep_sea_water_salinity)+
                         # s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_2a)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_2a)  #4606.968 #4426.925 (all interactions)
par(mfrow= c (1,1))
plot(Models_2a)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_2)

# coef(Models_1)

###
ItemName="Models_2a"
ModelName="BestMod_SurfInter_NoMidEKE"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_2a
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod3a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestMod_SurfInter_NoMidEKE",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_2a",
                          ModelName="BestMod_SurfInter_NoMidEKE",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  



###Model 2b Previous best 4d plus deep vars no mid EKE and Mid Interactions ----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_2b <- mgcv::gam(presence ~ s(log_sst)+ #shallow_sea_water_salinity
                         # s(mid_temp)+
                         s(mid_temp, mid_sea_water_salinity)+
                         # s(deep_temp)+
                         # s(deep_temp, deep_sea_water_salinity)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_2b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_2b)  # #4426.925 (all interactions)
par(mfrow= c (1,1))
plot(Models_2)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_2)

# coef(Models_1)

###
ItemName="Models_2b"
ModelName="BestMod_MidInter_NoMidEKE"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_2b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod3b_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestMod_MidInter_NoMidEKE",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_2b",
                          ModelName="BestMod_MidInter_NoMidEKE",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3b_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  


###Model 2c Previous best 4d plus deep vars no Deep EKE and Mid Interactions ----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            'shallow_northward_sea_water_velocity',
            'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity',
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            'deep_northward_sea_water_velocity',
            'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               'shallow_northward_sea_water_velocity',
               'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               'deep_northward_sea_water_velocity',
               'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_2c <- mgcv::gam(presence ~ s(log_sst)+ #shallow_sea_water_salinity
                         s(mid_temp)+
                         # s(mid_temp, mid_sea_water_salinity)+
                         # s(deep_temp)+
                         s(deep_temp, deep_sea_water_salinity)+
                         s(shallow_sea_water_salinity)+
                         s(mid_sea_water_salinity)+
                         # s(deep_sea_water_salinity)+
                         s(shallow_eastward_sea_water_velocity)+
                         s(shallow_northward_sea_water_velocity)+
                         s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         s(deep_eastward_sea_water_velocity)+ 
                         s(deep_northward_sea_water_velocity)+
                         s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_2c)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_2c)  #4576.399 #4426.925 (all interactions)
par(mfrow= c (1,1))
plot(Models_2c)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_2c)

# coef(Models_1)

###
ItemName="Models_2c"
ModelName="BestMod_DeepInter_NoMidEKE"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_2c
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod3c_Eval=EvalGAM(SDM_HD=SDM_HD,
#                   FilePath=FilePath,
#                   ModelName="BestMod_DeepInter_NoMidEKE",
#                   pseudo_data=pseudo_data_cl, ##change this to clean
#                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_2b",
                          ModelName="BestMod_DeepInter_NoMidEKE",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3c_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  


###Model 3 Simplify number 4 Previous best 4d plus deep vars no EKE, no Interactions, no nw, ----
##Note it didnt have dist to the different Isobaths.
var.name<-c('log_sst',
            'mid_temp',
            'deep_temp',
            'shallow_sea_water_salinity',
            'mid_sea_water_salinity',
            'deep_sea_water_salinity',
            'shallow_eastward_sea_water_velocity',
            # 'shallow_northward_sea_water_velocity',
            # 'log_shallow_EKE',
            # 'mid_eastward_sea_water_velocity',    
            'mid_northward_sea_water_velocity', ##removed because its neg
            # 'log_mid_EKE',
            'deep_eastward_sea_water_velocity',
            # 'deep_northward_sea_water_velocity',
            # 'log_deep_EKE',
            'log_sea_bottom_temperature',
            'sea_surface_height',
            # 'log_chloro',
            'sqrt_dist_coast',
            # 'depth',
            # 'log_slope', 
            # 'sqrt_dist_smnt', 
            'sqrt_dist_smntSML',
            # 'sqrt_dist_smntLRG',
            'ftagid', 
            'dum')
envar.names<-c('log_sst',
               'mid_temp',
               'deep_temp',
               'shallow_sea_water_salinity',
               'mid_sea_water_salinity',
               'deep_sea_water_salinity',
               'shallow_eastward_sea_water_velocity',
               # 'shallow_northward_sea_water_velocity',
               # 'log_shallow_EKE',
               # 'mid_eastward_sea_water_velocity',    
               'mid_northward_sea_water_velocity',
               # 'log_mid_EKE',
               'deep_eastward_sea_water_velocity',
               # 'deep_northward_sea_water_velocity',
               # 'log_deep_EKE',
               'log_sea_bottom_temperature',
               'sea_surface_height',
               # 'log_chloro',
               'sqrt_dist_coast',
               # 'depth',
               # 'log_slope', 
               # 'sqrt_dist_smnt', 
               'sqrt_dist_smntSML',
               # 'sqrt_dist_smntLRG',
               'ftagid')

Models_3 <- mgcv::gam(presence ~ s(log_sst)+ 
                        s(mid_temp)+
                        s(deep_temp)+
                        s(shallow_sea_water_salinity)+
                        s(mid_sea_water_salinity)+
                        s(deep_sea_water_salinity)+
                        # s(shallow_sea_water_salinity)+
                        # s(mid_sea_water_salinity)+
                        # s(deep_sea_water_salinity)+
                        s(shallow_eastward_sea_water_velocity)+
                        # s(shallow_northward_sea_water_velocity)+
                        # s(log_shallow_EKE)+
                        # s(mid_eastward_sea_water_velocity)+
                        s(mid_northward_sea_water_velocity)+
                        # s(log_mid_EKE)+
                        s(deep_eastward_sea_water_velocity)+ 
                        # s(deep_northward_sea_water_velocity)+
                        # s(log_deep_EKE)+  
                        s(log_sea_bottom_temperature)+
                        s(sea_surface_height)+
                        # s(log_chloro)+
                        s(sqrt_dist_coast)+
                        # s(depth)+
                        # s(log_slope)+
                        # s(sqrt_dist_smnt)+
                        s(sqrt_dist_smntSML)+
                        # s(sqrt_dist_smntLRG)+
                        s(ftagid,bs = "re", by=dum)
                      ,data=data_1b,family=binomial,REML = TRUE)


summary(Models_3)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_3)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_3)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_3)

# coef(Models_1)

###
ItemName="Models_3"
ModelName="BestMod_Simp5_NoInter_NoEKE_NoNW"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod5_Eval=EvalGAM(SDM_HD=SDM_HD,
                  FilePath=FilePath,
                  ModelName="BestMod_Simp5_NoInter_NoEKE_NoNW",
                  pseudo_data=pseudo_data_cl, ##change this to clean
                  ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_3",
                          ModelName="BestMod_Simp5_NoInter_NoEKE_NoNW",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod4<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                             # sp="Sperm Whale",
                             FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

# Mod3_MapsYear<-Fxn_YrMonth_Maps(Ori_tracks = Ori_tracks,
#                                 sp="Sperm Whale",
#                                 FilePath=FilePath)



###Model 4 Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4 <- mgcv::gam(presence ~ 
                        s(log_sst)+ #, by=Month
                        s(mid_temp)+
                        s(deep_temp)+
                        s(shallow_sea_water_salinity)+
                        # #s(shallow_sea_water_salinity, by=Month)+
                        s(mid_sea_water_salinity)+
                        s(deep_sea_water_salinity)+
                        # s(shallow_eastward_sea_water_velocity)+
                        # s(shallow_northward_sea_water_velocity)+
                        # s(log_shallow_EKE)+
                        # s(mid_eastward_sea_water_velocity)+
                        # s(mid_northward_sea_water_velocity)+
                        # s(log_mid_EKE)+
                        # s(deep_eastward_sea_water_velocity)+ 
                        # s(deep_northward_sea_water_velocity)+
                        # s(log_deep_EKE)+  
                        s(log_sea_bottom_temperature)+
                        s(sea_surface_height)+
                        # s(log_chloro)+
                        s(sqrt_dist_coast)+
                        # s(depth)+
                        # s(log_slope)+
                        # s(sqrt_dist_smnt)+
                        # s(sqrt_dist_smntSML)+
                        # s(sqrt_dist_smntLRG)+
                        s(ftagid,bs = "re", by=dum)
                      ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4"
ModelName="Simple_7_ImpDynVar_Coast"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4",
                          ModelName="Simple_7_ImpDynVar_Coast",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()


###Model 4a1 Simple vars 6a-i of imp swap Coast w Depth ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4a1 <- mgcv::gam(presence ~ 
                          s(log_sst)+ #, by=Month
                          s(mid_temp)+
                          s(deep_temp)+
                          s(shallow_sea_water_salinity)+
                          # #s(shallow_sea_water_salinity, by=Month)+
                          s(mid_sea_water_salinity)+
                          s(deep_sea_water_salinity)+
                          # s(shallow_eastward_sea_water_velocity)+
                          # s(shallow_northward_sea_water_velocity)+
                          # s(log_shallow_EKE)+
                          # s(mid_eastward_sea_water_velocity)+
                          # s(mid_northward_sea_water_velocity)+
                          # s(log_mid_EKE)+
                          # s(deep_eastward_sea_water_velocity)+ 
                          # s(deep_northward_sea_water_velocity)+
                          # s(log_deep_EKE)+  
                          s(log_sea_bottom_temperature)+
                          s(sea_surface_height)+
                          # s(log_chloro)+
                          # s(sqrt_dist_coast)+
                          s(depth)+
                          # s(log_slope)+
                          # s(sqrt_dist_smnt)+
                          # s(sqrt_dist_smntSML)+
                          # s(sqrt_dist_smntLRG)+
                          s(ftagid,bs = "re", by=dum)
                        ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4a1)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4a1)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4a1"
ModelName="Simple_7a1_ImpDynVar_SwapCwDepth"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4a1
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                        FilePath=FilePath,
#                        ModelName="Simple_7a1_ImpDynVar_SwapCwDepth",
#                        pseudo_data=pseudo_data_cl, ##change this to clean
#                        ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4a1",
                          ModelName="Simple_7a1_ImpDynVar_SwapCwDepth",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()


###Model 4a2 Simple vars 6a-i of imp swap Coast w Depth no Seabottom T----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4a2 <- mgcv::gam(presence ~ 
                          s(log_sst)+ #, by=Month
                          s(mid_temp)+
                          s(deep_temp)+
                          s(shallow_sea_water_salinity)+
                          # #s(shallow_sea_water_salinity, by=Month)+
                          s(mid_sea_water_salinity)+
                          s(deep_sea_water_salinity)+
                          # s(shallow_eastward_sea_water_velocity)+
                          # s(shallow_northward_sea_water_velocity)+
                          # s(log_shallow_EKE)+
                          # s(mid_eastward_sea_water_velocity)+
                          # s(mid_northward_sea_water_velocity)+
                          # s(log_mid_EKE)+
                          # s(deep_eastward_sea_water_velocity)+ 
                          # s(deep_northward_sea_water_velocity)+
                          # s(log_deep_EKE)+  
                          # s(log_sea_bottom_temperature)+
                          s(sea_surface_height)+
                          # s(log_chloro)+
                          # s(sqrt_dist_coast)+
                          s(depth)+
                          # s(log_slope)+
                          # s(sqrt_dist_smnt)+
                          # s(sqrt_dist_smntSML)+
                          # s(sqrt_dist_smntLRG)+
                          s(ftagid,bs = "re", by=dum)
                        ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4a2)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4a2)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4a2"
ModelName="Simple_7a2_ImpDynVar_SwapCwDepth_NoSbTemp"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4a2
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7a2_ImpDynVar_SwapCwDepth_NoSbTemp",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4a2",
                          ModelName="Simple_7a2_ImpDynVar_SwapCwDepth_NoSbTemp",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)


###Model 4b Simple vars 6a-i of imp_Swap BottomT w Depth  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4b <- mgcv::gam(presence ~ 
                         s(log_sst)+ #, by=Month
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         # #s(shallow_sea_water_salinity, by=Month)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         # s(shallow_eastward_sea_water_velocity)+
                         # s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         # s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         # s(deep_eastward_sea_water_velocity)+ 
                         # s(deep_northward_sea_water_velocity)+
                         # s(log_deep_EKE)+  
                         # s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         s(sqrt_dist_coast)+
                         s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4b)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4b)

# coef(Models_1)

###
ItemName="Models_4b"
ModelName="Simple_7_ImpDynVar_Coast_Depth"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_Depth",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4b",
                          ModelName="Simple_7_ImpDynVar_Coast_Depth",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()



###Model 4c Simple vars 6a-i of imp_Swap DistCoast w Slope  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  # 'depth',
  'log_slope',
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  'log_slope',
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4c <- mgcv::gam(presence ~ 
                         s(log_sst)+ #, by=Month
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         # #s(shallow_sea_water_salinity, by=Month)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         # s(shallow_eastward_sea_water_velocity)+
                         # s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         # s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         # s(deep_eastward_sea_water_velocity)+ 
                         # s(deep_northward_sea_water_velocity)+
                         # s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         # s(sqrt_dist_coast)+
                         # s(depth)+
                         s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4c)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4c)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4b)

# coef(Models_1)

###
ItemName="Models_4c"
ModelName="Simple_7c_ImpDynVar_Swap2Slope"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4c
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7c_ImpDynVar_Swap2Slope",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4c",
                          ModelName="Simple_7c_ImpDynVar_Swap2Slope",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()




###Model 4d Simple vars 6a-i of imp_Swap DistCoast w Iso1000m  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'sqrt_dist_1000m',
  # 'depth',
  # 'log_slope',
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'sqrt_dist_1000m',
  # 'depth',
  # 'log_slope',
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4d <- mgcv::gam(presence ~ 
                         s(log_sst)+ #, by=Month
                         s(mid_temp)+
                         s(deep_temp)+
                         s(shallow_sea_water_salinity)+
                         # #s(shallow_sea_water_salinity, by=Month)+
                         s(mid_sea_water_salinity)+
                         s(deep_sea_water_salinity)+
                         # s(shallow_eastward_sea_water_velocity)+
                         # s(shallow_northward_sea_water_velocity)+
                         # s(log_shallow_EKE)+
                         # s(mid_eastward_sea_water_velocity)+
                         # s(mid_northward_sea_water_velocity)+
                         # s(log_mid_EKE)+
                         # s(deep_eastward_sea_water_velocity)+ 
                         # s(deep_northward_sea_water_velocity)+
                         # s(log_deep_EKE)+  
                         s(log_sea_bottom_temperature)+
                         s(sea_surface_height)+
                         # s(log_chloro)+
                         # s(sqrt_dist_coast)+
                         s(sqrt_dist_1000m)+
                         # s(depth)+
                         # s(log_slope)+
                         # s(sqrt_dist_smnt)+
                         # s(sqrt_dist_smntSML)+
                         # s(sqrt_dist_smntLRG)+
                         s(ftagid,bs = "re", by=dum)
                       ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4d)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4d)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4b)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4d)

# coef(Models_1)

###
ItemName="Models_4d"
ModelName="Simple_7d_ImpDynVar_Swap2Iso1K"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4d
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7d_ImpDynVar_Swap2Iso1K",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4d",
                          ModelName="Simple_7d_ImpDynVar_Swap2Iso1K",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()



###Model 4_InteractionSurf Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntS <- mgcv::gam(presence ~ 
                             # s(log_sst)+ #, by=Month
                             s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                             s(mid_temp)+
                             s(deep_temp)+
                             # s(shallow_sea_water_salinity)+
                             s(mid_sea_water_salinity)+
                             s(deep_sea_water_salinity)+
                             # s(shallow_eastward_sea_water_velocity)+
                             # s(shallow_northward_sea_water_velocity)+
                             # s(log_shallow_EKE)+
                             # s(mid_eastward_sea_water_velocity)+
                             # s(mid_northward_sea_water_velocity)+
                             # s(log_mid_EKE)+
                             # s(deep_eastward_sea_water_velocity)+ 
                             # s(deep_northward_sea_water_velocity)+
                             # s(log_deep_EKE)+  
                             s(log_sea_bottom_temperature)+
                             s(sea_surface_height)+
                             # s(log_chloro)+
                             s(sqrt_dist_coast)+
                             # s(depth)+
                             # s(log_slope)+
                             # s(sqrt_dist_smnt)+
                             # s(sqrt_dist_smntSML)+
                             # s(sqrt_dist_smntLRG)+
                             s(ftagid,bs = "re", by=dum)
                           ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntS)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntS)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4_IntS"
ModelName="Simple_7_ImpDynVar_Coast_IntS"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntS
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_IntS",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntS",
                          ModelName="Simple_7_ImpDynVar_Coast_IntS",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsAll<-Fxn_AllMaps_Pred(df=data2_cl,
                          grid = PRED_DATA_cl,
                          Ori_tracks = Ori_tracks,
                          FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
# FolderName="Evaluations"
# output_dir <- paste(FilePath,FolderName, sep = "/")  

# plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

depth="shallow"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("log_sst", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nLog Temperature (°C)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)


###Model 4_InteractionSurf Simple vars 6a-i of imp  no seabottom----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntS_noSbT <- mgcv::gam(presence ~ 
                                   # s(log_sst)+ #, by=Month
                                   s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                                   s(mid_temp)+
                                   s(deep_temp)+
                                   # s(shallow_sea_water_salinity)+
                                   s(mid_sea_water_salinity)+
                                   s(deep_sea_water_salinity)+
                                   # s(shallow_eastward_sea_water_velocity)+
                                   # s(shallow_northward_sea_water_velocity)+
                                   # s(log_shallow_EKE)+
                                   # s(mid_eastward_sea_water_velocity)+
                                   # s(mid_northward_sea_water_velocity)+
                                   # s(log_mid_EKE)+
                                   # s(deep_eastward_sea_water_velocity)+ 
                                   # s(deep_northward_sea_water_velocity)+
                                   # s(log_deep_EKE)+  
                                   # s(log_sea_bottom_temperature)+
                                   s(sea_surface_height)+
                                   # s(log_chloro)+
                                   s(sqrt_dist_coast)+
                                   # s(depth)+
                                   # s(log_slope)+
                                   # s(sqrt_dist_smnt)+
                                   # s(sqrt_dist_smntSML)+
                                   # s(sqrt_dist_smntLRG)+
                                   s(ftagid,bs = "re", by=dum)
                                 ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntS_noSbT)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntS_noSbT)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4_IntS_noSbT"
ModelName="Simple_7_ImpDynVar_Coast_IntS_noSbT"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntS_noSbT
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_IntS_noSbT",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntS_noSbT",
                          ModelName="Simple_7_ImpDynVar_Coast_IntS_noSbT",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
# FolderName="Evaluations"
# output_dir <- paste(FilePath,FolderName, sep = "/")  

# plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

depth="shallow"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("log_sst", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nLog Temperature (°C)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)



###Model 4_InteractionSurf Simple vars 6a-i of imp k+5 ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntS_ks <- mgcv::gam(presence ~ 
                                # s(log_sst)+ #, by=Month
                                s(log_sst, shallow_sea_water_salinity, k=5)+ #, by=Month
                                s(mid_temp, k=5)+
                                s(deep_temp, k=5)+
                                # s(shallow_sea_water_salinity)+
                                s(mid_sea_water_salinity, k=5)+
                                s(deep_sea_water_salinity, k=5)+
                                # s(shallow_eastward_sea_water_velocity)+
                                # s(shallow_northward_sea_water_velocity)+
                                # s(log_shallow_EKE)+
                                # s(mid_eastward_sea_water_velocity)+
                                # s(mid_northward_sea_water_velocity)+
                                # s(log_mid_EKE)+
                                # s(deep_eastward_sea_water_velocity)+ 
                                # s(deep_northward_sea_water_velocity)+
                                # s(log_deep_EKE)+  
                                s(log_sea_bottom_temperature, k=5)+
                                s(sea_surface_height, k=5)+
                                # s(log_chloro)+
                                s(sqrt_dist_coast, k=5)+
                                # s(depth)+
                                # s(log_slope)+
                                # s(sqrt_dist_smnt)+
                                # s(sqrt_dist_smntSML)+
                                # s(sqrt_dist_smntLRG)+
                                s(ftagid,bs = "re", by=dum)
                              ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntS_ks)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntS_ks)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4_IntS_ks)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4_IntS_ks)

# coef(Models_1)

###
ItemName="Models_4_IntS_k5s"
ModelName="Simple_7_ImpDynVar_Coast_IntS_wK5s"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntS_ks
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                    FilePath=FilePath,
#                    ModelName="Simple_7_ImpDynVar_Coast_IntS_wK5s",
#                    pseudo_data=pseudo_data_cl, ##change this to clean
#                    ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntS_k5s",
                          ModelName="Simple_7_ImpDynVar_Coast_IntS_wK5s",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  # grid = PRED_DATA2_cl, ##what we had
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsAll<-Fxn_AllMaps_Pred(df=data2_cl,
                          grid = PRED_DATA_cl,
                          Ori_tracks = Ori_tracks,
                          FilePath=FilePath)

# MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
#                               FilePath=FilePath)
MapsClean_Mod6g<-Fxn_JointMap_wide(Ori_tracks = Ori_tracks,
                                   FilePath=FilePath)

##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
# FolderName="Evaluations"
# output_dir <- paste(FilePath,FolderName, sep = "/")  

# plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
# print(plot(c6, select = 9))
dev.off()
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

depth="shallow"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("log_sst", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1, #test .5
                 # phi = 10, #test 50, ori 10
                 phi = 20, 
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nLog Temperature (°C)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)

##Fix output
range_log_sst <- range(gam$model$log_sst, na.rm = TRUE)
range_salinity <- range(gam$model$shallow_sea_water_salinity, na.rm = TRUE)

range_log_sst
range_salinity

pred_vals <- predict(gam, type = "response")
range_response <- range(pred_vals, na.rm = TRUE)
range_response

##test two 
# depth="shallow"
# plot_name <- paste0(depth,"_interaction_smooth_plot", "_right_v2.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5
# int_p <- vis.gam(gam,
#                  view = c("log_sst", "shallow_sea_water_salinity"),
#                  type = "response",
#                  n.grid = 50,
#                  too.far = 0.1,
#                  phi = 25,          # A bit higher angle to reduce exaggeration
#                  theta = 40,        # Keep rotation
#                  # zlim = c(0, 1),    # Fix z-axis scale
#                  ticktype = "detailed",
#                  main = "Surface",
#                  xlab = "\n\nLog Temperature (°C)",
#                  ylab = "\n\nSalinity (PSU)",
#                  zlab = "\n\nResponse",
#                  cex.lab = 1.2,
#                  cex.axis = 0.9,
#                  cex.main = 1.5,
#                  font.main = 2
# )
# int_p
# print(int_p)

# library(itsadug)
# plot_smooth(gam, view = c("log_sst", "shallow_sea_water_salinity"), plot.type = "contour")
# 
# pvisgam(gam, view = c("log_sst", "shallow_sea_water_salinity"), cond = list(C = 5), plot.type = "persp")


###Model 4_InteractionSurf_a Simple vars 6a-i of imp swap C with Depth ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntS_a <- mgcv::gam(presence ~ 
                               # s(log_sst)+ #, by=Month
                               s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                               s(mid_temp)+
                               s(deep_temp)+
                               # s(shallow_sea_water_salinity)+
                               s(mid_sea_water_salinity)+
                               s(deep_sea_water_salinity)+
                               # s(shallow_eastward_sea_water_velocity)+
                               # s(shallow_northward_sea_water_velocity)+
                               # s(log_shallow_EKE)+
                               # s(mid_eastward_sea_water_velocity)+
                               # s(mid_northward_sea_water_velocity)+
                               # s(log_mid_EKE)+
                               # s(deep_eastward_sea_water_velocity)+ 
                               # s(deep_northward_sea_water_velocity)+
                               # s(log_deep_EKE)+  
                               s(log_sea_bottom_temperature)+
                               s(sea_surface_height)+
                               # s(log_chloro)+
                               # s(sqrt_dist_coast)+
                               s(depth)+
                               # s(log_slope)+
                               # s(sqrt_dist_smnt)+
                               # s(sqrt_dist_smntSML)+
                               # s(sqrt_dist_smntLRG)+
                               s(ftagid,bs = "re", by=dum)
                             ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntS_a)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntS_a)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
# plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
# gam.check(Models_4)

# coef(Models_1)

###
ItemName="Models_4_IntS_a"
ModelName="Simple_7_ImpDynVar_Coast_IntS_SwapCwDepth"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntS_a
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                    FilePath=FilePath,
#                    ModelName="Simple_7_ImpDynVar_Coast_IntS_SwapCwDepth",
#                    pseudo_data=pseudo_data_cl, ##change this to clean
#                    ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntS_a",
                          ModelName="Simple_7_ImpDynVar_Coast_IntS_SwapCwDepth",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
# FolderName="Evaluations"
# output_dir <- paste(FilePath,FolderName, sep = "/")  

# plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

depth="shallow"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("log_sst", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nLog Temperature (°C)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)


###Model 4_InteractionSurf_b Simple vars 6a-i of imp  swap C with Depth and no seabott Temp----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  # 'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  # 'sqrt_dist_coast',
  'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntS_b <- mgcv::gam(presence ~ 
                               # s(log_sst)+ #, by=Month
                               s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                               s(mid_temp)+
                               s(deep_temp)+
                               # s(shallow_sea_water_salinity)+
                               s(mid_sea_water_salinity)+
                               s(deep_sea_water_salinity)+
                               # s(shallow_eastward_sea_water_velocity)+
                               # s(shallow_northward_sea_water_velocity)+
                               # s(log_shallow_EKE)+
                               # s(mid_eastward_sea_water_velocity)+
                               # s(mid_northward_sea_water_velocity)+
                               # s(log_mid_EKE)+
                               # s(deep_eastward_sea_water_velocity)+ 
                               # s(deep_northward_sea_water_velocity)+
                               # s(log_deep_EKE)+  
                               # s(log_sea_bottom_temperature)+
                               s(sea_surface_height)+
                               # s(log_chloro)+
                               # s(sqrt_dist_coast)+
                               s(depth)+
                               # s(log_slope)+
                               # s(sqrt_dist_smnt)+
                               # s(sqrt_dist_smntSML)+
                               # s(sqrt_dist_smntLRG)+
                               s(ftagid,bs = "re", by=dum)
                             ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntS_b)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntS_b)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
# plot(Models_4)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
# gam.check(Models_4_IntS_b)

# coef(Models_1)

###
ItemName="Models_4_IntS_b"
ModelName="Simple_7_ImpDynVar_Coast_IntS_b_SwapCwDepth_NoSbTemp"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntS_b
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_IntS_b_SwapCwDepth_NoSbTemp",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntS_b",
                          ModelName="Simple_7_ImpDynVar_Coast_IntS_b_SwapCwDepth_NoSbTemp",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
# FolderName="Evaluations"
# output_dir <- paste(FilePath,FolderName, sep = "/")  

# plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

depth="shallow"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("log_sst", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nLog Temperature (°C)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)



###Model 4_InteractionMid Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntM <- mgcv::gam(presence ~ 
                             s(log_sst)+ #, by=Month
                             # s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                             # s(mid_temp)+
                             s(mid_temp, mid_sea_water_salinity)+
                             s(deep_temp)+
                             s(shallow_sea_water_salinity)+
                             # s(mid_sea_water_salinity)+
                             s(deep_sea_water_salinity)+
                             # s(shallow_eastward_sea_water_velocity)+
                             # s(shallow_northward_sea_water_velocity)+
                             # s(log_shallow_EKE)+
                             # s(mid_eastward_sea_water_velocity)+
                             # s(mid_northward_sea_water_velocity)+
                             # s(log_mid_EKE)+
                             # s(deep_eastward_sea_water_velocity)+ 
                             # s(deep_northward_sea_water_velocity)+
                             # s(log_deep_EKE)+  
                             s(log_sea_bottom_temperature)+
                             s(sea_surface_height)+
                             # s(log_chloro)+
                             s(sqrt_dist_coast)+
                             # s(depth)+
                             # s(log_slope)+
                             # s(sqrt_dist_smnt)+
                             # s(sqrt_dist_smntSML)+
                             # s(sqrt_dist_smntLRG)+
                             s(ftagid,bs = "re", by=dum)
                           ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntM)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntM)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
# plot(Models_4_IntM)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4_IntM)

# coef(Models_1)

###
ItemName="Models_4_IntM"
ModelName="Simple_7_ImpDynVar_Coast_IntM"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntM
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_IntM",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntM",
                          ModelName="Simple_7_ImpDynVar_Coast_IntM",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()

library(mgcViz); library(mgcv)

depth="mid"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300)
int_p<-vis.gam(gam, view=c("mid_temp","mid_sea_water_salinity"), type="response", n.grid=50, too.far=0.1, phi=20, theta=215, ticktype = "detailed", #phi 10
               main="Mid depth",
               xlab = "\n\nTemperature (°C)", 
               ylab = "\n\nSalinity (PSU)", 
               zlab = "\n\nResponse",
               mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
               mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
               cex.lab = 1.2, # Axis label size
               cex.axis = 0.9, # Axis tick label size
               font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
               cex.main=1.5 ##title size
               # mgp = c(4, 4, 0)#, # Adjust the space for x, y, z axis labels had 4, 4, 0
) ##i likey phi=20, theta 220, too.far=0.05
int_p
print(int_p)
dev.off()


###Model 4_InteractionDeep Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntD <- mgcv::gam(presence ~ 
                             s(log_sst)+ #, by=Month
                             # s(log_sst, shallow_sea_water_salinity)+ #, by=Month
                             s(mid_temp)+
                             # s(mid_temp, mid_sea_water_salinity)+
                             # s(deep_temp)+
                             s(deep_temp, deep_sea_water_salinity)+
                             s(shallow_sea_water_salinity)+
                             s(mid_sea_water_salinity)+
                             # s(deep_sea_water_salinity)+
                             # s(shallow_eastward_sea_water_velocity)+
                             # s(shallow_northward_sea_water_velocity)+
                             # s(log_shallow_EKE)+
                             # s(mid_eastward_sea_water_velocity)+
                             # s(mid_northward_sea_water_velocity)+
                             # s(log_mid_EKE)+
                             # s(deep_eastward_sea_water_velocity)+ 
                             # s(deep_northward_sea_water_velocity)+
                             # s(log_deep_EKE)+  
                             s(log_sea_bottom_temperature)+
                             s(sea_surface_height)+
                             # s(log_chloro)+
                             s(sqrt_dist_coast)+
                             # s(depth)+
                             # s(log_slope)+
                             # s(sqrt_dist_smnt)+
                             # s(sqrt_dist_smntSML)+
                             # s(sqrt_dist_smntLRG)+
                             s(ftagid,bs = "re", by=dum)
                           ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntD)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntD)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
# plot(Models_4_IntM)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4_IntD)

# coef(Mod# coef(Mod# coef(Models_1)

###
ItemName="Models_4_IntD"
ModelName="Simple_7_ImpDynVar_Coast_IntD"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntD
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
                   FilePath=FilePath,
                   ModelName="Simple_7_ImpDynVar_Coast_IntD",
                   pseudo_data=pseudo_data_cl, ##change this to clean
                   ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntD",
                          ModelName="Simple_7_ImpDynVar_Coast_IntD",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
# tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()

library(mgcViz); library(mgcv)
depth="deep"
plot_name <- paste0(depth,"_interaction_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300)
int_p<-vis.gam(gam, view=c("deep_temp","deep_sea_water_salinity"), type="response", n.grid=50, too.far=0.1, phi=10, theta=40, ticktype = "detailed",
               main="Deep",
               xlab = "\n\nTemperature (°C)", 
               ylab = "\n\nSalinity (PSU)", 
               zlab = "\n\nResponse",
               mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
               mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
               cex.lab = 1.2, # Axis label size
               cex.axis = 0.9, # Axis tick label size
               font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
               cex.main=1.5 ##title size
) ##i likey phi=20, theta 220, too.far=0.05
int_p
print(int_p)
dev.off()



###Model 4a_Interaction Coast and Temp Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntCoTem <- mgcv::gam(presence ~ 
                                 # s(log_sst)+ 
                                 s(mid_temp)+
                                 s(deep_temp)+
                                 s(shallow_sea_water_salinity)+
                                 # #s(shallow_sea_water_salinity, by=Month)+
                                 s(mid_sea_water_salinity)+
                                 s(deep_sea_water_salinity)+
                                 # s(shallow_eastward_sea_water_velocity)+
                                 # s(shallow_northward_sea_water_velocity)+
                                 # s(log_shallow_EKE)+
                                 # s(mid_eastward_sea_water_velocity)+
                                 # s(mid_northward_sea_water_velocity)+
                                 # s(log_mid_EKE)+
                                 # s(deep_eastward_sea_water_velocity)+ 
                                 # s(deep_northward_sea_water_velocity)+
                                 # s(log_deep_EKE)+  
                                 s(log_sea_bottom_temperature)+
                                 s(sea_surface_height)+
                                 # s(log_chloro)+
                                 # s(sqrt_dist_coast)+
                                 s(sqrt_dist_coast, log_sst)+
                                 # s(depth)+
                                 # s(log_slope)+
                                 # s(sqrt_dist_smnt)+
                                 # s(sqrt_dist_smntSML)+
                                 # s(sqrt_dist_smntLRG)+
                                 s(ftagid,bs = "re", by=dum)
                               ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntCoTem)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntCoTem)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4_IntCoTem)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4_IntCoTem)

# coef(Models_1)

###
ItemName="Models_4_IntCoTem"
ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfTemp"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntCoTem
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                    FilePath=FilePath,
#                    ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfTemp",
#                    pseudo_data=pseudo_data_cl, ##change this to clean
#                    ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntCoTem",
                          ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfTemp",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
dev.off()

library(mgcViz); library(mgcv)

depth="shallow"
plot_name <- paste0(depth,"_interactionCoastTemp_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("sqrt_dist_coast", "log_sst"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nSqrt Distance to coast (km)",
                 ylab = "\n\nLog Temperature (°C)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)


###Model 4a_Interaction Coast and Surf Sal Simple vars 6a-i of imp  ----
var.name<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid', 
  'dum')
envar.names<-c(
  'log_sst',
  'mid_temp',
  'deep_temp',
  'shallow_sea_water_salinity',
  'mid_sea_water_salinity',
  'deep_sea_water_salinity',
  # 'shallow_eastward_sea_water_velocity',
  # 'shallow_northward_sea_water_velocity',
  # 'log_shallow_EKE',
  # 'mid_eastward_sea_water_velocity',    
  # 'mid_northward_sea_water_velocity',
  # 'log_mid_EKE',
  # 'deep_eastward_sea_water_velocity',
  # 'deep_northward_sea_water_velocity',
  # 'log_deep_EKE',
  'log_sea_bottom_temperature',
  'sea_surface_height',
  # 'log_chloro',
  'sqrt_dist_coast',
  # 'depth',
  # 'log_slope', 
  # 'sqrt_dist_smnt', 
  # 'sqrt_dist_smntSML',
  # 'sqrt_dist_smntLRG',
  'ftagid')

Models_4_IntCoSal <- mgcv::gam(presence ~ 
                                 # s(log_sst)+ 
                                 s(mid_temp)+
                                 s(deep_temp)+
                                 # s(shallow_sea_water_salinity)+
                                 s(mid_sea_water_salinity)+
                                 s(deep_sea_water_salinity)+
                                 # s(shallow_eastward_sea_water_velocity)+
                                 # s(shallow_northward_sea_water_velocity)+
                                 # s(log_shallow_EKE)+
                                 # s(mid_eastward_sea_water_velocity)+
                                 # s(mid_northward_sea_water_velocity)+
                                 # s(log_mid_EKE)+
                                 # s(deep_eastward_sea_water_velocity)+ 
                                 # s(deep_northward_sea_water_velocity)+
                                 # s(log_deep_EKE)+  
                                 s(log_sea_bottom_temperature)+
                                 s(sea_surface_height)+
                                 # s(log_chloro)+
                                 # s(sqrt_dist_coast)+
                                 s(sqrt_dist_coast, shallow_sea_water_salinity)+
                                 # s(depth)+
                                 # s(log_slope)+
                                 # s(sqrt_dist_smnt)+
                                 # s(sqrt_dist_smntSML)+
                                 # s(sqrt_dist_smntLRG)+
                                 s(ftagid,bs = "re", by=dum)
                               ,data=data_1b,family=binomial,REML = TRUE)

summary(Models_4_IntCoSal)
# concurvity(Models_1b) ##not needed atm, but helps see the bell btwn each variable
# vis.concurvity(Models_1b) ##see correlation values
AIC(Models_4_IntCoSal)  #4648.654; 4426.925 (mod 3 w interaction)
par(mfrow= c (1,1))
plot(Models_4_IntCoSal)
# plot(Models_0, residuals=TRUE, pch=1)

par(mfrow= c (2,2)) 
gam.check(Models_4_IntCoSal)

# coef(Models_1)

###
ItemName="Models_4_IntCoSal"
ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfSal"
FilePath=paste0(Results, "/FinalModel/CleanMods_NoNAs/SimpleVarMods/", ModelName)
ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")

ifelse(!dir.exists(paste(FilePath)), dir.create(paste(FilePath)), "Folder already exists")
# saveRDS(Models_1, file = paste0(FilePath, "/",ModelName, "_Output.rds"))
# gam<-Models_4_IntCoSal
# data_gam<-data_1b
# save(gam,ItemName, ModelName, model_ID,var.name,envar.names,data_gam,FilePath, file = paste0(FilePath, "/",ModelName,"_Output.RData"))

load(paste0(FilePath, "/", ModelName,"_Output.RData"))

# Mod6a_Eval=EvalGAM(SDM_HD=SDM_HD,
#                    FilePath=FilePath,
#                    ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfSal",
#                    pseudo_data=pseudo_data_cl, ##change this to clean
#                    ori_data = original_data)

EnvarPlot<-Fxn_VarImpPlot(ItemName="Models_4_IntCoSal",
                          ModelName="Simple_7_ImpDynVar_Coast_IntCoastSurfSal",
                          FilePath=FilePath,
                          DataType="Normal")

##perhaps use data with less na's
Mod3_Maps<-pred_maps_gam(#df=data, #had data, switched to data2_cl
  df=data2_cl, #had data, switched to data2_cl
  # grid = PRED_DATA,
  # grid = PRED_DATA2,
  # grid = PRED_DATA_cl,
  grid = PRED_DATA2_cl,
  gam = gam,
  Ori_tracks = Ori_tracks,
  sp="Sperm Whale",
  FilePath=FilePath)

MapsClean_Mod6g<-Fxn_JointMap(Ori_tracks = Ori_tracks,
                              # sp="Sperm Whale",
                              FilePath=FilePath)
##Visualize things
Gamplots<-Fxn_GAM_plot(FilePath = FilePath)

c6 <- getViz(gam)
FolderName="Evaluations"
output_dir <- paste(FilePath,FolderName, sep = "/")  

plot_name <- paste0("smooth_plot_", 2, "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
# print(plot(c6, select = 2))
# print(plot(c6, select = 1))
# print(plot(c6, select = 3))
print(plot(c6, select = 7))
dev.off()

library(mgcViz); library(mgcv)

depth="shallow"
plot_name <- paste0(depth,"_interactionCoastSal_smooth_plot", "_right.tiff")
tiff(filename = file.path(output_dir, plot_name), width = 6, height = 6, units = "in", res = 300) ##had h=5 and w=5

int_p <- vis.gam(gam,
                 view = c("sqrt_dist_coast", "shallow_sea_water_salinity"),
                 type = "response",
                 n.grid = 50,
                 too.far = 0.1,
                 phi = 10,
                 theta = 40,
                 ticktype = "detailed",
                 main="Surface",
                 xlab = "\n\nSqrt Distance to coast (km)",
                 ylab = "\n\nSalinity (PSU)",
                 zlab = "\n\nResponse",  # Push z-axis label upward
                 mar = c(5, 5, 3, 3),# Margins: bottom, left, top, right
                 mgp = c(4, 3, 0),   # Adjust label vs. tick spacing
                 cex.lab = 1.2, # Axis label size
                 cex.axis = 0.9, # Axis tick label size #or0.8
                 font.main = 2, ## title font (1=plain, 2=bold, 3=italic)
                 cex.main=1.5 ##title size
)
int_p
print(int_p)

