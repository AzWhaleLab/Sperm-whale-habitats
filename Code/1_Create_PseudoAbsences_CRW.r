                              ### Create Pseudo-absences using CRW

##Code to generate pseudo-absences for tracking data for Sperm whales in the Azores using Correlated Random Walks
##Code was originally developed by Hazen et al (2017, 2021) and then updated by Pérez-Jorge et al. (2020). This code then were further updated, developed and adapted for this study. 
##Authors: Mareike D. Duffing Romero, Elliott Hazen, Sergi Pérez-Jorge. 

##Load files -----
##Packages 
# require packages
require(adehabitatLT) 
require(maps)       # for map.where
require(mapdata)    # for worldHires
require(sp); require(maptools); library(geosphere); library(raster)

SDM_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod"; setwd(SDM_HD)
Results<-"~/Projects/SpermWhale/3_SpeciesDistMod/Pseudo-absences/Results"
SP_HD<-"~/Projects/SpermWhale/3_SpeciesDistMod/Sperm_whale_tracking_Sergi"##directory of spatial data

##Load Tracking data 
original_tracks_pre<- readRDS("Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.rds") #as RDS file
original_tracks_pre<- read.csv("Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.csv") ##or as csv

original_tracks=original_tracks_pre #n=13 whales
unique(original_tracks$id) ##mdr added v2
print(length(unique(original_tracks$id))) ##mdr added v2
original_tracks<-as.data.frame(original_tracks) #mdr added

##Tidy up tracking data and remove whale 244386 
tags<-as.data.frame(original_tracks[c('id', 'date', 'lon', 'lat')])
require(dplyr)
tags<-filter(tags, !id=="244386") ##mdr added
length(unique(tags$id)) ##mdr added
detach("package:dplyr", unload = TRUE) ##mdr added
data<-tags

data$date<-as.POSIXct(data$date, format="%d/%m/%Y %H:%M",tz='gmt')
nlev_id<-nlevels(factor(data$id)) ; print(nlev_id)
nlevels_id<-levels(factor(data$id)); print(nlevels_id)

# setup map per tag
#mex.map = map('worldHires', 'Gabon', fill=T, col='transparent', plot=F)
mex.map = map('worldHires', fill=T, col='transparent', plot=F, ylim = c(-12,7))

mex.IDs = sapply(strsplit(mex.map$names, ":"), function(x) x[1])
mex.sp = map2SpatialPolygons(mex.map, IDs=mex.IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))
# mex.sp<-st_as_sf() #maybe MDR

##open shapefile for creating pseudo-absences within a given space range
study_area_PMA_final<- shapefile(paste0(SP_HD, "/", "Shapefiles/study_area_PMA_final_3h_final.shp")) 

##MDR removed code about some shapefiles we don't need
map<-map('worldHires',c('Azores'))
map.id = sapply(strsplit(map$names, ":"), function(x) x[1])
mex.sp = map2SpatialPolygons(map,IDs=map.id, proj4string=CRS("+proj=longlat +datum=WGS84")) #added the ID thing here

points(-5.65,56.7233,col=2,pch=18)
map.axes()
map.scale()

# Creation of an object of class "ltraj"
tr = as.ltraj(cbind(data$lon, data$lat), date=data$date, id=data$id)

library(future)
future::plan(multisession, workers = availableCores() - 2)  #run things in parallel background to help reduce overuse of computer #MDR added

### Run simulations for Pseudo-absences----
## Create pseudo-absences
### sensitivity analysis (40 repetitions)

# sen_rep<-as.integer(1:2) ##change this to 2, but test 1 real quick
sen_rep<-as.integer(1:40) ##change this to 2, but test 1 real quick

sen_rep<-as.data.frame(sen_rep)
nlev_sen_rep<-nlevels(factor(sen_rep$sen_rep)) 
nlevels_sen_rep<-levels(factor(sen_rep$sen_rep)) 

##Mdr added this
require(hms) #V2
StartTime=Sys.time() ; StartTime
print(paste0("Start Time:", StartTime))

for (s in 1:nlev_sen_rep){
  print(paste0("Running sensitivity analysis rep: ", s)) #mdr v2, nlev_sen_rep
  
  for (m in 1:nlev_id){
    
    #select by ID
    tr1 = tr[[m]]
    head(tr1)
    #x11()
    #plot(tr1)
    
    # get indicies to rows of trajectory that are not NA for sampling in simulation
    i.tr1.nona = which(!is.na(tr1[['dist']]) & !is.na(tr1[['rel.angle']]))
    
    data_id<-subset(data, data$id==nlevels_id[m])
    animal<-unique(data_id$id) ##mdr added v2
    
    print(paste0("Running simulation for ID: ", animal))
    ##Calculate angle of the original track 80702
    
    angles<-as.matrix(subset(data_id, select=c("lon", "lat")))
    
    angles_total<-0
    
    angle<-0
    
    for(i in 1) {
      
      angle<-bearing(angles[1,], angles[nrow(data_id),],a=6378137, f=1/298.257223563)
      
      if(angle > 0){angles_total<-angle
      
      }
      else {angles_total<-angle + 360
      
      }}
    
    angle_id <-angles_total
    
    
    ### Calcualte distance of the original track 80702
    
    data_id_SP<-data_id
    
    sp::coordinates(data_id_SP) <- ~lon+lat
    sp::proj4string(data_id_SP) <-CRS("+proj=longlat +datum=WGS84 +no_defs")
    
    coordinates(data_id) <- ~lon+lat
    
    Fdist_id<-Line(coordinates(data_id))
    
    Flines_id <-Lines(Fdist_id, ID=nlevels_id[m])
    
    Length_id<-LinesLength(Flines_id, longlat = FALSE)
    
    # setup data frame for simulation
    n.data = nrow(data_id)
    ##Chose your Number of simulations
    #n.sim = n.tag # or change this to 1000, or n.tag+1000
    # n.sim =10 # or change this to 1000, or n.tag+1000
    # n.sim =200 # or change this to 1000, or n.tag+1000 
    n.sim =100 # or change this to 1000, or n.tag+1000 
   
    
    for (k in 1:n.sim){
      
      print(sprintf("  k'th simulation: %d",k))
      sim = data.frame(x=numeric(n.data), y=numeric(n.data), t=as.POSIXct(rep(NA,n.data)), tagid=rep(data_id$id))     # new location
      
      # create column to keep track of sim iteration number
      sim$iteration = rep(k)
      
      # populate initial location
      sim[1,'x'] = data_id$lon[1]
      sim[1,'y'] = data_id$lat[1]
      sim[1,'t'] = tr1$date[1]
      angle      = tr1[2,'abs.angle'] # starting angle
      dtime      = tr1[1,'date']      # starting time
      
      #debug plot check
      #par(mar=c(3.5, 2, 1, 1)) # margin = c(bottom, left, top, right)
      # map('worldHires', xlim=c(-36,-24), ylim=c(36,41))         # atlantic-centric projection
      map('worldHires', xlim=c(-34,-23), ylim=c(34,43))         # atlantic-centric projection
      map.axes()
      lines(data_id$lon, data_id$lat,col='grey')                          # plot grey lines for original track
      points(sim[1,'x'],sim[1,'y'], col='blue', pch=2, cex=2)     # plot initial point like ltraj symbology
      title(main = sprintf("k'th CRW Simulation for tagid %d", k))
      
      
      # for (j in 2:n.sim){ # j = 2
      for (j in 2:n.data){ # j = 2
        on.land = T # start with assumption on land to force getting new location
        while (on.land == T){
          # get random index of trajectory row, outside NA values
          i = sample(i.tr1.nona, 1)
          
          # get distance and angle from original data
          dist  = tr1[i,'dist']
          angle = angle + tr1[i,'rel.angle']
          
          # calculate a new x and y
          x = sim[j-1,'x'] + dist * cos(angle)
          y = sim[j-1,'y'] + dist * sin(angle)
          
          # check if on land
          pt = SpatialPoints(matrix(c(x,y),nrow=1))
          sp::proj4string(pt) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
          sp::proj4string(study_area_PMA_final) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
          place = over(study_area_PMA_final,pt)
          
          #debug plot check
          if (is.na(place)){
            # assign valid on.water location to sim data frame
            sim[j,'x'] = x
            sim[j,'y'] = y
            sim[j,'t'] = sim[j-1,'t'] + tr1[i,'dt']
            
            # update plot with points
            lines(sim[(j-1):j, c('x','y')], col='black')
            points(x, y, col='blue', pch=20)
            
            # set to F so bumps out of while loop and to next simulated location
            on.land = F
          } else {
            points(x, y, col='red', pch=20)
          }
        } # end while (on.land == T){
        
      } # end for (j in 2:n.sim)
      
      # append full tag simulation to sim.alltags
      if (exists('sim.alltags')){
        sim.alltags = rbind(sim.alltags, sim)
      } else {
        sim.alltags = sim
        
      } # end for (k in 2:n.sim)
      
    }
    
    sim.alltags$ID <- cumsum(!duplicated(sim.alltags[4:5]))
    
    #setwd("C:/R/State-space_models/shapefiles/BPH")# data segments >5km
    ##export table
    
    write.table(sim.alltags, "sim.alltags.csv", sep=",")
    
    alltags<-sim.alltags
    
    #rm(sim.alltags)
    
    sp::coordinates(alltags) <- ~x+y
    sp::proj4string(alltags) <-CRS("+proj=longlat +datum=WGS84 +no_defs")
  
    
    nlev<-nlevels(factor(alltags$ID)) 
    nlevels<-levels(factor(alltags$ID)) 
    
    Fdist <- list()
    
    for (l in 1:nlev){
      
      #select by ID
      hh4<-subset(alltags, alltags$ID==l)
      
      #coordinates(hh4) <- ~x+y
      
      #remove row with less than 1 min
      Fdist[[l]] <-Line(coordinates(hh4))
    }
    
    ## Create Lines
    Flines <- list()
    
    for (l in 1:nlev){
      
      #remove row with less than 1 min
      Flines[[l]] <-Lines(Fdist[l], ID=nlevels[l])
    }
    
    ### Create spatiallines
    
    Flines_SL <- SpatialLines(Flines)
    
    ## sample data: line lengths
    require(sf) ##MDR added this 
    # Convert Flines_SL to an sf object
    Flines_sf <- st_as_sf(Flines_SL)
    
    # Check for duplicate IDs #chat gpt recs
    duplicated_ids <- duplicated(Flines_sf$ID)
    if (any(duplicated_ids)) {
      print("Duplicate IDs found:")
      print(Flines_sf$ID[duplicated_ids])
      Flines_sf$ID <- make.unique(as.character(Flines_sf$ID))
    }
    
    # Calculate lengths
    df <- data.frame(len = st_length(Flines_sf))
    
    # Set unique row names
    rownames(df) <- Flines_sf$ID
    
    ## SpatialLines to SpatialLinesDataFrame
    PMA_lines <- SpatialLinesDataFrame(Flines_SL, data = df)
    
    PMA_lines$ID <- seq.int(nrow(PMA_lines))
    
    ##save shapefile
    #writeOGR(BPH_lines, dsn="C:/R/State-space_models/shapefiles/BPH",layer= "BPH_lines", driver="ESRI Shapefile", overwrite_layer = T)
    
    ##Calculate angle of the simulated tracks
    
    nlev<-nlevels(factor(alltags$ID)) 
    nlevels<-levels(factor(alltags$ID)) 
    
    angles_total_alltags<-0
    
    angle_alltags<-0
    
    #angles<-as.matrix(subset(data_80702, select=c("lon", "lat")))
    
    library(geosphere)
    
    for(l in 1:nlev) {
      
      angles<-subset(alltags, alltags$ID==nlevels[l])
      #coordinates(angles) <- ~x+y
      
      angle_alltags<-bearing(angles[1,], angles[nrow(angles),],a=6378137, f=1/298.257223563)
      
      if(angle_alltags > 0){
        
        angles_total_alltags[l]<-angle_alltags
        
      }
      else {angles_total_alltags[l]<-angle_alltags + 360}
      
    }
    
    angle_sim_id <-as.data.frame(angles_total_alltags)
    
    ###Dataframe with all info on iteration number, angles and distance
    
    ID<-unique(alltags$ID)
    angle_sim<-angle_sim_id$angles_total_alltags
    distance_sim<-df$len
    
    PMA_id<-as.data.frame(cbind(ID, angle_sim,distance_sim))
    
    
    ## Calculate flag value for each track
    ##Flag = 2x((distance_whale-distance_sim)/distance_whale) + (angle_whale-angle_sim)/90
    
    distance_whale<-Length_id
    angle_whale<-angle_id
    
    Flag<-0
    distance_flag<-0
    angle_flag<-0
    
    for(l in 1:nlev) {
      
      distance_flag[l]<-abs((distance_whale-distance_sim[l])/distance_whale) 
      
      angle_flag[l]<-abs((angle_whale-angle_sim[l])/90)
      
      Flag[l]<-(2 * (distance_flag[l]))+angle_flag[l]
      
    }
    
    boxplot(Flag)
    
    PMA_id<-cbind(PMA_id,Flag)
    
    ##Create a column identifying those tracks which flag value are on the upper quantile
    #upper_quantile=0 means that it is on the upper quantile
    #upper_quantile=1 means that it is below the upper quantile
    
    PMA_id$upper_quantile<-0
    
    for(i in 1:nrow(PMA_id)) {
      
      if(PMA_id$Flag[i]<=quantile(PMA_id$Flag,.75)){
        
        PMA_id$upper_quantile[i]=1
      }
      else{(PMA_id$upper_quantile[i]=0)}
    }
    
    
    library(plyr)
    
    count(PMA_id$upper_quantile)
    
    ## identify CRW tracks in the upper quartile
    
    Flag_new <- Flag[Flag<quantile(Flag,.75)] 
    
    boxplot(Flag_new)
    
    ## eliminate lines of the upper quantile in BPH_lines
    
    drops <- which(Flag>quantile(Flag,.75))
    
    PMA_lines_new <- PMA_lines[!PMA_lines$ID %in% drops,] 
    
    ##save shapefile
    terra::writeOGR(PMA_lines_new, dsn="C:/R/Projects/Sperm_whale_tracking/Pseudo_absences/New_tracks_original",layer= "PMA_lines_new", driver="ESRI Shapefile", overwrite_layer = T)
    
    ### eliminate simulations of the upper quantiles in alltags (Spatil Points Data frame)
    
    alltags_new <- alltags[!alltags$ID %in% drops,] 
    
    ##save shapefile

    ##Identify tracks that crossed land (Azores)
  
    sp::proj4string(PMA_lines_new) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
    # sp::proj4string(world) <-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") #don´t need it (mdr)
    
    #ole<-gIntersection(world,  BPH_lines_new, byid=T, id=T,  drop_lower_td=FALSE, unaryUnion_if_byid_false=TRUE, checkValidity=T)
    
    #ole<-gIntersection(world,  PMA_lines_new, byid=T)
    
    #if(is.null(ole)==TRUE){
    
    # tt<-0
    
    #}
    
    #if(is.null(ole)==FALSE){
    
    # tt<-sapply(slot(ole, "lines"), function(x) slot(x, "ID"))
    
    #}
    
    #tt_2<-as.numeric(substr(tt, start=4, stop=7))
    
    #PMA_lines_new_final <- PMA_lines_new[!PMA_lines_new$ID %in% tt_2,] 
    
    ##save shapefile
    # terra::writeOGR(PMA_lines_new_final, dsn=paste(Results,"/New_tracks_original"),layer= paste("PMA_lines_new_final_",nlevels_sen_rep[s],"_",nlevels_id[m], sep=""), driver="ESRI Shapefile", overwrite_layer = F)
    
    select<-sample(PMA_lines_new$ID,size=2)
    
    PMA_lines_new_final_id <- PMA_lines_new[PMA_lines_new$ID %in% select,] 
    
    #save shapefile
    # terra::writeOGR(PMA_lines_new_final_id, dsn=paste(Results, "/New_tracks_original",layer= paste("PMA_lines_new_final_id_",nlevels_sen_rep[s],"_",nlevels_id[m], sep=""), driver="ESRI Shapefile", overwrite_layer = F)
    
    ### Spatial points dataframe
    
    #alltags_new_final <- alltags_new[!alltags_new$ID %in% tt_2,] 
    
    ##save shapefile
    #terra::writeOGR(alltags_new_final, dsn="C:/R/Projects/Sperm_whale_tracking/Pseudo_absences/New_tracks_original",layer= paste("alltags_new_final_",nlevels_sen_rep[s],"_",nlevels_id[m], sep=""), driver="ESRI Shapefile", overwrite_layer = F)
    
    alltags_new_final_id <- alltags_new[alltags_new$ID %in% select,]
    
    ##save shapefile
    #writeOGR(alltags_new_final_id, dsn="C:/R/Projects/Sperm_whale_tracking/Pseudo_absences/New_tracks_original",layer= paste("alltags_new_final_id_",nlevels_sen_rep[s],"_",nlevels_id[m], sep=""), driver="ESRI Shapefile", overwrite_layer = F)
    
    hh5<-as.data.frame(alltags_new_final_id)
    
    hh5$sen_rep<-nlevels_sen_rep[s]
    
    if (exists('hh6')){
      hh6 = rbind(hh6, hh5)
    } else {
      hh6 = hh5
      
    }
    rm(sim.alltags)
  }
  
}

##Save outputs
setwd(paste0(Results, "/", "New_tracks_original"))

write.table(hh6, "sim_CRW_all.csv", sep=",")
require(hms) #V2
EndSave=Sys.time()
RunTimeAll=EndSave - StartTime;
TimeString_Run2=paste0("CRW Modeling Start Time: ",
                       StartTime,
                       ";  End Time after saving files: ",
                       EndSave,
                       ";  Duration ",
                       as_hms(round(RunTimeAll,2))) ##TimeString
print(TimeString_Run2)
# future::plan(future::sequential)  #return to single core, reset the computer

### Map your data ----
library(ggplot2)
library(rnaturalearth)
hh6<-all_clean_pre
hh6<-read.csv(paste0(Results, "/", "New_tracks_original/sim_CRW_all.csv"))
sim_data<-filter(all_clean_pre, type=="CRW") ##simulated data
# sim_data$type="CRW"
# ori_indi=tags
ori_indi<-filter(all_clean_pre, type=="original") ##original track data
# ori_indi$type="Original"
world <- ne_countries(scale = "medium", returnclass = "sf"); class(world)
##remember x=lon, y=lat
lon.range <- c(min(sim_data$x)-0.5,max(sim_data$x)+0.5) 
lat.range <- c(min(sim_data$y)-0.5,max(sim_data$y)+0.5)

lon.range <- c(min(sim_data$lon)-0.5,max(sim_data$lon)+0.5) 
lat.range <- c(min(sim_data$lat)-0.5,max(sim_data$lat)+0.5)

p_sim<-ggplot()+
  # geom_sf(data = world)+ #if we wanted to add world map  
  # geom_polygon(data = study_area_PMA_final)+ #if we wanted to add world map
  # geom_polygon(data = mex.sp)+ 
  # geom_point(data = sim_data, aes(x=x, y=y), color="blue", size=0.8)+
  geom_point(data = sim_data, aes(x=lon, y=lat), color="blue", size=0.8)+
  # geom_path(data = ori_indi, aes(x=lon, y=lat), color="black")+
  geom_point(data = ori_indi, aes(x=lon, y=lat), color="red", size=0.8)+
  
  
  xlim(lon.range) + #
  ylim(lat.range) +
  # ggtitle("Original")+
  labs(subtitle = " ")+
  facet_wrap(~tagid)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text( margin = margin( b = 0.75, t = 0),size=11 ),
    strip.background = element_blank(),
    plot.subtitle = element_text(size=8)
  )
p_sim
