            ###Functions for running GAMS
#Load packaages
# library(purrr) ##i dont have it
library(tidyverse) ##i dont think i need it
library(Hmisc)
library(raster)
library(dplyr)
library(lubridate)
library(mgcv)
library(stats)
library(dismo)
# library(gbm)  ##not needed, only for boosted regression trees
library(DHARMa)
library(mgcViz)
library(ggplot2)
# library(wesanderson)
library(RColorBrewer)
library(viridis) 
library(biomod2)
library(sf)
# library(geoR)

require(pryr)
loadMemory<-mem_used() 

##Version Control----
##v0.5: update values add comparison model, added mgviz function to print plots, removed some previous plots (as hashtag) because they make the code freeze
## For maps: 
  ##v0.5 and 0.6 updated to bring all the original data
  ##For v0.7 need to: 
  ##a) add an if statement and bring in predictions that have been computed before function. 
  ##b)subset original tracks by year and month
  ##c) add pred standard error dataframe
  ##d) change coloring scheme
  ## Cowplot
  ##Keep in mind same range of aspect ratio and lat/longs.
##v0.8 
  ##Function EvalGAM: add data_gam for data_1 (df you used for running the gam ori and pseudo by rep)
  ##Function_EvalGAM: add code for looking into NAs/inf by environmental variable and keep running model plus save the error messages for the replicates
##0.8 ## as a function for the bar chart https://stackoverflow.com/questions/17273949/clustered-bar-plot-in-r-using-ggplot2

##v0.9: Function EvalGam:Fixed code for the error log, added assign variables to global enviro, fixed plot printing/saving without having R plot freaking out. Added another layer of NAs summary by replicated. Added steps in our code. 
    ##Functions for maps: added Ori_tracks for the dataframe that contains the original tracks
##v1.0: 
  ##Fxn_EvalGam add a section of code to make sure there are no issues when rbinding our pseudo and original data.
  ##Developed new: Fxn_VarImpPlot to have a nice plot of variable importance maybe we can add
  ####Fix factor levels in case for Fxn_PredMaps
##v1.1 Fxn_VarImpPlot: rename column name of var_importance to match spat calc for resolving merging. added a caption for the plot for anything that is spatial calculated to see difference
##v1.2 
  ###Fxn_JointMap new map for our 3 years of pred year. This function may need to change if we have more years of data. Has a section of subsetting the prediction dataframe to match our Ori_tracks Years and months.
  ###added x axis to fxn_VarImp
##v1.3 Add calc AUC for first original model that we feed into the Eval
      ##Calc mean AIC and sd for all sims and saved output
      ## added ModelName for outputs of Evaluations
##v1.4 wide figure for maps Fxn_JointMaps_wide and added a loop inside for simple singular maps
##v1.4 increase size of axis for fxn_gam plots
##V1.5 Create new mapping for function for the months of our interest Feb-Oct 2018-2024 with new prediction data and SE maps Fxn_AllMaps_Pred
##v1.5 create a fxn to keep vars of gam in our pred and remove some infinite/na vals as needed. Fxn_clean_pred_for_gam
##v1.6 added some tweaks to Fxn_AllMaps_Pred including: print number of difference in rows, added animal tracks to plots. (last saved Jan 29, 2026 1415 hrs)
##v1.6Temporarily turn off tracks to show up on the maps of Fxn_AllMaps_Pred Feb 12. Bur suggest for next version to have an if_statement to have this on the map or not.
##v1.7 Change the loop for months/years maps to have a better title under the Fxn_AllMaps_Pred. 
##v1.8 Fxn_ Gam plots: increase size again for the ticks and labels

##Vis Concurvity ----
vis.concurvity <- function(b, type="estimate"){
  cc <- concurvity(b, full=FALSE)[[type]]
  
  diag(cc) <- NA
  cc[lower.tri(cc)]<-NA
  
  layout(matrix(1:2, ncol=2), widths=c(5,1))
  opar <- par(mar=c(5, 6, 5, 0) + 0.1)
  # main plot
  image(z=cc, x=1:ncol(cc), y=1:nrow(cc), ylab="", xlab="",
        axes=FALSE, asp=1, zlim=c(0,1))
  axis(1, at=1:ncol(cc), labels = colnames(cc), las=2)
  axis(2, at=1:nrow(cc), labels = rownames(cc), las=2)
  # legend
  opar <- par(mar=c(5, 0, 4, 3) + 0.1)
  image(t(matrix(rep(seq(0, 1, len=100), 2), ncol=2)),
        x=1:3, y=1:101, zlim=c(0,1), axes=FALSE, xlab="", ylab="")
  axis(4, at=seq(1,101,len=5), labels = round(seq(0,1,len=5),1), las=2)
  par(opar)
} ## end of Fxn_Concurvity

###save Model
##Evaluation Functions----
##Function Evaluation Model----
EvalGAM<-function(SDM_HD,FilePath, ModelName, pseudo_data, ori_data#, 
                  #data_gam ##df you used for running the gam ori and pseudo by rep
                  )
  {
  print("Load Files")
  load(paste0(FilePath, "/", ModelName,"_Output.RData"))
  form<-gam$formula
  Model_ID_ori=model_ID; print("Model ID for pseudo-absence:"); print(Model_ID_ori) #v0.5
  data_1<-data_gam #v0.8
  require(hms)
  StartTime=Sys.time() ; StartTime
  print(paste0("Start Time:", StartTime))
  
  ###Predictions By best model----
  print("Step 1: Predictions by best model") #v0.9 update to add step
  Models_pred<-data.frame(predict(gam,se.fit=T,newdata=data_1,
                                  type="response",backtransform=F))
  
  Models_pred_fit<-Models_pred$fit
  Models_pred_fit
  
  ### CALCULATE C-INDEX For best model ----
  ###Bring in the variables of importance in my best model 
  ###There is a section here for evaluating my data for the model and across simulation.
  cindex_Models_pred<-rcorr.cens(Models_pred_fit,data_1$presence)#0.8  -->AUC
  cindex_Models_pred
  cindex_Models_pred_Ori<-cindex_Models_pred #v1.3
  
  ##Calculate importance of Variables in our model----
  print("Step 2: Calculate Variable Importance") #v0.9 update to add step
  library(biomod2)
  
  e<-bm_VariablesImportance(bm.model = gam, #Models_4d
                            expl.var = data_1[, var.name], ##check it with a 
                            method = "full_rand",
                            nb.rep=1)
  
  e<-e[e$expl.var!="dum",]
  e$norm.imp<-(e$var.imp/sum(e$var.imp))*100
  print(e)
  
  ##predictions
  Models_pred<-data.frame(predict(gam,se.fit=T,newdata=data_1,
                                  type="response",backtransform=F))
  
  Models_pred_fit<-Models_pred$fit
  # bplot<-boxplot(Models_pred_fit)
  # bplot
  ### CALCULATE C-INDEX  By Simulation of our data set######
  ###There is a section here for evaluating my data for the model and across simulation.
  print("Step 3: Calculate C Index by Simulations") #v0.9 update to add step
  cindex_Models_pred<-rcorr.cens(Models_pred_fit,data_1$presence)#0.8 
  cindex_Models_pred
  
  #pseudo_data<-pseudo_data[!pseudo_data$sen_rep==2,]
  ##Run the evaluation for each simulation 
  nlev_id<-nlevels(factor(pseudo_data$sen_rep)) 
  nlevels_id<-levels(factor(pseudo_data$sen_rep)) 
  # l=1
  ##changes to loop below to make sure we can run simulations and skip some vars
  na_summary <- data.frame() #v0.8
  models_list <- list() #v0.8
  # error_log <- data.frame() #v0.8
  error_log <- list() #v0.9
  
  # Initialize storage objects outside the loop #update v0.9
  if (!exists("na_summary")) na_summary <- data.frame()
  if (!exists("na_summary_all")) na_summary_all <- list()
 
  ##Note Removed the previous original loop thats in version 0-0.7 and updated to below
  
  for (l in 1:nlev_id) {
    
    tryCatch({ #v0.8 trycatch issue
      
      data_1 <- subset(pseudo_data, pseudo_data$sen_rep == nlevels_id[l])
      ##Extra code to review columns and maintain same formats (v1.0)
      # Ensure same column order
      data_1 <- data_1[, colnames(original_data)]
      
        ##Ensure all classes match
        # for (col in colnames(data_1)) {
        #   class_target <- class(original_data[[col]])
        #   class_data1  <- class(data_1[[col]])
        #   
        #   # Handle mismatches
        #   if (class_target != class_data1) {
        #     message(paste("Converting column", col, "from", class_data1, "to", class_target))
        #     
        #     # Basic conversion logic — you can expand this if needed
        #     if (class_target == "numeric") data_1[[col]] <- as.numeric(data_1[[col]])
        #     if (class_target == "character") data_1[[col]] <- as.character(data_1[[col]])
        #     if (class_target == "factor") data_1[[col]] <- factor(data_1[[col]], levels = levels(original_data[[col]]))
        #     if (class_target == "POSIXct") data_1[[col]] <- as.POSIXct(data_1[[col]])
        #   }
        # } #v1.0
      
      ##Combine data 
      data_1 <- rbind(data_1, original_data)
      
      ##Double check NAs
          ##Count NAs in environmental variables for this replicate
          na_counts <- sapply(envar.names, function(var) sum(is.na(data_1[[var]])))
          ## Tidy it up into a data.frame
          na_summary_loop <- data.frame(
            sen_rep = nlevels_id[l],
            variable = names(na_counts),
            na_count = as.integer(na_counts)
          )
          ## Append to a master NA summary dataframe
          na_summary <- rbind(na_summary, na_summary_loop)
          
          # ---- NA Summary for entire form model (bad rows) ----
          bad_rows <- which(!complete.cases(data_1[, all.vars(form)]))
          if (length(bad_rows) > 0) {
            na_summary_all[[as.character(nlevels_id[l])]] <- list(
              count = length(bad_rows),
              rows = data_1[bad_rows, ]
            )
            print(paste("Replicate", nlevels_id[l], "has", length(bad_rows), "bad rows with NA/NaN/Inf."))
          } #v0.9
      
      # Fit the GAM
      Models <- mgcv::gam(formula = form, data = data_1, family = binomial, REML = TRUE)
      
      # Store the model
      models_list[[as.character(nlevels_id[l])]] <- Models
      
      # Predictions
      Models_pred <- data.frame(predict(Models, se.fit = TRUE, newdata = data_1,
                                        type = "response", backtransform = FALSE))
      Models_pred_fit <- Models_pred$fit
      
      # C-index
      cindex_Models_pred <- Hmisc::rcorr.cens(Models_pred_fit, data_1$presence)
      cindex <- as.data.frame(cindex_Models_pred)$cindex_Models_pred[1]
      
      if (exists("cindex_final")) {
        cindex_final <- rbind(cindex_final, cindex)
      } else {
        cindex_final <- cindex
      }
      
      # Model summary
      summary_model <- summary(Models)
      
      # AIC
      AIC_val <- AIC(Models)
      if (exists("AIC_final")) {
        AIC_final <- rbind(AIC_final, AIC_val)
      } else {
        AIC_final <- AIC_val
      }
      
      # R-squared
      r.sq <- summary_model$r.sq
      if (exists("r.sq_final")) {
        r.sq_final <- rbind(r.sq_final, r.sq)
      } else {
        r.sq_final <- r.sq
      }
      
      # Smoothing table
      s.table <- summary_model$s.table
      if (exists("s.table_final")) {
        s.table_final <- rbind(s.table_final, s.table)
      } else {
        s.table_final <- s.table
      }
      
      print(paste("Replicate", nlevels_id[l], "processed successfully."))
      
    }, error = function(e) {
      message(paste("Error in replicate", nlevels_id[l], ":", conditionMessage(e)))
      error_log[[length(error_log) + 1]] <- list(
        sen_rep = nlevels_id[l],
        message = conditionMessage(e)
      )##updated v0.9
    })
  } ##end of loop v0.8
  
  ## Save error log
  if (length(error_log) > 0) {
    error_log_df <- do.call(rbind, lapply(error_log, as.data.frame))
    write.csv(error_log_df, file = file.path(FilePath, "gam_error_log.csv"), row.names = FALSE)
    print(paste("Error log saved to", file.path(FilePath, "gam_error_log.csv")))
  } #v0.8 and updated for v0.9
  
  ## Save fitted models
  saveRDS(models_list, file = file.path(FilePath, "fitted_gam_models.rds"))
  print(paste("Fitted models saved to", file.path(FilePath, "fitted_gam_models.rds"))) #v0.8
  
  ##Save Na summaries 
  ## Summarize overall NA issues
  total_issues <- sum(sapply(na_summary_all, function(x) x$count))
  print(paste("Total replicates with NA issues:", length(na_summary_all)))
  print(paste("Total rows with NA/NaN/Inf across all replicates:", total_issues))
  ## Save the by-variable summary
  write.csv(na_summary, file = file.path(FilePath, "na_summary_by_variable.csv"), row.names = FALSE) #v0.8
  print(paste("Saved NA summary by variable to", file.path(FilePath, "na_summary_by_variable.csv")))
  ## Save the detailed NA list with row data
  saveRDS(na_summary_all, file = file.path(FilePath, "na_summary_all_by_senrep.rds"))
  print(paste("Saved NA summary with row data to", file.path(FilePath, "na_summary_all_by_senrep.rds")))
  
  ##Merge dataframes from loop
  AIC_final_pre<-AIC_final #v0.5
  mean_AIC<-mean(AIC_final) #v1.3
  sd_AIC<-sd(AIC_final) #v1.3
  AIC_allSims_sum_df<-as.data.frame(cbind(mean_AIC, sd_AIC)) #v1.3
  
  mean_cindex_final<-mean(cindex_final)
  sd_cindex_final<-sd(cindex_final)
  
  cindex_final_mean_sd<-as.data.frame(cbind(mean_cindex_final,sd_cindex_final))
  mean_r.sq_final<-mean(r.sq_final)
  sd_r.sq_final<-sd(r.sq_final)

  print("Finished running gams by simulation and saving outputs") #v.0.8
  
  ###Summarize replicates for analyzing importance of variables----
  print("Step 3b: Summarizing reps and calc variable importantance by reps") #v0.9 update to add step
  #each is equal to the number of variables on the final model
  nrows<-order(nlevels_id)
  nrows ##mdr added

  env.count<-length(envar.names) ##v0.5
  env.count
  # ID_rep_a<-rep(1:40, each = 4, len = nrow(AIC_final)) 
  ID_rep_a<-rep(1:40, each = env.count, len = nrow(AIC_final)) 
  AIC_final<-cbind(AIC_final,ID_rep_a)
  
  # ID_rep_s<-rep(1:40, each = 8, len = nrow(s.table_final)) 
  ID_rep_s<-rep(1:40, each = env.count, len = nrow(s.table_final)) 
  s.table_final<-cbind(s.table_final,ID_rep_s)
  
  s.table_final_2<-as.data.frame (s.table_final)
  
  ###variables
  ## de las repeticiones ver cuales son las mejores variables por repeticion

  # variable<-rep(c(var.name
  # ),each = 1,len=nrow(s.table_final)) ##row 852 TEST this
  
  variable<-rep(c(envar.names
  ),each = 1,len=nrow(s.table_final)) #right v0.5
  
  s.table_final_2<-cbind(variable,s.table_final_2)

  s.table_final_2$significance <- rep ("NA", NROW(s.table_final))
  
  #s.table_final_2<-cbind(s.table_final,significance)
  
  for(i in (seq(1,NROW(s.table_final_2)))){
    
    if(s.table_final_2$`p-value`[i]<=0.001){
      
      s.table_final_2$significance[i]="SI"}
    
    else{ s.table_final_2$significance[i]="NO"
    }
  }
  
  count<- as.data.frame(table(s.table_final_2$variable,s.table_final_2$significance))
  count_final<-count[count$Var2=="SI",]
  
  # write.csv(count_final,file=paste0(Output_HD, "/", FolderName, "/","count_final.csv"))
  
  ### Select best model to create predictions
  
  AIC_final_DF<-as.data.frame(AIC_final)
  
  AIC_order <- AIC_final_DF[order(AIC_final_DF$V1),]
  
  ### Select data of the best model----
  print("Step 4: Best Model Assessment") #v0.9 update to add step
  # model_ID<-3  ##Model randomly selected to define the best model
  model_ID<-model_ID  ##Model randomly selected to define the best model
  print(model_ID)
  
  x<-1:40
  model_ID<-sample(x,size=1)# model_ID<-8
  print(model_ID)
  data_1<-subset(pseudo_data, pseudo_data$sen_rep==model_ID)
  data_1<-rbind(data_1,original_data)
  model_ID_v2<-model_ID
  
  data_1<-rbind(data_1,original_data)
  
  Models <- mgcv::gam(formula = form, 
                      data=data_1,family=binomial,REML = TRUE)
  
  summary<-summary(Models)#-REML = 312.57
  print(summary)
  print(AIC(Models))
  
  # capture.output(summary, file = paste0(Output_HD, "/", FolderName, "/","summary_PMA_Best_model.txt"))
  
  ##predictions
  Models_pred<-data.frame(predict(Models,se.fit=T,newdata=data_1,
                                  type="response",backtransform=F))
  Models_pred_fit<-Models_pred$fit
  
    ##Run quick model comparison from original pseudo-selection #v0.5
    print("Compare best model with original simulation")
    print(Model_ID_ori)
    
    data_2<-subset(pseudo_data, pseudo_data$sen_rep==Model_ID_ori)
    data_2<-rbind(data_2,original_data)
  
    Models_comp <- mgcv::gam(formula = form, 
                        data=data_2,family=binomial,REML = TRUE)
    summary2<-summary(Models)
    print(summary2)
    print(AIC(Models_comp))
  
    comp<-getViz(Models)
    print(plot(comp), ask = FALSE)
    
  ##Save files part 1: ----
  print("Create Output Folder")
  Output_HD<-FilePath
  FolderName="Evaluations"
  setwd(Output_HD)
  ifelse(!dir.exists(paste(getwd(),"/", FolderName, sep = "")), dir.create(paste(getwd(), "/", FolderName, sep = "")), "Folder already exists")
  
  
  ##need e, bplot, other items
  ##v1.3 adding ModelName
  write.csv(AIC_final_pre,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_AIC_df_pre.csv")) #v0.5
  write.csv(AIC_allSims_sum_df,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_AIC_sum_df.csv")) #v1.3
  write.csv(AIC_final, file=paste0(Output_HD, "/", FolderName, "/",ModelName, "_AIC_Final_AllSimsDF.csv"))
  
  write.csv(e,file=paste0(Output_HD, "/", FolderName, "/", "EnviroVariables_ImportanceDF.csv")) #left as is
  
  write.csv(cindex_final,file=paste0(Output_HD, "/", FolderName, "/",ModelName, "_cindex_final.csv"))
  write.csv(cindex_final_mean_sd,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_cindex_final_mean_sd.csv"))
  write.csv(mean_r.sq_final,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_mean_r.sq_final.csv"))
  write.csv(sd_r.sq_final,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_sd_r.sq_final.csv"))
  write.csv(r.sq_final,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_r.sq_final.csv"))
  
  capture.output(cindex_Models_pred_Ori,file = paste0(Output_HD, "/", FolderName, "/","AUC_OriFirst_model.txt")) #v1.3 for origingal data
  capture.output(summary, file = paste0(Output_HD, "/", FolderName, "/",ModelName,"_summary_PMA_Best_model.txt"))
  write.csv(count_final,file=paste0(Output_HD, "/", FolderName, "/",ModelName,"_EnviroVar_count_final.csv")) ##rename from count final to Enviro v0.7
  
  ###Assign variables to Global env #v0.9
  # Export selected objects to the global environment
  assign("AIC_final_pre", AIC_final_pre, envir = .GlobalEnv)
  assign("AIC_sum_df", AIC_allSims_sum_df, envir = .GlobalEnv) #V1.3
  assign("Models", Models, envir = .GlobalEnv)
  assign("AIC_final", AIC_final, envir = .GlobalEnv)
  assign("cindex_final", cindex_final, envir = .GlobalEnv)
  assign("s.table_final", s.table_final, envir = .GlobalEnv)
  assign("r.sq_final", r.sq_final, envir = .GlobalEnv)
  assign("EnviroCount_final", count_final, envir = .GlobalEnv)
  assign("EnviroVariables_Imp", e, envir = .GlobalEnv)
  assign("models_list", models_list, envir = .GlobalEnv)
  assign("error_log", error_log, envir = .GlobalEnv)
  print("Key results exported to global environment")
  
  
  #### CALCULATE C-INDEX Best Models-----
  cindex_Models_pred<-rcorr.cens(Models_pred_fit,data_1$presence)#0.8 
  capture.output(cindex_Models_pred,file = paste0(Output_HD, "/", FolderName, "/","cindex_PMA_Best_model.txt"))
  
  #### CALCULATE AIC and plot best models-----
  print("Step 4b: Calculate last AIC and make outputs for final Models")#v0.9 update to add step
  AIC<-AIC(Models)
  capture.output(AIC,  file = paste0(Output_HD, "/", FolderName, "/","AIC_PMA_Best_model.txt"))
  
  # Output_HD_EV<-paste(Results, "FinalModel", ModelName,FolderName, sep = "/")
  Output_HD_EV<-paste(FilePath, "Evaluations",sep = "/") #re-updated v0.8
  setwd(Output_HD_EV)
  
  AIC<-AIC(Models)
  
  capture.output(AIC,  file = paste0(Output_HD, "/", FolderName, "/","AIC_PMA_Best_model.txt"))
  
  tryCatch({ #fixed line 376
    jpeg(file.path(FilePath, paste0("gam_check_", ModelName, "_best_model.jpeg")), 
         width = 1200, height = 1000, res = 150)
    
    par(mfrow = c(2, 2), mar = c(5, 4, 2, 1))  # Safe margin setup for plotting
    gam.check(Models)  # This is your final model
    dev.off()
    message("gam.check plot saved successfully.")
    
  }, error = function(e) {
    message(paste("gam.check plot failed:", conditionMessage(e)))
    # Move forward without halting the function
  })
  
  tryCatch({
    print("Reviewing concurvity")
    # Print concurvity summary
    conc <- concurvity(Models, full = TRUE)
    print(conc)
    # Save concurvity plot
    jpeg(file.path(FilePath, paste0("concurvity_", ModelName, "_best_model.jpeg")), 
         width = 1200, height = 1000, res = 150)
    
    par(mfrow = c(1, 1), mar = c(8, 7, 4, 4))
    vis.concurvity(Models)
    dev.off()
    
    message("Concurvity plot saved successfully.")
    
  }, error = function(e) {
    message(paste("Concurvity check failed:", conditionMessage(e)))
  }) #v0.9
  
  
  print("Step 5: Getting visuals")
  d<-getViz(Models) #get visuals from Models
  assign("ModelVisuals", d, envir = .GlobalEnv) #v0.9
  
  setwd(SDM_HD)
  Results=list(AIC_final_DF, summary, Models, d) 

  print("Step 6 Final: Saving Outputs")
  save(Models,Models_comp,summary,summary2,ItemName, ModelName, Model_ID_ori, model_ID_v2,var.name, envar.names,e,d,file = paste0(Output_HD_EV, "/",ModelName,"_Evaluation_Output.RData")) #v0.5
  
  require(hms)
  EndSave=Sys.time()
  RunTimeAll=EndSave - StartTime;
  TimeString_RunAll=paste0("GAM Evaluation Start Time: ",
                           StartTime,
                           ";  End Time after saving files: ",
                           EndSave,
                           ";  Duration ",
                           as_hms(round(RunTimeAll,2))) ##TimeString
  print(TimeString_RunAll)
  
  ##moved this below
  print("Printing Plots of Gams for best model")
  print(plot(d), ask = FALSE)
  
  print(paste0("Done Running Evaluations for:", ModelName)) #update to add modelName v0.9
  return(Results) 
  ##create list of items for return item
} ##end Evaluation Code

library(mgcViz)
##https://eric.univ-lyon2.fr/iec/material/5-MgcViz.pdf
##https://awllee.github.io/sc1/tidyverse/mgcviz_case_study
##https://rdrr.io/cran/mgcViz/man/plot.gamViz.html
##Function GAM Plots ----
Fxn_GAM_plot <- function(FilePath) {
  require(mgcViz)
  print("Load Files")
  load(paste0(FilePath, "/", ModelName,"_Output.RData"))

  c6<-getViz(gam)
  class(c6)
  pl3 <- plot(c6) ;class(pl3)
  
  
  # print(plot(c6, allTerms = TRUE), pages = 2)
  print(plot(c6), ask = FALSE)
  
  par(ask = FALSE)
  ##this works
  pl3 <- pl3 + 
    l_fitLine() + 
    l_ciLine()+
    l_ciPoly(fill = "lightblue", alpha = 0.3) +
    geom_hline(yintercept = 0,  color = "red") +
    l_rug() +# #adds the PA bar at the bottom v0.7 part 2
    labs(y = "Partial effect") + #v1.8 
    theme(
      axis.text = element_text(size = 19),   #v1.4 was 16; v1.8 reupdated   # tick labels 
      axis.title = element_text(size = 21)  #v1.4 was 18; v1.8 reupdated  # axis titles
    )
  pl3$plots
  
  # print(pl3, pages = 3, allTerms=TRUE) ##for all plots
  print(pl3, pages = 3, ask = FALSE) ##for all plots
  
  print("Saving GAM plots")
  FolderName="Evaluations"
  # Save each individual plot
  
  output_dir <- paste(FilePath, FolderName, sep = "/")  # Specify your output directory
  
  for (i in seq_along(pl3$plots)) {
    # Create a filename for each plot
    plot_name <- paste0("smooth_plot_", i, ".tiff")
    
    # Open a graphics device
    tiff(filename = file.path(output_dir, plot_name), width = 6, height = 4, units = "in", res = 300)
    
    # Print the plot to the graphics device
    print(pl3$plots[[i]])
    
    # Close the graphics device
    dev.off()
  } ##Yes
  
  Results=list(plot = pl3)
  return(Results)
} ##end function plot gams

##Function Make Maps ----
##sp is the species
##need to add directories here, and packages requirements
##need to plot original plots

##v0.5 and 0.6 updated to bring all the original data
##For v0.7 need to: 
    ##a) add an if statement and bring in predictions that have been computed before function. 
    ##b)subset original tracks by year and month
    ##c) add pred standard error dataframe
    ##d) change coloring scheme
    ## Cowplot
      ##Keep in mind same range of aspect ratio and lat/longs.

##V0.9 add an item of Ori tracks
##v1.4 created new function Fxn_JointMap_wide for wider lat maps and made tweets to the sizes of plots when merging these
##v1.5 function to make maps with SE for Feb-Oct for all 3 years Fxn_AllMaps_Pred
##

pred_maps_gam <- function(df, grid, gam, sp, FilePath, 
                          Ori_tracks #v0.9 addition 
                          ) {
  require(ggplot2)
  require(dplyr)
  require(lubridate)
  # Load shapefiles and data
  Az <- sf::st_read(paste0(Sergi_HD, "/Shapefiles/", "Azores.shp"))
  # Ori_tracks <- readRDS(file = paste0(SDM_HD, "/Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.rds")) ##Need to update
  Ori_tracks$month<-month(Ori_tracks$date)
  Ori_tracks$year<-year(Ori_tracks$date)
  
  print("Create Output Map Folder")
  Output_HD <- FilePath
  FolderName <- "OutputMaps"
  ifelse(!dir.exists(paste(Output_HD, "/", FolderName, sep = "")), dir.create(paste(Output_HD, "/", FolderName, sep = "")), "Folder already exists")
  mapsdir <- paste(Output_HD, FolderName, sep = "/")
  
  final_pred_df <- data.frame()  # Initialize final dataframe outside the year loop
  
  # Loop through years
  for (y in unique(grid$year)) {
    grid_year <- grid[grid$year == y, ]
    print(y)
    Ori_tracks_year<-Ori_tracks[Ori_tracks$year == y, ] #v0.7
    
    # Initialize yearly dataframe
    m_pred_df <- data.frame()
    
    # Loop through months
    for (m in unique(grid_year$month)) {
      grid_month <- grid_year[grid_year$month == m, ]
      print(m)
      
      Ori_tracks_month<-Ori_tracks_year[Ori_tracks_year$month == m, ] #v0.7
      
      # ##Fix factor levels in case (v1.0)
      # factor_vars <- names(Filter(is.factor, original_data))  # <- or use whatever dataset trained the model
      # for (v in factor_vars) {
      #   if (v %in% colnames(grid_month)) {
      #     grid_month[[v]] <- factor(grid_month[[v]], levels = levels(original_data[[v]]))
      #   }
      # }
      
      # Prediction for the month
      pred <- data.frame(predict(gam, se.fit = TRUE, newdata = grid_month, type = "response")) #original
      # pred <- data.frame(predict(gam, se.fit = TRUE, newdata = grid_month, type = "response", na.action=na.omit)) #updated v1.0 possibly
      
      # Extract lon, lat, and predicted values
      pred_map <- data.frame(
        lon = grid_month$lon,
        lat = grid_month$lat,
        pred = pred$fit,
        pred_se=pred$se.fit, #v0.7 add se
        date = as.Date(grid_month$date)
      )
      pred_map$year <- year(pred_map$date)
      pred_map$month <- month(pred_map$date)
      
      # Append monthly data to yearly dataframe
      m_pred_df <- rbind(m_pred_df, pred_map)
      
      # Generate and save monthly map
      gmap <- ggplot() +
        geom_tile(data = pred_map, aes(x = lon, y = lat, fill = pred)) +
        geom_point(data = Ori_tracks_month, aes(x = lon, y = lat), size = 0.75) + #switched to ori_track_month v0.7
        geom_sf(data = Az) +
        # scale_fill_fermenter(n.breaks = 10, palette = "RdYlBu", limits = c(0, 1), na.value = NA) +
        scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
        guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
        ggtitle(paste(sp, y, m, sep = "-")) +
        xlab("Longitude") + ylab("Latitude") +
        theme_classic()+
        theme(
          legend.text = element_text(size = 8),
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(0.5, "cm")
        )
      
      ggsave(gmap, file = paste0(mapsdir, "/pred_gam_", y, "_", m, ".tiff"), width = 6.88, height = 4.25, bg = "white", compression = "lzw", dpi = 300)
    }
    
    # Append yearly data to final dataframe
    final_pred_df <- rbind(final_pred_df, m_pred_df)
  }
  
  # Save the final dataframe
  saveRDS(final_pred_df, paste(FilePath, "/final_pred_df.RDS", sep = ""))
  final_pred_df<-readRDS(file = paste(FilePath, "/final_pred_df.RDS", sep = ""))
  
  # Make map by Year
  for (y in unique(final_pred_df$year)) {
    pred_year <- final_pred_df[final_pred_df$year == y, ]
    print("Make Maps for yearly prediction:")
    print(y)
    Ori_tracks_year<-Ori_tracks[Ori_tracks$year == y, ] #v0.7
    
    gmap3 <- ggplot() +
      geom_tile(data = pred_year, aes(x = lon, y = lat, fill = pred)) +
      geom_point(data = Ori_tracks_year, aes(x = lon, y = lat), size = 0.75) +
      geom_sf(data = Az) +
      # scale_fill_fermenter(n.breaks = 10, palette = "RdYlBu", limits = c(0, 1), na.value = NA) +
      scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
      guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
      ggtitle(paste(sp, "distribution in",y, sep = " ")) +
      facet_wrap(~month, ncol = 3, nrow = 4) +  # Set 4 rows
      xlab("Longitude") + ylab("Latitude") +
      theme_classic() +
      theme(
        plot.title = element_text(size = 22, face = "bold"),  # Increase title size
        axis.title = element_text(size = 18, face = "bold"),  # Increase axis title size
        axis.text = element_text(size = 14),                  # Increase axis text size
        legend.text = element_text(size = 10),                # Larger legend text
        legend.key.height = unit(1.5, "cm"),                  # Legend height
        legend.key.width = unit(0.75, "cm"),                  # Legend width
        strip.text = element_text(size = 14),                 # Facet label text
        panel.spacing = unit(1, "lines")#,                     # Panel spacing
        # aspect.ratio = 1                                      # Square aspect ratio
      )
  
    # ggsave(gmap3, file = paste0(mapsdir, "/", "ByYear_", y, "_pred_gam.tiff"), width = 12, height = 12, bg = "white", compression = "lzw", dpi = 400)
    ggsave(gmap3, file = paste0(mapsdir, "/", "ByYear_", y, "_pred_gam.tiff"), 
           width = 16, height = 16, bg = "white", compression = "lzw", dpi = 600)
  }
  
  return(list(plot = gmap3, final_DF = final_pred_df))
  # return(list(final_DF = final_pred_df))
}


Fxn_YrMonth_Maps<- function(sp, FilePath, 
                            Ori_tracks) {
  require(ggplot2)
  require(dplyr)
  require(lubridate)
  
  # Load shapefiles and data
  Az <- sf::st_read(paste0(Sergi_HD, "/Shapefiles/", "Azores.shp"))
  # Ori_tracks <- readRDS(file = paste0(SDM_HD, "/Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.rds"))
  Ori_tracks$month<-month(Ori_tracks$date)
  Ori_tracks$year<-year(Ori_tracks$date)
  
  final_pred_df<-readRDS(file = paste(FilePath, "/final_pred_df.RDS", sep = ""))

# Make map by Year
for (y in unique(final_pred_df$year)) {
  pred_year <- final_pred_df[final_pred_df$year == y, ]
  print("Make Maps for yearly prediction:")
  print(y)
  Ori_tracks_year<-Ori_tracks[Ori_tracks$year == y, ] #v0.7
  
  gmap3 <- ggplot() +
    geom_tile(data = pred_year, aes(x = lon, y = lat, fill = pred)) +
    geom_point(data = Ori_tracks_year, aes(x = lon, y = lat), size = 0.75) +
    geom_sf(data = Az) +
    # scale_fill_fermenter(n.breaks = 10, palette = "RdYlBu", limits = c(0, 1), na.value = NA) +
    scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
    guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
    ggtitle(paste(sp, "distribution in",y, sep = " ")) +
    facet_wrap(~month, ncol = 3, nrow = 4) +  # Set 4 rows
    xlab("Longitude") + ylab("Latitude") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22, face = "bold"),  # Increase title size
      axis.title = element_text(size = 18, face = "bold"),  # Increase axis title size
      axis.text = element_text(size = 14),                  # Increase axis text size
      legend.text = element_text(size = 10),                # Larger legend text
      legend.key.height = unit(1.5, "cm"),                  # Legend height
      legend.key.width = unit(0.75, "cm"),                  # Legend width
      strip.text = element_text(size = 14),                 # Facet label text
      panel.spacing = unit(1, "lines")#,                     # Panel spacing
      # aspect.ratio = 1                                      # Square aspect ratio
    )
  
  
  # ggsave(gmap3, file = paste0(mapsdir, "/", "ByYear_", y, "_pred_gam.tiff"), width = 12, height = 12, bg = "white", compression = "lzw", dpi = 400)
  ggsave(gmap3, file = paste0(mapsdir, "/", "ByYear_", y, "_pred_gam.tiff"), 
         width = 16, height = 16, bg = "white", compression = "lzw", dpi = 600)
}

return(list(plot = gmap3))
}

add_facet_labels <- function(df, date_col = "date") {
  # Ensure the date column is present
  library(dplyr)
  library(lubridate)
  if (!date_col %in% names(df)) stop("Date column not found in dataframe.")
  
  df <- df %>%
    mutate(
      year = year(.data[[date_col]]),
      month_name = format(.data[[date_col]], "%B"),
      label = paste0(month_name, " ", year)
    )
  
  # Create ordered levels for label
  ordered_levels <- df %>%
    distinct(label, .keep_all = TRUE) %>%
    arrange(.data[[date_col]]) %>%
    pull(label)
  
  df <- df %>%
    mutate(
      label = factor(label, levels = ordered_levels),
      facet_row = as.character(year)
    )
  
  return(df)
}


Fxn_JointMap<- function(#sp, 
                        FilePath, Ori_tracks){ #v1.2
  require(ggplot2)
  require(dplyr)
  require(lubridate)
  require(zoo)
  # require(ggh4x) #could be useful
  require(cowplot)
# Load shapefiles and data
Az <- sf::st_read(paste0(Sergi_HD, "/Shapefiles/", "Azores.shp"))
# Ori_tracks <- readRDS(file = paste0(SDM_HD, "/Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.rds"))
Ori_tracks$month<-month(Ori_tracks$date)
Ori_tracks$year<-year(Ori_tracks$date)
Ori_tracks$YearMonth<-format(Ori_tracks$date, format = "%Y-%m")
Ori_tracks$YearMonth2<-as.yearmon(Ori_tracks$date, format = "%Y-%m") #zoo package

final_pred_df<-readRDS(file = paste(FilePath, "/final_pred_df.RDS", sep = ""))
final_pred_df$YearMonth<-format(final_pred_df$date, format = "%Y-%m")
final_pred_df$YearMonth2<-as.yearmon(final_pred_df$date, format = "%Y-%m") #zoo package

##Subset and create new dataframe to match the years of data
# pred_year<- final_pred_df %>%
#             filter(YearMonth %in% Ori_tracks$YearMonth) #option1 but want unique

# Ensure YearMonth columns are consistent and only keep unique months from ori_tracks
valid_months <- unique(Ori_tracks$YearMonth)

# Subset final_pred_df to only include those YearMonth values
pred_year <- final_pred_df %>%
  filter(YearMonth %in% valid_months)

pred_year <- pred_year %>%
  mutate(
    label = format(date, "%B %Y"),
    label = factor(label, levels = unique(label[order(date)])),
    year_group = case_when(
      year(date) == 2018 ~ "Row 1: 2018",
      year(date) == 2023 ~ "Row 2: 2023",
      year(date) == 2024 & month(date) %in% 2:5 ~ "Row 3: 2024 (Part 1)",
      year(date) == 2024 & month(date) %in% 6:9 ~ "Row 4: 2024 (Part 2)"
    )
  )

Ori_tracks <- Ori_tracks %>%
  mutate(
    label = format(date, "%B %Y"),
    year_group = case_when(
      year(date) == 2018 ~ "Row 1: 2018",
      year(date) == 2023 ~ "Row 2: 2023",
      year(date) == 2024 & month(date) %in% 2:5 ~ "Row 3: 2024 (Part 1)",
      year(date) == 2024 & month(date) %in% 6:9 ~ "Row 4: 2024 (Part 2)"
    )
  )

##Subset Data and create seperate maps 

plot_data_2018 <- pred_year %>%
  filter(year == 2018) %>%
  arrange(date) %>%  # arrange ensures chronological order
  mutate(label = factor(label, levels = unique(label)))

# Same for Ori_tracks
track_data_2018 <- Ori_tracks %>%
  filter(year == 2018) %>%
  mutate(label = factor(label, levels = unique(plot_data_2018$label)))

# Then use these in your plot
plot2018 <- ggplot() +
  geom_tile(data = plot_data_2018, aes(x = lon, y = lat, fill = pred)) +
  geom_point(data = track_data_2018, aes(x = lon, y = lat), size = 0.75) +
  geom_sf(data = Az) +
  scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
  guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
  scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
  facet_wrap(~label, ncol = 4) +  
  xlab("") + ylab("") +
  theme_classic() +
  theme(
    axis.text.y  = element_text(size = 10),
    # axis.text.x = element_blank(),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.75, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

plot2018

##2023
plot_data_2023 <- pred_year %>%
  filter(year == 2023) %>%
  arrange(date) %>%  # arrange ensures chronological order
  mutate(label = factor(label, levels = unique(label)))

track_data_2023 <- Ori_tracks %>%
  filter(year == 2023) %>%
  mutate(label = factor(label, levels = unique(plot_data_2023$label)))

# Then use these in your plot
plot2023 <- ggplot() +
  geom_tile(data = plot_data_2023, aes(x = lon, y = lat, fill = pred)) +
  geom_point(data = track_data_2023, aes(x = lon, y = lat), size = 0.75) +
  geom_sf(data = Az) +
  scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
  guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
  scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
  facet_wrap(~label, ncol = 4) +  
  xlab("") + ylab("") +
  theme_classic() +
  theme(
    axis.text.y  = element_text(size = 10),
    # axis.text.x = element_blank(),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.75, "cm"),
    strip.text = element_text(size = 12, face="bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

plot2023

##2024
plot_data_2024 <- pred_year %>%
  filter(year == 2024) %>%
  arrange(date) %>%  # arrange ensures chronological order
  mutate(label = factor(label, levels = unique(label)))

track_data_2024 <- Ori_tracks %>%
  filter(year == 2024) %>%
  mutate(label = factor(label, levels = unique(plot_data_2024$label)))

# Then use these in your plot
plot2024 <- ggplot() +
  geom_tile(data = plot_data_2024, aes(x = lon, y = lat, fill = pred)) +
  geom_point(data = track_data_2024, aes(x = lon, y = lat), size = 0.75) +
  geom_sf(data = Az) +
  scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
  guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
  scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
  facet_wrap(~label, ncol = 4) +  
  xlab("") + ylab("") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    # axis.text.x = element_text(angle = 35),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.75, "cm"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.spacing = unit(1, "lines"),
    legend.position = "none"
  )

plot2024

# plot2018_f<-plot_grid(plot2018, NULL, NULL, rel_widths = c(1, 0.25, 0.25), ncol = 3)
plot2018_f<-plot_grid(plot2018, NULL, NULL, rel_widths = c(1, 0.425, 0.425), ncol = 3) #tried 0.5 too big
plot2023_f<-plot_grid(plot2023, NULL, rel_widths = c(1, 0.25), ncol = 2)

# plot2018_23_f_pre<-plot_grid(plot2018_f, plot2023_f, rel_heights = c(1, 1), rel_widths = c(0.75, 1.5), ncol=1) #V1
plot2018_23_f_pre<-plot_grid(plot2018_f, plot2023_f, rel_heights = c(1, 1), rel_widths = c(1, 1), ncol=1) #V2

# plot2018_23_f<-plot_grid(NULL, plot2018_23_f_pre, rel_widths = c(0.04,1), ncol=2) #0.01

plot2024_f<-plot_grid(plot2024, rel_widths = c(1), ncol = 1)
# plot2024_f<-plot_grid(plot2024, rel_widths = c(1.5), ncol = 1) #width 1.2


# plot_all<-plot_grid(plot2018_f, plot2023_f, plot2024_f, ncol = 1, rel_widths = c(1,1,1), rel_heights = c(1, 1, 2))
# plot_all<-plot_grid(plot2018_23_f, plot2024_f, ncol = 1, rel_widths = c(1,1), rel_heights = c(2, 2))
plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(1,1.75), rel_heights = c(2, 2)) #1.5 for plot2024

plot_all

# Save it
ggsave(paste0(FilePath,"/OutputMaps/final_map.tiff"), plot_all, width = 15, height = 12, dpi = 600, bg="white", compression = "lzw")

ggsave(paste0(FilePath,"/OutputMaps/final_map_2018-2023.tiff"), plot2018_23_f_pre, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw") #width 14.5

ggsave(paste0(FilePath,"/OutputMaps/final_map_2024.tiff"), plot2024_f, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw")
plot2024_f2<-plot_grid(plot2024, rel_widths = c(1.75), ncol = 1) #width 1.2

# ggsave(paste0(FilePath,"/OutputMaps/final_map_2024_b.tiff"), plot2024_f2, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw")


##in progress 
legend_plot <- ggplot() +
  geom_tile(data = filter(pred_year, year=="2018"), aes(x = lon, y = lat, fill = pred)) +
  # #scale_fill_distiller(palette = "Spectral",limits = c(0, 1),na.value = NA, name = "Prediction") +
  scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA, name = "Probability") +
  guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
  facet_wrap(~label, ncol = 5) +  
  theme_classic() +
  theme(
    legend.direction = "vertical",
    legend.justification = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

legend <- cowplot::get_legend(legend_plot)

final_plot<-plot_grid(plot_all, legend, ncol = 2, rel_widths  = c(1.5, 0.1))
final_plot

return(list(plot = plot_all, legend_plot, final_plot))
}

##v1.4
Fxn_JointMap_wide<- function(#sp, 
  FilePath, Ori_tracks){ #v1.2
  require(ggplot2)
  require(dplyr)
  require(lubridate)
  require(zoo)
  # require(ggh4x) #could be useful
  require(cowplot)
  # Load shapefiles and data
  Az <- sf::st_read(paste0(Sergi_HD, "/Shapefiles/", "Azores.shp"))
  # Ori_tracks <- readRDS(file = paste0(SDM_HD, "/Data/SpermWhale_Locals_AllYrsTag_Predicted_6hrSSM_v2.rds"))
  Ori_tracks$month<-month(Ori_tracks$date)
  Ori_tracks$year<-year(Ori_tracks$date)
  Ori_tracks$YearMonth<-format(Ori_tracks$date, format = "%Y-%m")
  Ori_tracks$YearMonth2<-as.yearmon(Ori_tracks$date, format = "%Y-%m") #zoo package
  
  final_pred_df<-readRDS(file = paste(FilePath, "/final_pred_df.RDS", sep = ""))
  final_pred_df$YearMonth<-format(final_pred_df$date, format = "%Y-%m")
  final_pred_df$YearMonth2<-as.yearmon(final_pred_df$date, format = "%Y-%m") #zoo package
  
  valid_months <- unique(Ori_tracks$YearMonth)
  
  ###Prep and subset data----
  # Subset final_pred_df to only include those YearMonth values
  pred_year <- final_pred_df %>%
    filter(YearMonth %in% valid_months)
  
  pred_year <- pred_year %>%
    mutate(
      label = format(date, "%B %Y"),
      label = factor(label, levels = unique(label[order(date)])),
      year_group = case_when(
        year(date) == 2018 ~ "Row 1: 2018",
        year(date) == 2023 ~ "Row 2: 2023",
        year(date) == 2024 & month(date) %in% 2:5 ~ "Row 3: 2024 (Part 1)",
        year(date) == 2024 & month(date) %in% 6:9 ~ "Row 4: 2024 (Part 2)"
      )
    )
  
  Ori_tracks <- Ori_tracks %>%
    mutate(
      label = format(date, "%B %Y"),
      year_group = case_when(
        year(date) == 2018 ~ "Row 1: 2018",
        year(date) == 2023 ~ "Row 2: 2023",
        year(date) == 2024 & month(date) %in% 2:5 ~ "Row 3: 2024 (Part 1)",
        year(date) == 2024 & month(date) %in% 6:9 ~ "Row 4: 2024 (Part 2)"
      )
    )
  
  
  ###Subset Data and create seperate maps ----
  
  plot_data_2018 <- pred_year %>%
    filter(year == 2018) %>%
    arrange(date) %>%  # arrange ensures chronological order
    mutate(label = factor(label, levels = unique(label)))
  
  # Same for Ori_tracks
  track_data_2018 <- Ori_tracks %>%
    filter(year == 2018) %>%
    mutate(label = factor(label, levels = unique(plot_data_2018$label)))
  
  # Then use these in your plot
  plot2018 <- ggplot() +
    geom_tile(data = plot_data_2018, aes(x = lon, y = lat, fill = pred)) +
    geom_point(data = track_data_2018, aes(x = lon, y = lat), size = 0.75) +
    geom_sf(data = Az) +
    scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
    guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
    scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
    facet_wrap(~label, ncol = 4) +  
    xlab("") + ylab("") +
    theme_classic() +
    theme(
      axis.text.y  = element_text(size = 10),
      # axis.text.x = element_blank(),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.75, "cm"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  plot2018
  
  ##2023
  plot_data_2023 <- pred_year %>%
    filter(year == 2023) %>%
    arrange(date) %>%  # arrange ensures chronological order
    mutate(label = factor(label, levels = unique(label)))
  
  track_data_2023 <- Ori_tracks %>%
    filter(year == 2023) %>%
    mutate(label = factor(label, levels = unique(plot_data_2023$label)))
  
  # Then use these in your plot
  plot2023 <- ggplot() +
    geom_tile(data = plot_data_2023, aes(x = lon, y = lat, fill = pred)) +
    geom_point(data = track_data_2023, aes(x = lon, y = lat), size = 0.75) +
    geom_sf(data = Az) +
    scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
    guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
    scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
    facet_wrap(~label, ncol = 4) +  
    xlab("") + ylab("") +
    theme_classic() +
    theme(
      axis.text.y  = element_text(size = 10),
      # axis.text.x = element_blank(),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.75, "cm"),
      strip.text = element_text(size = 12, face="bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  plot2023
  
  ##2024
  plot_data_2024 <- pred_year %>%
    filter(year == 2024) %>%
    arrange(date) %>%  # arrange ensures chronological order
    mutate(label = factor(label, levels = unique(label)))
  
  track_data_2024 <- Ori_tracks %>%
    filter(year == 2024) %>%
    mutate(label = factor(label, levels = unique(plot_data_2024$label)))
  
  # Then use these in your plot
  plot2024 <- ggplot() +
    geom_tile(data = plot_data_2024, aes(x = lon, y = lat, fill = pred)) +
    geom_point(data = track_data_2024, aes(x = lon, y = lat), size = 0.75) +
    geom_sf(data = Az) +
    scale_fill_distiller(n.breaks = 10, palette = "Spectral", limits = c(0, 1), na.value = NA) +
    guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
    scale_y_continuous(breaks = seq(36, 40, by = 1)) +  # <-- set tick marks
    facet_wrap(~label, ncol = 4) +  
    xlab("") + ylab("") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 10),
      # axis.text.x = element_text(angle = 35),
      # axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.75, "cm"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  plot2024
  
  ###Prep and Merge Plots----
  ##Prep to merge
  plot2018_f<-plot_grid(plot2018, NULL, NULL, rel_widths = c(1, 0.5, 0.5), ncol = 3) #v1.4
  plot2023_f<-plot_grid(plot2023, NULL, rel_widths = c(1.05, 0.42), ncol = 2) #v1.4 draft-->vielleicht 1, 0.425; 0.42
  
  
  # plot2018_23_f_pre<-plot_grid(plot2018_f, plot2023_f, rel_heights = c(1, 1), rel_widths = c(0.75, 1.5), ncol=1) #V1
  plot2018_23_f_pre<-plot_grid(plot2018_f, plot2023_f, rel_heights = c(1, 1), rel_widths = c(1, 1), ncol=1) #V2
  
  # plot2018_23_f<-plot_grid(NULL, plot2018_23_f_pre, rel_widths = c(0.04,1), ncol=2) #0.01
  
  plot2024_f<-plot_grid(plot2024, rel_widths = c(1), ncol = 1)
  # plot2024_f<-plot_grid(plot2024, rel_widths = c(1.5), ncol = 1) #width 1.2
  
  
  # plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(1.2,1.75), rel_heights = c(2, 2)) #1.5 for plot2024 #was rel (1, 1.7)
  
  # plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(1.25,1.25), rel_heights = c(2, 2)) #1.5 for plot2024 #testing 1.4 #or 1.25 and 1.25 ##good but change
  
  plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(1.85,1.25), rel_heights = c(2, 2)) #1.5 for plot2024 #testing 1.4 #or 1.25 and 1.25 #good ish
  # plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(2.05,1.25), rel_heights = c(2, 2)) #1.5 for plot2024 #testing 1.4 #or 1.25 and 1.25 
  # plot_all<-plot_grid(plot2018_23_f_pre, plot2024_f, ncol = 1, rel_widths = c(1.85,1), rel_heights = c(2, 2)) #1.5 for plot2024 #testing 1.4 #or 1.25 and 1.25
  
  plot_all
  
  # Save it
  ggsave(paste0(FilePath,"/OutputMaps/final_map.tiff"), plot_all, width = 15, height = 12, dpi = 600, bg="white", compression = "lzw")
  
  ggsave(paste0(FilePath,"/OutputMaps/final_map_2018-2023.tiff"), plot2018_23_f_pre, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw") #width 14.5
  
  ggsave(paste0(FilePath,"/OutputMaps/final_map_2024.tiff"), plot2024_f, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw")
  
  plot2024_f2<-plot_grid(plot2024, rel_widths = c(1.75), ncol = 1) #width 1.2
  
  ggsave(paste0(FilePath,"/OutputMaps/final_map_2024_b.tiff"), plot2024_f2, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw")
  
  plot2018_23_f2<-plot_grid(plot2018_23_f_pre, rel_widths = c(2.25), ncol=1) #V2 #1.85 #v1.4
  ggsave(paste0(FilePath,"/OutputMaps/final_map_20218-23_b.tiff"), plot2018_23_f2, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw") #changing width to
  
  ###
  legend_plot <- ggplot() +
    # Probability raster
    geom_tile(
      data = filter(pred_year, year == "2018"),aes(x = lon, y = lat, fill = pred)) +
      scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1),na.value = NA, name = "Probability") +
    
    # Animal track points
    geom_point(data = track_data_2024, aes(x = lon, y = lat, shape = "Animal \ntrack"),  # mapped for legend
      size = 0.75, color = "black"  # fixed color, but shape is in legend
    ) +
    scale_shape_manual(
      name = NULL,
      values = c("Animal \ntrack" = 16)  # solid circle
    ) +
    
    guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE), shape = guide_legend(override.aes = list(color = "black", size = 1.5))
    ) +facet_wrap(~label, ncol = 5, scales = "fixed") +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "right",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    )
  
  # Extract just the unified legend if you need it separately
  legend_combined <- cowplot::get_legend(legend_plot) #renamed combined plot to legend plot
  
  final_plot<-plot_grid(plot_all, legend_combined, ncol = 2, rel_widths  = c(1, 0.1))
  # ggsave(paste0(FilePath,"/OutputMaps/final_map_wLegend.tiff"), plot_f, width = 15, height = 6, dpi = 600, bg="white", compression = "lzw")
 
  ###Loops for single maps----
  ##v1.4 Loop for clean simple models
  Output_HD <- FilePath
  FolderName <- "OutputMaps/Simple"
  ifelse(!dir.exists(paste(Output_HD, "/", FolderName, sep = "")), dir.create(paste(Output_HD, "/", FolderName, sep = "")), "Folder already exists")
  mapsdir <- paste(Output_HD, FolderName, sep = "/")
  
  # assign("AIC_final_pre", AIC_final_pre, envir = .GlobalEnv)
  pred_year
  Ori_tracks
  
  for (m in unique(pred_year$YearMonth)) {
    pred_map <- pred_year[pred_year$YearMonth == m, ]
    print(m)
    
    Ori_tracks_month<-Ori_tracks[Ori_tracks$YearMonth == m, ] #v0.7
    Monat<-unique(pred_map$label)
    require(ggtext)
    gmap <- ggplot() +
      geom_tile(data = pred_map, aes(x = lon, y = lat, fill = pred)) +
      geom_point(data = Ori_tracks_month, aes(x = lon, y = lat), size = 0.75) + 
      geom_sf(data = Az) +
      scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
      guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
      ggtitle(paste(Monat)) + ##changed this to label
      xlab("") + ylab("") +
      theme_classic()+
      # theme(plot.title= element_text(size = 14, face="bold", hjust = 0.5),#had 12 
      theme(plot.title= element_text(size = 20, face="bold", hjust = 0.5),#can change back to 16, #had 12 or 14
            # axis.text = element_text(size = 10))+ 
            axis.text = element_text(size = 16))+ #can change back to 12
      theme(legend.position = "blank")
    
    ggsave(gmap, file = paste0(mapsdir, "/pred_gam_", m,"_cl", ".tiff"), width = 6, height = 4.25, bg = "white", compression = "lzw", dpi = 300) #ori #or use Monat
  }
  
  return(list(plot = plot_all, legend_plot, final_plot))
}

##v1.5 
Fxn_AllMaps_Pred <- function(df, ##dataframe with original data
                             grid, ##Prediction grid
                             # sp, 
                             FilePath, 
                             Ori_tracks #v0.9 addition 
) { ##New function under v1.5
  require(ggplot2)
  require(dplyr)
  require(lubridate)
  # Load shapefiles and data
  Az <- sf::st_read(paste0(Sergi_HD, "/Shapefiles/", "Azores.shp"))
  Ori_tracks$month<-month(Ori_tracks$date)
  Ori_tracks$year<-year(Ori_tracks$date)
  Ori_tracks$Monat<-format(Ori_tracks$date, format = "%B")
  
  load(paste0(FilePath, "/", ModelName,"_Output.RData")) ##Avoid extra stuff and loads our model without needing to put it at first #v1.5
  
  print("Create Output Map Folder")
  Output_HD <- FilePath
  FolderName <- "OutputMaps/FinalMaps" ##Changed this to have the finalMaps subfolder for this function.v1.5
  ifelse(!dir.exists(paste(Output_HD, "/", FolderName, sep = "")), dir.create(paste(Output_HD, "/", FolderName, sep = "")), "Folder already exists")
  mapsdir <- paste(Output_HD, FolderName, sep = "/") 
  
  ##Prepping the grid:
  grid_ori=grid
  length_ori<-nrow(grid_ori)#; length_ori<-as.numeric(length_ori)
  print(paste0("Number of rows for original grid: ", length_ori)) #v1.6 addition to double check some things.
  
  grid<-Fxn_clean_pred_for_gam(pred_df = grid,
                               gam = gam)
  
  length_gridcl<-nrow(grid)#; length_ori<-as.numeric(length_ori)
  print(paste0("Number of rows for original grid: ", length_gridcl))
  dif_grids<-length_ori-length_gridcl
  print(paste0("Number of rows dropped: ", dif_grids))
  

  grid_sub <- grid %>%
    filter(
      year %in% c(2018, 2023, 2024),
      month %in% c(2:10)
    ) #option A, but if we upgrade the function we can call the Years and Months, see below
  
  # grid_sub$ftagid<-as.factor(grid_sub$ftagid) ##see if this helps
  
  ##Option B: if we add the Years and Months before hand
  # Years=c(2018, 2023, 2024)
  # Months=2:10
  # grid_sub2 <- grid %>%
  #           filter(year %in% Years,
  #                  month %in% Months)
  
  final_pred_df <- data.frame()  # Initialize final dataframe outside the year loop
  
  # Loop through years
  for (y in unique(grid_sub$year)) {
    grid_year <- grid_sub[grid_sub$year == y, ]
    print(y)
    Ori_tracks_year<-Ori_tracks[Ori_tracks$year == y, ] #v0.7
    
    # Initialize yearly dataframe
    m_pred_df <- data.frame()
    
    # Loop through months
    for (m in unique(grid_year$month)) {
      grid_month <- grid_year[grid_year$month == m, ]
      print(m)
      
      Ori_tracks_month<-Ori_tracks_year[Ori_tracks_year$month == m, ] #v0.7
      
      # Prediction for the month
      pred <- data.frame(predict(gam, se.fit = TRUE, newdata = grid_month, type = "response")) #original
      # pred <- data.frame(predict(gam, se.fit = TRUE, newdata = grid_month, type = "response", na.action=na.omit)) #updated v1.0 possibly
      
      # Extract lon, lat, and predicted values
      pred_map <- data.frame(
        lon = grid_month$lon,
        lat = grid_month$lat,
        pred = pred$fit,
        pred_se=pred$se.fit, #v0.7 add se
        date = as.Date(grid_month$date)
      )
      pred_map$year <- year(pred_map$date)
      pred_map$month <- month(pred_map$date)
      
      
      # Append monthly data to yearly dataframe
      m_pred_df <- rbind(m_pred_df, pred_map)
      
      pred_map2=pred_map #v1.7
      pred_map2$Monat<-format(pred_map2$date, format = "%B") #addition v1.5 post first rendition -->v1.7
      Monat2<-unique(pred_map2$Monat)
      
      # Generate and save monthly map
      gmap <- ggplot() +
        geom_tile(data = pred_map2, aes(x = lon, y = lat, fill = pred)) + ##changed this to pred_map2 v1.7
        geom_point(data = Ori_tracks_month, aes(x = lon, y = lat), size = 0.75) + #switched to ori_track_month v0.7
        geom_sf(data = Az) +
        # scale_fill_fermenter(n.breaks = 10, palette = "RdYlBu", limits = c(0, 1), na.value = NA) +
        scale_fill_distiller(name = "Prediction",n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
        guides(fill = guide_coloursteps(show.limits = TRUE, even.steps = FALSE)) +
        # ggtitle(paste(#sp,
        #   y, m, sep = "-")) +
        ggtitle(paste(Monat2, y, sep = " "))+
        xlab("Longitude") + ylab("Latitude") +
        theme_classic()+
        theme(
          plot.title = element_text(size=12, face = "bold", hjust = 0.5), #v1.7
          strip.text = element_text(size = 10, face = "bold"), #new post addition v1.5
          legend.text = element_text(size = 8),
          legend.key.height = unit(1.5, "cm"),
          legend.key.width = unit(0.5, "cm")
        )
      
      ggsave(gmap, file = paste0(mapsdir, "/pred_gam_", y, "_", m, ".tiff"), width = 6.88, height = 4.25, bg = "white", compression = "lzw", dpi = 300)
    }
    
    # Append yearly data to final dataframe
    final_pred_df <- rbind(final_pred_df, m_pred_df)
  }
  # Save the final dataframe
  saveRDS(final_pred_df, paste(mapsdir, "/final_pred_df_AllYrs_Wide.RDS", sep = ""))
  final_pred_df<-readRDS(file = paste(mapsdir, "/final_pred_df_AllYrs_Wide.RDS", sep = ""))
  
  ##Create final map
  month_levels <- month.name[2:10]  # February to October
  
  # final_pred_df$Monat<-format(final_pred_df$date, format = "%B")
  final_pred_df$Monat <- factor(
    format(final_pred_df$date, "%B"),levels = month_levels)
  
  Ori_tracks$Monat <- factor(format(Ori_tracks$date, "%B"),
                             levels = month_levels)
  
  
  # scale_fill_distiller(n.breaks = 10, palette = "Spectral",limits = c(0, 1), na.value = NA) + ##v0.6
  
  predmap <- ggplot() +
    geom_tile(data = final_pred_df, aes(x = lon, y = lat, fill = pred)) +
    # geom_point(data = Ori_tracks, aes(x = lon, y = lat, 
                                      # color="Animal tracks" ##added post creation v1.5 to test
                                      # ), size = 0.5) + #size was 0.76
    geom_sf(data = Az) +
    scale_fill_distiller(
      name = "Prediction", 
      n.breaks=10, 
      palette = "Spectral",
      limits = c(0, 1),
      na.value = NA
    ) +
    # scale_color_manual(
    #   name = NULL,
    #   values = c("Animal tracks" = "black")
    # ) + #rename and color to black ##Temp turn off v1.6 Feb2026
    
    # guides(
    #   fill  = guide_coloursteps(order = 1),
    #   color = guide_legend(
    #     order = 2,
    #     override.aes = list(size = 3)
    #   ) ##make bigger dots v1.5 post
    # ) +  ##v1.7 -->i've turned this off, but need to check this fxn because it was affecting the maps and dropped the 0 & 1 values from the legend.
    scale_x_continuous(
      breaks = seq(-32, -23, by = 2), #changed 3 to 2
      labels = function(x) abs(x)
    ) +
    scale_y_continuous(
      breaks = seq(35, 40, by = 1), #changed 1 to 2
      labels = function(y) abs(y)
    )+
    facet_grid(Monat ~ year) +
    xlab("Longitude (°W)") + ylab("Latitude (°N)") +
    theme_classic()+
    theme(
      axis.title = element_text(size=13, face="bold"), 
      axis.text = element_text(size=11),
      strip.text = element_text(size = 12),
      legend.title = element_text(size=12, face="bold"),
      legend.text = element_text(size = 11),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  predmap
  
  ggsave(filename =paste0(mapsdir,"/final_Predmap_AllYears.tiff"), predmap, 
         width = 9, height = 12,#width = 7.5, height = 10,
         dpi = 600, bg="white", compression = "lzw") #width 14.5
  
  ##Standard error map
  sdmap <- ggplot() +
    geom_tile(data = final_pred_df, aes(x = lon, y = lat, fill = pred_se)) +
    # geom_point(data = Ori_tracks, aes(x = lon, y = lat), size = 0.5) + #size was 0.76
    geom_sf(data = Az) +
    scale_fill_distiller(
      name = "SE Prediction", 
      n.breaks=10, 
      palette = "Spectral",
      limits = c(0, 1),
      na.value = NA
    ) +
    # scale_color_manual(
    #   name = NULL,
    #   values = c("Animal tracks" = "black")
    # ) + #rename and color to black ##Temp turn off v1.6 Feb2026
    scale_x_continuous(
      breaks = seq(-32, -23, by = 2), #changed 3 to 2
      labels = function(x) abs(x)
    ) +
    scale_y_continuous(
      breaks = seq(35, 40, by = 1), #changed 1 to 2
      labels = function(y) abs(y)
    )+
    facet_grid(Monat ~ year) +
    xlab("Longitude (°W)") + ylab("Latitude (°N)") +
    theme_classic()+
    theme(
      axis.title = element_text(size=13, face="bold"), 
      axis.text = element_text(size=11),
      strip.text = element_text(size = 12),
      legend.title = element_text(size=12, face="bold"),
      legend.text = element_text(size = 11),
      legend.key.height = unit(1.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    )
  
  sdmap
  
  ggsave(filename =paste0(mapsdir,"/final_SE_Predmap_AllYears.tiff"), sdmap, 
         width = 9, height = 12,#width = 7.5, height = 10,
         dpi = 600, bg="white", compression = "lzw") #width 14.5
  
  
  return(list(predplot=predmap,sdplot=sdmap, final_DF = final_pred_df))
} #End of Fxn_AllMaps_Pred for making the new maps for 9 months across 3 years for Sperm whales
##Functions to plot enviro importance ----

Fxn_VarImpPlot <- function(ItemName,
                           ModelName,
                           FilePath, 
                           DataType # "Normal" or "SpatialCalc", 
                           
){
  require(ggplot2)
  require(dplyr)
  require(openxlsx)
  
  # Load variable name references
  cleanVar_df <- read.xlsx(paste0(SDM_HD, "/EnviroVars_StockList.xlsx"), sheet = "AllVars")
  cleanVar_df_spat <- read.xlsx(paste0(SDM_HD, "/EnviroVars_StockList.xlsx"), sheet = "ForMeansSD")
  
  # Load importance data
  load(paste0(FilePath, "/", ModelName,"_Output.RData")) #to pull out type v1.1
  var_imp <- read.csv(file = paste0(FilePath, "/Evaluations/EnviroVariables_ImportanceDF.csv")) %>%
    rename(R_VarName = expl.var)
  var_imp$R_VarName <- as.factor(var_imp$R_VarName)
  
  # Initialize plot variable
  plot_result <- NULL
  
  if (DataType == "Normal") {
    
    var_imp_df <- left_join(var_imp, cleanVar_df, by = "R_VarName") %>%
      mutate(norm.imp = as.numeric(norm.imp)) %>%
      filter(!is.na(Enviro_Var))
    
    plot_result <- ggplot(data = var_imp_df, 
                          aes(x = reorder(Enviro_Var, norm.imp), y = norm.imp)) +
      geom_bar(stat = 'identity', fill = "steelblue") +
      coord_flip() +
      labs(x = '', y = "Percent of Importance", #v1.2 added x axis, but remember it is inverted
           title = "Variable importance") +
      theme_classic() +
      theme(plot.title = element_text(size = 9, hjust = 0.5),
            axis.text = element_text(size = 7),
            axis.title = element_text(size=8), #v1.2
            # plot.subtitle = element_text(size = 7.5, hjust = 0)
            plot.caption = element_text(size = 7.5, hjust = 0) #v1.1
      )
    
    print(plot_result)
    
    ggsave(plot_result, filename = paste0(FilePath, "/Evaluations/EnviroImp_barplot.tiff"),
           width = 5.5, height = 3.2, bg = "white", compression = "lzw", dpi = 300)
    
  } else if (DataType == "SpatialCalc") {
    
    # var_imp$R_VarName_SpatCalc <- paste0(var_imp$R_VarName, "_Mean")
    var_imp<-var_imp%>%rename(R_VarName_SpatCalc =  R_VarName ) #v1.1
    var_imp_df_spat <- left_join(var_imp, cleanVar_df_spat, by = "R_VarName_SpatCalc") %>%
      mutate(norm.imp = as.numeric(norm.imp)) %>%
      filter(!is.na(Enviro_Var))
    
    plot_result <- ggplot(data = var_imp_df_spat, 
                          aes(x = reorder(Enviro_Var, norm.imp), y = norm.imp)) +
      geom_bar(stat = 'identity', fill = "steelblue") +
      coord_flip() +
      labs(x = '', y = "Percent of Importance", #v1.2 added x axis
           title = "Variable importance", 
           # subtitle = paste0("Data for spatial calculation: ", type)
           caption  = paste0("Data for spatial calculation: ", type) #v1.1
           ) +
      theme_classic() +
      theme(plot.title = element_text(size = 9, hjust = 0.5),
            axis.text = element_text(size = 7),
            axis.title= element_text(size=8), #v1.2
            # plot.subtitle = element_text(size = 7.5, hjust = 0)
            plot.caption = element_text(size = 7.5, hjust = 0) #v1.1
            )
    
    print(plot_result)
    ggsave(plot_result, filename = paste0(FilePath, "/Evaluations/EnviroImp_barplot.tiff"),
           width = 5.5, height = 3.2, bg = "white", compression = "lzw", dpi = 300)
    
  } else {
    stop("Invalid DataType. Please use 'Normal' or 'SpatialCalc'.")
  }
  
  return(plot_result)
} #end Fxn_VarImpPlot v1.0


##Supplementary Functions----
##Function to clean pred df based on gams v1.5
Fxn_clean_pred_for_gam <- function(pred_df, gam) {
  
  gam_vars <- setdiff(names(gam$model), c("presence", "(Intercept)"))
  log_vars <- gam_vars[grepl("^log_", gam_vars)]
  
  pred_df <- pred_df %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(log_vars),
        ~ ifelse(is.infinite(.), NA, .)
      )
    ) %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(gam_vars),
        ~ !is.na(.)
      )
    )
  
  return(pred_df)
} ##Function to clean up vars and match the gams v1.5
