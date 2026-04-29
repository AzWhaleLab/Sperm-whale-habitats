<img width="1276" height="622" alt="flukup_aerial2 copy" src="https://github.com/user-attachments/assets/68c565f2-08ff-42a4-a7dc-bb58a2bec235" />

# Sperm-whale-habitats
This repo is meant to help in the development of species distribution models (SDM) of sperm whales using tracking data.
Code supporting the publication “Ecological drivers of movements and habitats of a deep diver: 
the elusive sperm whale”. 

The goal of this repo is to use tracking data to develop Species Distribution Models (SDMs) by following a similar workflow as Hazen et al 2017 and Pérez-Jorge et al 2020.  The workflow (see below), includes creating pseudo-absences, extract environmental data for both tracking data and prediction maps, and applying Generalized Additive Mixed Models (GAMMs). 

**Required data:**
Tracking data with lat, long, time. 
Location data used for modeling are accessible in the data folder: these are only the state-space-modeled tracks (the Correlated Random Walk (CRW) at a 6-hour interval). 
The authors of this data are Mareike D. Duffing Romero, Rui Prieto, Sergi Pérez-Jorge, Mónica Silva

For standardizing tracking data by given time interval and method, we used the ‘aniMotum’ R package using codes from https://ianjonsen.github.io/aniMotum/index.html.

**Code authors:** Mareike D. Duffing Romero (OKEANOS/UAç), Elliott Hazen (UCSC), Sergi Pérez-Jorge (OKEANOS/UAç). 

**Functions:** 

<ins>Environmental Data extraction:</ins> function to extract environmental data from netcdf files including netcdf files that are in 3D or 4D or for extracting variables from multiple netcdf files that are differentiated by having a single date per file. Other functions include create spatial means and SD within a 3x3 grid cell. 

<ins>Model Evaluation:</ins> functions that allow to evaluate GAMM output models, create prediction maps, assess variable contribution and create gam plots, call Function_GAM_EvalMaps_version. R

**Scripts that call functions for SDM workflow:** 
1) Create pseudo-absences using CRW approach. Calls 1_Create PsuedoAbsences_CRW.R
2) Extract environmental data from netcdf files and include in the data frame with tracks and pseudo-absences. Calls 2a_Download_environmental data.R. For extracting environmental variables and then creating a prediction dataframe to develop a dataframe with all variables, call 2b_Download_predictiondata.R
3) Data exploration and assessment of collinearity of your variables, call 3_DataExploration.R
4) Development of SDMs via GAMMS for main models used and pertinent evaluation codes used for this project. 

**Relevant manuscripts:**
Hazen, E. L., Palacios, D. M., Forney, K. A., Howell, E. A., Becker, E., Hoover, A. L., Irvine, L., DeAngelis, M., Bograd, S. J., Mate, B. R., & Bailey, H. (2017). WhaleWatch: A dynamic management tool for predicting blue whale density in the California Current. Journal of Applied Ecology, 54(5), 1415–1428. https://doi.org/10.1111/1365-2664.12820

Hazen, E. L., Abrahms, B., Brodie, S., Carroll, G., Welch, H., & Bograd, S. J. (2021). Where did they not go? Considerations for generating pseudo-absences for telemetry-based habitat models. Movement Ecology, 9(1), 5. https://doi.org/10.1186/s40462-021-00240-2

Pérez‐Jorge, S., Tobeña, M., Prieto, R., Vandeperre, F., Calmettes, B., Lehodey, P., & Silva, M. A. (2020). Environmental drivers of large‐scale movements of baleen whales in the mid‐North Atlantic Ocean. Diversity and Distributions, 26(6), 683–698. https://doi.org/10.1111/ddi.13038
