### AllIGATR_create v0.3.0
### Addison Rice (github.com/addison-rice)
### 2024-05-17
###
###
### Algorithmically optimizes curves to fit GDGT distributions in modern and ancient data and temperatures in modern data


library(tidyverse)
library(rmoo)
library(doParallel)
library(rdist)



all_gdgt_data <- read_csv("../data/model_fit/rattanasriampaipong2022_compilation.csv", 
                          col_types = cols(Latitude = col_double(), 
                                           Longitude = col_double(), modernWaterDepth = col_double(), 
                                           sampleDepth = col_double(), `Sample Age (Ma)` = col_double(), 
                                           reported_Temp = col_double(), fGDGT_0 = col_double(), 
                                           fGDGT_1 = col_double(), fGDGT_2 = col_double(), 
                                           fGDGT_3 = col_double(), fGDGT_cren = col_double(), 
                                           `fGDGT_cren'` = col_double(), BITindex = col_double())) %>%
  rename(type = dataType_level2) %>%
  mutate(MI = (`fGDGT_1` + `fGDGT_2` + `fGDGT_3`) / (`fGDGT_1` + `fGDGT_2` + `fGDGT_3` + `fGDGT_cren` + `fGDGT_cren'`) > 0.3,
         total_isos = `fGDGT_0` + `fGDGT_1` + `fGDGT_2` + `fGDGT_3` + `fGDGT_cren` + `fGDGT_cren'`,
         gdgt0 = 100*`fGDGT_0`/total_isos,
         gdgt1 = 100*`fGDGT_1`/total_isos,
         gdgt2 = 100*`fGDGT_2`/total_isos,
         gdgt3 = 100*`fGDGT_3`/total_isos,
         cren = 100*`fGDGT_cren`/total_isos,
         creniso = 100*`fGDGT_cren'`/total_isos) 


all_gdgt_selected_cols <- select(all_gdgt_data, c("type", "reported_Temp", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso", "MI", "BITindex", "Source"))


### load Devika's dataset, just isos
oh_gdgt_iso_data <- read.csv("../data/model_fit/DAS_ohgdgt_surfacesediment_data_for_alligatr.csv") %>%
  rename(reported_Temp = Annual.mean.SST...C.a,
         Temp_0_200 = Annual.mean.temperature.0.200.m...C.a, # Annual mean temperature 0-200 m (°C)a
         BITindex = BITd) %>%
  mutate(total_isos = GDGT.0 + GDGT.1 + GDGT.2 + GDGT.3 + cren + cren.,
         gdgt0 = 100*GDGT.0/total_isos,
         gdgt1 = 100*GDGT.1/total_isos,
         gdgt2 = 100*GDGT.2/total_isos,
         gdgt3 = 100*GDGT.3/total_isos,
         cren = 100*cren/total_isos,
         creniso = 100*cren./total_isos,
         MI = (gdgt1 + gdgt2 + gdgt3) / (gdgt1 + gdgt2 + gdgt3 + cren + creniso) > 0.3,
         type = "Core top") %>%
  select(c("type", "reported_Temp", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso", "MI", "BITindex", "Region", "Temp_0_200"))











### optimization function
gdgt_optim <- function(reference_dirs, sed_gdgt = sedimentary_data_all, temp_gdgt = coretops_subset, temp_depth = "SST"){
  temp_range <- -50:600 # temperatures up to 60C
  temp_range <- temp_range/10
  
  x0 = reference_dirs
  input_model <- data.frame(mu = numeric(length = length(x0)/5),
                            sig = numeric(length = length(x0)/5),
                            max = numeric(length = length(x0)/5),
                            min = numeric(length = length(x0)/5),
                            weights = numeric(length = length(x0)/5)) 
  
  input_model <- input_model %>%
    mutate(mu = x0[1:(length(x0)/5)],
           sig = x0[(1+length(x0)/5):(2*length(x0)/5)],
           max = x0[(1+2*length(x0)/5):(3*length(x0)/5)],
           min = x0[(1+3*length(x0)/5):(4*length(x0)/5)],
           weights = x0[(1+4*length(x0)/5):(length(x0))])
  
  model_data <- as.data.frame(temp_range) %>%
    rename(reported_Temp = temp_range) %>%
    mutate(gdgt0 = input_model$min[1] + (input_model$max[1] - input_model$min[1])*exp(-(temp_range-input_model$mu[1])^2/(2*input_model$sig[1]^2)),
           gdgt1 = input_model$min[2] + (input_model$max[2] - input_model$min[2])*exp(-(temp_range-input_model$mu[2])^2/(2*input_model$sig[2]^2)),
           gdgt2 = input_model$min[3] + (input_model$max[3] - input_model$min[3])*exp(-(temp_range-input_model$mu[3])^2/(2*input_model$sig[3]^2)),
           gdgt3 = input_model$min[4] + (input_model$max[4] - input_model$min[4])*exp(-(temp_range-input_model$mu[4])^2/(2*input_model$sig[4]^2)),
           cren = input_model$min[5] + (input_model$max[5] - input_model$min[5])*exp(-(temp_range-input_model$mu[5])^2/(2*input_model$sig[5]^2)),
           creniso = input_model$min[6] + (input_model$max[6] - input_model$min[6])*exp(-(temp_range-input_model$mu[6])^2/(2*input_model$sig[6]^2)),
           sumGDGTs = gdgt0+gdgt1+gdgt2+gdgt3+cren+creniso,
           gdgt0 = 100*gdgt0/sumGDGTs, # normalize
           gdgt1 = 100*gdgt1/sumGDGTs,
           gdgt2 = 100*gdgt2/sumGDGTs,
           gdgt3 = 100*gdgt3/sumGDGTs,
           cren = 100*cren/sumGDGTs,
           creniso = 100*creniso/sumGDGTs)
  
  
  
  
  
  ### optimize on distance for sedimentary data
  dist_optim <- function(model_test){
    
    model_data_dist <- model_data %>%
      mutate(gdgt0 = gdgt0*0.2,  # lower weighting for distance calculation
             cren = cren*0.2) %>%
      select(c("gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))
    
    sed_gdgt <- sed_gdgt %>%
      mutate(gdgt0 = gdgt0*0.2,  # lower weighting for distance calculation
             cren = cren*0.2) %>%
      select(c("gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))
    
    data_model_distance <- rdist::cdist(sed_gdgt, model_data_dist)
    sed_gdgt$distance <- apply(data_model_distance, 1, FUN = min)
    
    rmse_dist <- sqrt(mean((sed_gdgt$distance)^2))
    
    return(rmse_dist)
  }
  
  
  
  
  
  ### optimize on temperature for curated core tops
  temp_optim <- function(input_model){
    temp_gdgt$out_temps <- apply(temp_gdgt,1,function(x){
      model_data$reported_Temp[which.min((as.numeric(x["gdgt0"])*input_model$weights[1] - as.numeric(model_data$gdgt0)*
                                            input_model$weights[1])^2 +
                                           (as.numeric(x["gdgt1"])*input_model$weights[2] - as.numeric(model_data$gdgt1)*
                                              input_model$weights[2])^2 +
                                           (as.numeric(x["gdgt2"])*input_model$weights[3] - as.numeric(model_data$gdgt2)*
                                              input_model$weights[3])^2 +
                                           (as.numeric(x["gdgt3"])*input_model$weights[4] - as.numeric(model_data$gdgt3)*
                                              input_model$weights[4])^2 +
                                           (as.numeric(x["cren"])*input_model$weights[5] - as.numeric(model_data$cren)*
                                              input_model$weights[5])^2 +
                                           (as.numeric(x["creniso"])*input_model$weights[6] - as.numeric(model_data$creniso)*
                                              input_model$weights[6])^2)]})
    
    if(temp_depth == "SST"){
      rmse_temp <- sqrt(mean((temp_gdgt$out_temps - temp_gdgt$reported_Temp)^2))
    } else{
      rmse_temp <- sqrt(mean((temp_gdgt$out_temps - temp_gdgt[[temp_depth]])^2))
    }
    
    return(rmse_temp)
  }
  
  
  
  
  ### reduce distributions too far from the closed sum effect
  closed_sum <- sqrt(mean((model_data$sumGDGTs - 100)^2))
  
  ### avoid over-parameterizing to closed sum
  if(closed_sum < 5){
    closed_sum = (closed_sum+5)/2
    }
  
  
  
  ### call the functions
  return(c(dist_optim(input_model), temp_optim(input_model), closed_sum))
}





### model restrictions (for SSTs). Each line is GDGT-0,-1,-2,-3,cren,cren'
model_mins_sst <- c(-15,10,25,32,35,30, # minimum temperature optimum
                    5,5,5,2,5,5, # minimum sigma (peak width)
                    35,5,5,1,60,5, # minimum value for peak maximum
                    0,0,0,0,0,0, # minimum value for curve minimum
                    0,0,0,0,0,0) # minimum weighting
model_maxs_sst <- c(10,30,45,50,100,100, # maximum temperature optimum
                    50,20,20,20,50,50, # maximum sigma (peak width)
                    100,12,12,10,100,30, # maximum value for peak maximum
                    70,2,2,2,70,2, # maximum value for curve minimum
                    1,1,1,1,1,1) # maximum weighting
# model_visfit <- c(0, 22, 30, 36, 50, 50,
#                   18, 11, 11, 8, 16, 16,
#                   55, 7, 7, 3, 80, 16,
#                   1, 0.2, 0.2, 0.25, 45, 0.1,
#                   0.2, 1, 1, 1, 0.2, 1)
model_seeds_sst <- as.matrix(read.csv("../Models/seeds/sst_seeds.csv"))



### model restrictions (for 0-200m temps). Each line is GDGT-0,-1,-2,-3,cren,cren'
model_mins_deep <- c(-15,5,15,25,30,30, # minimum temperature optimum
                     2,2,2,2,2,2, # minimum sigma (peak width)
                     35,5,5,1,60,5, # minimum value for peak maximum
                     0,0,0,0,0,0, # minimum value for curve minimum
                     0,0,0,0,0,0) # minimum weighting
model_maxs_deep <- c(10,30,40,55,80,80, # maximum temperature optimum
                     50,20,20,20,50,50, # maximum sigma (peak width)
                     100,12,12,10,100,30, # maximum value for peak maximum
                     70,2,2,2,70,2, # maximum value for curve minimum
                     1,1,1,1,1,1) # maximum weighting
# model_visfit_deep <- c(0, 18, 23, 28, 40, 40,
#                        13, 8, 8, 7, 12, 12,
#                        55, 7, 7, 3, 80, 16,
#                        1, 0.2, 0.2, 0.25, 45, 0.1,
#                        0.2, 1, 1, 1, 0.2, 1)
model_seeds_deep <- as.matrix(read.csv("../Models/seeds/deep_seeds.csv"))

gdgts = c("gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso")
params = c("mu", "sig", "max", "min", "weights")

all_params <- expand.grid(gdgts = gdgts, params = params) %>%
  mutate(colnames = paste(gdgts, params, sep = "_"))


save_model <- function(gdgt_model, model_name){
  model_pareto_fits <- gdgt_model@fitness %>%
    matrix(ncol = 3)
  model_pareto_fits <- as.data.frame(model_pareto_fits) %>%
    mutate(V1 = round(V1, 2),
           V2 = round(V2, 2),
           V3 = round(V3, 3)) %>%
    rename(Distance = V1,
           Temperature = V2,
           ClosedSum = V3)
  model_pareto_fits$ParetoFrontRank <- gdgt_model@front
  
  
  all_models <- gdgt_model@population %>%
    matrix(ncol = 30) %>%
    as.data.frame()
  colnames(all_models) <- all_params$colnames
  
  all_models <- cbind(all_models, model_pareto_fits)
  
  write.csv(all_models, str_c("../Models/model_parameters/",model_name,".csv"))
}


### wrapper function
alligatr_run_model <- function(coretop_data, sediment_data, modelname){
  gdgt_model_sst <- rmoo(type = "real-valued",
                         fitness = gdgt_optim,
                         algorithm = "NSGA-II",
                         lower = model_mins_sst,
                         upper = model_maxs_sst,
                         sed_gdgt = sediment_data,
                         temp_gdgt = coretop_data,
                         nvars = 30,
                         popSize = 2000,
                         nObj = 3,
                         maxiter = 50,
                         suggestions = model_seeds_sst,
                         monitor = FALSE,
                         parallel = TRUE)
  
  save_model(gdgt_model_sst, str_c(modelname, "_sst"))
  
  gdgt_model_deep <- rmoo(type = "real-valued",
                         fitness = gdgt_optim,
                         algorithm = "NSGA-II",
                         lower = model_mins_deep,
                         upper = model_maxs_deep,
                         temp_depth = "Temp_0_200",
                         sed_gdgt = sediment_data,
                         temp_gdgt = subset(coretop_data, Temp_0_200 > -10),
                         nvars = 30,
                         popSize = 2000,
                         nObj = 3,
                         maxiter = 50,
                         suggestions = model_seeds_deep,
                         monitor = FALSE,
                         parallel = TRUE)
  
  save_model(gdgt_model_deep, str_c(modelname, "_deep"))
}


### curated datasets to use for the model
sedimentary_data_all <- subset(all_gdgt_selected_cols, gdgt0*gdgt1*gdgt2*gdgt3*cren*creniso >0 &  
                                 type %in% c("Core top (>100m)", "Core top (0-100m)", "Mesozoic", "Early Cenozoic", "Late Cenozoic") & 
                                 MI == FALSE & gdgt0 < 60 & creniso < 40 & BITindex < 0.5 & !is.na(gdgt0))

### curated core top subset. reported_Temp > -10 removes NA values
coretops_subset <- subset(oh_gdgt_iso_data, MI == FALSE & gdgt0 < 60 & reported_Temp > -10 & BITindex < 0.5)







### no Red Sea core tops, all sediments
alligatr_run_model(coretop_data = subset(coretops_subset, !Region == "Red Sea"), 
                   sediment_data = sedimentary_data_all, 
                   modelname = "noRScoretop_allseds")

### no Red Sea, no cold core tops, all sediments
alligatr_run_model(coretop_data = subset(coretops_subset, !Region == "Red Sea" & reported_Temp > 5), 
                   sediment_data = sedimentary_data_all, 
                   modelname = "noRSnocoldcoretop_allseds")

### no Red Sea core tops, no Mesozoic sediments
alligatr_run_model(coretop_data = subset(coretops_subset, !Region == "Red Sea"), 
                   sediment_data = subset(sedimentary_data_all, !type == "Mesozoic"), 
                   modelname = "noRScoretop_noMesoseds")




# calculate GDGT2/3 values
sedimentary_data_all$gdgt2over3 <- sedimentary_data_all$gdgt2/sedimentary_data_all$gdgt3
coretops_subset$gdgt2over3 <- coretops_subset$gdgt2/coretops_subset$gdgt3




### core top and seds with GDGT-2/GDGT-3 < 10, no Red Sea
alligatr_run_model(coretop_data = subset(coretops_subset, gdgt2over3 < 10 & !Region == "Red Sea"), 
                   sediment_data = subset(sedimentary_data_all, gdgt2over3 < 10), 
                   modelname = "noRScoretop_seds_2over3under10")

### core top and seds with GDGT-2/GDGT-3 < 10, no Red Sea, no Mesozoic
alligatr_run_model(coretop_data = subset(coretops_subset, gdgt2over3 < 10 & !Region == "Red Sea"), 
                   sediment_data = subset(sedimentary_data_all, gdgt2over3 < 10 & !type == "Mesozoic"), 
                   modelname = "noRScoretop_noMesoseds_2over3under10")

### core top and seds with GDGT-2/GDGT-3 < 10, no Red Sea, no cold core tops
alligatr_run_model(coretop_data = subset(coretops_subset, gdgt2over3 < 10 & !Region == "Red Sea" & reported_Temp > 5), 
                   sediment_data = subset(sedimentary_data_all, gdgt2over3 < 10), 
                   modelname = "noRSnocoldcoretop_seds_2over3under10")


