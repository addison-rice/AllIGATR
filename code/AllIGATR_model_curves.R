### AllIGATR_model_curves v0.1.0
### Addison Rice (github.com/addison-rice)
### 2024-05-24
###
###
### Get model GDGT curves from the parameters output by AllIGATR_create


library(tidyverse)


### take the models from AllIGATR_create as input and calculate the associated GDGT distributions at 0.1C intervals
model_GDGTs <- function(model_list, temp_range = -50:600, model_name = NULL){
  reported_Temp <- temp_range/10
  
  model_data <- merge(model_list, as.data.frame(reported_Temp)) %>%
    mutate(model_gdgt0 = gdgt0_min + (gdgt0_max[1] - gdgt0_min)*exp(-(reported_Temp - gdgt0_mu)^2/(2*gdgt0_sig^2)),
           model_gdgt1 = gdgt1_min + (gdgt1_max[1] - gdgt1_min)*exp(-(reported_Temp - gdgt1_mu)^2/(2*gdgt1_sig^2)),
           model_gdgt2 = gdgt2_min + (gdgt2_max[1] - gdgt2_min)*exp(-(reported_Temp - gdgt2_mu)^2/(2*gdgt2_sig^2)),
           model_gdgt3 = gdgt3_min + (gdgt3_max[1] - gdgt3_min)*exp(-(reported_Temp - gdgt3_mu)^2/(2*gdgt3_sig^2)),
           model_cren = cren_min + (cren_max[1] - cren_min)*exp(-(reported_Temp - cren_mu)^2/(2*cren_sig^2)),
           model_creniso = creniso_min + (creniso_max[1] - creniso_min)*exp(-(reported_Temp - creniso_mu)^2/(2*creniso_sig^2)),
           model_sumGDGTs = model_gdgt0+model_gdgt1+model_gdgt2+model_gdgt3+model_cren+model_creniso,
           model_gdgt0 = 100*model_gdgt0/model_sumGDGTs,
           model_gdgt1 = 100*model_gdgt1/model_sumGDGTs,
           model_gdgt2 = 100*model_gdgt2/model_sumGDGTs,
           model_gdgt3 = 100*model_gdgt3/model_sumGDGTs,
           model_cren = 100*model_cren/model_sumGDGTs,
           model_creniso = 100*model_creniso/model_sumGDGTs,
           model_tex86 = (model_gdgt2+model_gdgt3+model_creniso)/(model_gdgt1+model_gdgt2+model_gdgt3+model_creniso))
  write.csv(model_data, str_c("../Models/model_gdgts/", model_name, ".csv"))
  
}




### open the files generated from AllIGATR_create
models <- c("noRScoretop_allseds", "noRScoretop_noMesoseds", "noRScoretop_noMesoseds_2over3under10",  
            "noRScoretop_seds_2over3under10", "noRSnocoldcoretop_allseds", "noRSnocoldcoretop_seds_2over3under10")

model_fits <- data.frame()
all_model_fits <- model_fits

for(model in 1:length(models)){
  for(depth in c("_sst", "_deep")){
    model_fits <- unique(read.csv(str_c("../Models/model_parameters/", models[model], depth, ".csv"))[2:35]) %>%
      arrange(Temperature, Distance, ParetoFrontRank, ClosedSum) %>%  ### priotitize temperature fit first
      subset(gdgt3_max > 3 & gdgt3_max < 4) %>%                       ### restrict the GDGT-3 maximum to better fit ancient data
      slice_head(n = 150)                                             ### take the top 150 models
    model_fits$model_id <- seq.int(nrow(model_fits))                  ### assign a number
    model_fits <- model_fits %>%
      mutate(model_type = str_c(models[model], depth),
             model_id = str_c(model_type, "_", model_id))
    
    all_model_fits <- rbind(all_model_fits, model_fits)               ### compile the best models from each fitting scenario
    model_GDGTs(model_fits, model_name = str_c(models[model], depth)) ### create a csv with the GDGT distributions for -5 to 60C for the top 150 models for each fitting scenario
    
  }
  
}

write.csv(all_model_fits, "../Models/all_model_fits_top150.csv")
write.csv(subset(all_model_fits, model_type %in% c("noRScoretop_seds_2over3under10_sst", "noRScoretop_seds_2over3under10_deep")), 
          "../Models/noRS_low2over3_model_fits.csv")


model_param_summary <- select(all_model_fits, !model_id) %>%
  group_by(model_type) %>%
  summarize(across(c(1:34), list(mean = mean, min = min, max = max)))
write.csv(model_param_summary, "../results/model_param_summary.csv")

model_summary <- all_model_fits %>%
  group_by(model_type) %>%
  summarise(minTempRMSE = min(Temperature),
            maxTempRMSE = max(Temperature),
            meanTempRMSE = mean(Temperature),
            minDistRMSE = min(Distance),
            maxDistRMSE = max(Distance),
            meanDistRMSE = mean(Distance))
write.csv(model_summary, "../results/model_performance_summary.csv")
