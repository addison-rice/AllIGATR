### AllIGATR_temps v0.1.0
### Addison Rice (github.com/addison-rice)
### 2024-05-24
###
###
### define high temperature outlier models for noRS_low2over3 models


library(tidyverse)
library(rdist)


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
  
}





### match GDGT data to model temperatures output from the model_GDGTs function
temp_match <- function(input_data, all_model_data, model_name){
  model_data <- subset(all_model_data, model_id == model_name)
  model_temps <- apply(input_data,1,function(x){
    model_data$reported_Temp[mean(which.min((as.numeric(x["gdgt0"])*model_data$gdgt0_weights -
                                               as.numeric(model_data$model_gdgt0)*model_data$gdgt0_weights)^2 +
                                              (as.numeric(x["gdgt1"])*model_data$gdgt1_weights -
                                                 as.numeric(model_data$model_gdgt1)*model_data$gdgt1_weights)^2 +
                                              (as.numeric(x["gdgt2"])*model_data$gdgt2_weights -
                                                 as.numeric(model_data$model_gdgt2)*model_data$gdgt2_weights)^2 +
                                              (as.numeric(x["gdgt3"])*model_data$gdgt3_weights -
                                                 as.numeric(model_data$model_gdgt3)*model_data$gdgt3_weights)^2 +
                                              (as.numeric(x["cren"])*model_data$cren_weights -
                                                 as.numeric(model_data$model_cren)*model_data$cren_weights)^2 +
                                              (as.numeric(x["creniso"])*model_data$creniso_weights -
                                                 as.numeric(model_data$model_creniso)*model_data$creniso_weights)^2))
    ]})
  model_data_dist <- model_data %>%
    mutate(gdgt0 = model_gdgt0*0.2,  # lower weighting for distance calculation
           cren = model_cren*0.2) %>%
    rename(gdgt1 = model_gdgt1,
           gdgt2 = model_gdgt2,
           gdgt3 = model_gdgt3,
           creniso = model_creniso) %>%
    select(c("gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))
  
  sed_gdgt <- input_data %>%
    mutate(gdgt0 = gdgt0*0.2,  # lower weighting for distance calculation
           cren = cren*0.2) %>%
    select(c("gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))
  
  data_model_distance <- rdist::cdist(sed_gdgt, model_data_dist)
  distance <- apply(data_model_distance, 1, FUN = min)
  
  out_data <- input_data
  out_data[[model_id]] <- unlist(model_temps)
  out_data[[str_c(model_id, "_dist")]] <- unlist(distance) 
  
  return(out_data)
}



### first load dataset and calculate temperatures for sedimentary data
all_gdgt_selected_cols <- read_csv("../data/model_fit/rattanasriampaipong2022_compilation.csv", 
                          col_types = cols(Latitude = col_double(), 
                                           Longitude = col_double(), modernWaterDepth = col_double(), 
                                           sampleDepth = col_double(), `Sample Age (Ma)` = col_double(), 
                                           reported_Temp = col_double(), fGDGT_0 = col_double(), 
                                           fGDGT_1 = col_double(), fGDGT_2 = col_double(), 
                                           fGDGT_3 = col_double(), fGDGT_cren = col_double(), 
                                           `fGDGT_cren'` = col_double(), BITindex = col_double())) %>%
  rename(type = dataType_level2,
         sampleID = sampleName) %>%
  mutate(MI = (`fGDGT_1` + `fGDGT_2` + `fGDGT_3`) / (`fGDGT_1` + `fGDGT_2` + `fGDGT_3` + `fGDGT_cren` + `fGDGT_cren'`) > 0.3,
         total_isos = `fGDGT_0` + `fGDGT_1` + `fGDGT_2` + `fGDGT_3` + `fGDGT_cren` + `fGDGT_cren'`,
         gdgt0 = 100*`fGDGT_0`/total_isos,
         gdgt1 = 100*`fGDGT_1`/total_isos,
         gdgt2 = 100*`fGDGT_2`/total_isos,
         gdgt3 = 100*`fGDGT_3`/total_isos,
         cren = 100*`fGDGT_cren`/total_isos,
         creniso = 100*`fGDGT_cren'`/total_isos) %>%
  select(c("type", "reported_Temp", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso", "MI", "BITindex", "Source", "sampleID"))

all_gdgt_selected_cols$type[which(all_gdgt_selected_cols$type %in% c("Core top (>100m)", "Core top (0-100m)"))] <- "Core top"
all_gdgt_selected_cols$type[which(all_gdgt_selected_cols$type %in% c("Culture - AOA"))] <- "Culture"
all_gdgt_selected_cols$type[which(all_gdgt_selected_cols$type %in% c("Early Cenozoic", "Late Cenozoic"))] <- "Cenozoic"



### load model fits
all_model_fits <- read.csv("../Models/noRS_low2over3_model_fits.csv")
### calculate model GDGTs
all_model_gdgts <- model_GDGTs(all_model_fits)



### calculate SSTs and 0-200m temperatures for the full dataset
All_gdgt_alligatr <- all_gdgt_selected_cols
for(model_id in unique(subset(all_model_gdgts, model_type == "noRScoretop_seds_2over3under10_sst")$model_id)){
  All_gdgt_alligatr <- temp_match(All_gdgt_alligatr, subset(all_model_gdgts, model_type == "noRScoretop_seds_2over3under10_sst"), model_id)
}
All_gdgt_alligatr$ID <- seq.int(nrow(All_gdgt_alligatr))
All_gdgt_alligatrLong <- All_gdgt_alligatr %>%
  pivot_longer(cols = starts_with("noRScoretop_seds_2over3under10_sst_"), names_to = "param", values_to = "value")
All_gdgt_alligatrLong$temp_dist <- "temp"
All_gdgt_alligatrLong$temp_dist[endsWith(All_gdgt_alligatrLong$param, "_dist")] <- "dist"
All_gdgt_alligatrLong$model_num <- str_extract(All_gdgt_alligatrLong$param, "(?<=sst_)(\\d+)")
All_gdgt_alligatrLong <- select(All_gdgt_alligatrLong, -c("param")) %>%
  pivot_wider(names_from = temp_dist, values_from = value)

all_gdgt_tooth <- all_gdgt_selected_cols
for(model_id in unique(subset(all_model_gdgts, model_type == "noRScoretop_seds_2over3under10_deep")$model_id)){
  all_gdgt_tooth <- temp_match(all_gdgt_tooth, subset(all_model_gdgts, model_type == "noRScoretop_seds_2over3under10_deep"), model_id)
}
all_gdgt_tooth$ID <- seq.int(nrow(all_gdgt_tooth))
all_gdgt_toothLong <- all_gdgt_tooth %>%
  pivot_longer(cols = starts_with("noRScoretop_seds_2over3under10_deep_"), names_to = "param", values_to = "value")
all_gdgt_toothLong$temp_dist <- "temp"
all_gdgt_toothLong$temp_dist[endsWith(all_gdgt_toothLong$param, "_dist")] <- "dist"
all_gdgt_toothLong$model_num <- str_extract(all_gdgt_toothLong$param, "(?<=deep_)(\\d+)")
all_gdgt_toothLong <- select(all_gdgt_toothLong, -c("param")) %>%
  pivot_wider(names_from = temp_dist, values_from = value)



### Assess high temperature outlier models using histograms
gdgt_fit_hist_seds <- subset(All_gdgt_alligatrLong, type %in% c("Cenozoic", "Core top", "Mesozoic") & MI == FALSE  & gdgt0 < 60 & creniso < 40 & BITindex < 0.5 & !is.na(gdgt0)) %>%
  ggplot(aes(x = temp,
             fill = type)) +
  geom_histogram() +
  facet_wrap(~type, ncol = 1, scales = "free_y")
ggsave("../results/gdgt_fit_hist_seds_sst.pdf", useDingbats = FALSE, width = 4, height = 6)

gdgt_fit_hist_seds_tooth <- subset(all_gdgt_toothLong, type %in% c("Cenozoic", "Core top", "Mesozoic") & MI == FALSE  & gdgt0 < 60 & creniso < 40 & BITindex < 0.5 & !is.na(gdgt0)) %>%
  ggplot(aes(x = temp,
             fill = type)) +
  geom_histogram() +
  facet_wrap(~type, ncol = 1, scales = "free_y")
ggsave("../results/gdgt_fit_hist_seds_tooth.pdf", useDingbats = FALSE, width = 4, height = 6)



### eliminate models that often yield high temperatures (>50C for SST, >35C for 0-200m; exclude 60C due to model limit edge effects)
gdgt_models_hightempcount <- subset(All_gdgt_alligatrLong, type %in% c("Cenozoic", "Core top", "Mesozoic") & MI == FALSE  & gdgt0 < 60 & creniso < 40 & BITindex < 0.5 & !is.na(gdgt0) & temp > 50) %>%
  group_by(model_num) %>%
  summarize(countover50 = n())
hightempmodels <- gdgt_models_hightempcount$model_num[which(gdgt_models_hightempcount$countover50 > 100)]
write.csv(hightempmodels, "../Models/hightempmodels_noRScoretop_seds_2over3under10_sst.csv")

gdgt_tooth_models_hightempcount <- subset(all_gdgt_toothLong, type %in% c("Cenozoic", "Core top", "Mesozoic") & MI == FALSE  & gdgt0 < 60 & creniso < 40 & BITindex < 0.5 & !is.na(gdgt0) & temp > 40) %>%
  group_by(model_num) %>%
  summarize(countover40 = n())
hightemptoothmodels <- gdgt_tooth_models_hightempcount$model_num[which(gdgt_tooth_models_hightempcount$countover40 > 100)]
write.csv(hightemptoothmodels, "../Models/hightempmodels_noRScoretop_seds_2over3under10_deep.csv")

