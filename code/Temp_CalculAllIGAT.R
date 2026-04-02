### AllIGATR_temps v0.0.1
### Addison Rice (github.com/addison-rice)
### 2026-01-05
###
###
### Calculate AllIGATR temperatures for GDGT data.
###
### This script runs using tidyverse. Make sure it's installed! 
### 
###
### This code takes a csv as input, and produces a csv as output.
### You can use the template provided to ensure that column names will work
### with the code, or you can use the commented-out block on lines 66-72
### to rename your columns.
### 
### Make sure your input has a column with a unique sample ID!!!!
### 
### Define your input and output file names and file paths at the start, and
### the code should run without further issues. It does take some time and will
### use quite a bit of memory, so start with a few GDGT data at first to make
### sure you won't crash your computer. 
###
### There will be two output files. One will have the mean AllIGATR model 
### temperatures (both SST and 0-200m temperatures) for each sample ID. The 
### other contains the temperature given for each of the 150 models. If you 
### don't want the full results taking up space on your computer, comment out 
### line 223.
### 
### 
### Please cite Rice et al. (2026) when using this code in any publication.


### make sure file names end in ".csv"
input_filepath <- "../data/"
input_filename <- "myGDGTdata.csv"

output_filepath <- "../results/"
output_summary_filename <- "myGDGTdata_AllIGATR_summary.csv"
output_full_filename <- "myGDGTdata_AllIGATR_fullresults.csv"


### make sure that these files exist on your computer!
all_model_fits <- read.csv("../Models/noRS_low2over3_model_fits.csv")
hightempmodels <- read.csv("../Models/hightempmodels_noRScoretop_seds_2over3under10_sst.csv")
hightemptoothmodels <- read.csv("../Models/hightempmodels_noRScoretop_seds_2over3under10_deep.csv")





### If you used the csv template, there is no need to adjust code below this line!








library(tidyverse)


### load data
data <- read.csv(str_c(input_filepath, input_filename)) %>%
  ### If needed, edit the myXXname variables to match your input file and uncomment the following lines
  # rename(gdgt0 = mygdgt0name,
  #        gdgt1 = mygdgt1name,
  #        gdgt2 = mygdgt2name,
  #        gdgt3 = mygdgt3name,
  #        cren = mycrenname,
  #        creniso = mycrenisoname,
  #        sampleID = mysampleIDname) %>%
  select(c("sampleID", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))













### define model GDGT distributions for a 
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
           model_creniso = 100*model_creniso/model_sumGDGTs)
  return(model_data)
}






expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))


### match model temperatures to data and summarize alligatr temperatures and the data-model distance
temp_match <- function(input_data){
  
  ### recalculate GDGT fractions to ensure no errors
  calc_data <- select(input_data, c("sampleID", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso")) %>%
    mutate(gdgt0 = as.numeric(gdgt0),
           gdgt1 = as.numeric(gdgt1),
           gdgt2 = as.numeric(gdgt2),
           gdgt3 = as.numeric(gdgt3),
           cren = as.numeric(cren),
           creniso = as.numeric(creniso),
           across(everything(), ~replace(.x, is.nan(.x), 0)), ### assume nan or na values are 0
           across(everything(), ~replace(.x, is.na(.x), 0)),
           sumGDGTs = gdgt0+gdgt1+gdgt2+gdgt3+cren+creniso,
           gdgt0 = 100*gdgt0/sumGDGTs,
           gdgt1 = 100*gdgt1/sumGDGTs,
           gdgt2 = 100*gdgt2/sumGDGTs,
           gdgt3 = 100*gdgt3/sumGDGTs,
           cren = 100*cren/sumGDGTs,
           creniso = 100*creniso/sumGDGTs)
  
  
  
  long_data_all <- data.frame()
  output_data_all <- data.frame()
  
  
  for(modeltype in unique(all_model_data$model_type)){
    model_data <- subset(all_model_data, model_type == modeltype)
    
    data_model <- expand.grid.df(model_data, calc_data) %>%
      mutate(temp_dist = ((gdgt0_weights*(gdgt0 - model_gdgt0))^2 + 
                            (gdgt1_weights*(gdgt1 - model_gdgt1))^2 +
                            (gdgt2_weights*(gdgt2 - model_gdgt2))^2 +
                            (gdgt3_weights*(gdgt3 - model_gdgt3))^2 +
                            (cren_weights*(cren - model_cren))^2 +
                            (creniso_weights*(creniso - model_creniso))^2)^(1/2),
             flatweight_dist = ((0.2*(gdgt0 - model_gdgt0))^2 + 
                                 ((gdgt1 - model_gdgt1))^2 +
                                 ((gdgt2 - model_gdgt2))^2 +
                                 ((gdgt3 - model_gdgt3))^2 +
                                 (0.2*(cren - model_cren))^2 +
                                 ((creniso - model_creniso))^2)^(1/2))
    
    
    data_match <- data_model %>%
      mutate(row = str_c(sampleID, "-", model_id)) %>%
      group_by(row) %>%
      summarise(dist = min(flatweight_dist, na.rm = TRUE),
                temp = mean(reported_Temp[which.min(temp_dist)], na.rm = TRUE)) %>%
      ungroup() %>%
      separate_wider_delim(row, "-", names = c("sampleID", "model_id")) %>%
      mutate(model_depth = case_when(str_detect(model_id, c("sst")) ~ "sst",
                                     .default = "deep"),
             model_num = str_extract(model_id, "(?<=_)[^_]*$"),
             outlier_model = case_when((model_depth == "sst" & model_num %in% hightempmodels$x) | 
                                       (model_depth == "deep" & model_num %in% hightemptoothmodels$x) ~ 1,
                                       TRUE ~ 0)) %>%
      ungroup()
    
    long_data <- left_join(input_data, data_match, by = "sampleID")
    
    summary_data <- subset(long_data, temp > -5 & temp < 60) %>% ### remove model endpoint matches
      group_by(across(all_of(c("model_depth", "sampleID", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso")))) %>%
      summarize(countmodels = n(), 
                alligatr_Temp_allmodels = mean(temp, na.rm = TRUE),
                stdevTemp_allmodels = sd(temp, na.rm = TRUE),
                alligatr_Dist_allmodels = mean(dist, na.rm = TRUE),
                stdevDist_allmodels = sd(dist, na.rm = TRUE),
                alligatr_Temp_nooutliermodels = mean(temp[which(outlier_model == 0)], na.rm = TRUE),
                stdevTemp_nooutliermodels = sd(temp[which(outlier_model == 0)], na.rm = TRUE),
                alligatr_Dist_nooutliermodels = mean(dist[which(outlier_model == 0)], na.rm = TRUE),
                stdevDist_nooutliermodels = sd(dist[which(outlier_model == 0)], na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(model_type = modeltype)
    
    output_data <- left_join(input_data, summary_data, by = c("sampleID", "gdgt0", "gdgt1", "gdgt2", "gdgt3", "cren", "creniso"))
    
    long_data_all <- rbind(long_data_all, long_data)
    output_data_all <- rbind(output_data_all, output_data)
    
  }
     
  
  
  
  return(list(long_data_all, output_data_all))
  
}


all_model_data <- model_GDGTs(all_model_fits)

results <- temp_match(data)

write.csv(results[[2]], str_c(output_filepath, output_summary_filename))
write.csv(results[[1]], str_c(output_filepath, output_full_filename))


