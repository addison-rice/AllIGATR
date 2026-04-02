# AllIGATR

AllIGATR stands for ALL Isoprenoid GDGT (glycerol dibiphytanyl glycerol tetraether) Adaptation to Temperature in R, and is a calibration of marine GDGT distributions to temperature for use in paleoclimate applications. This code (specifically, Temp_CalculAllIGAT.R) can be used to calculate proxy SST (sea surface temperature) and TO0TH (depth-integrated Temperatures of 0 to Two Hundred meters). The code used to create temperature adaptation curves based on core top and ancient GDGT data is also provided. For full details on model development and performance, check out the preprinted article (CITATION).

How to use this code
Downloading the full repository and keeping the folder structure is simplest because the code uses relative file locations.

Packages 
All of the .R files use tidyverse. For calculating proxy temperature, this is the only package needed.
For creating and validating AllIGATR curves, the code also uses rdist, doParallel, and rmoo (version 0.3.0).

Calculating proxy temperatures
Proxy temperatures are calculated using Temp_CalculAllIGAT.R. This code takes an input .csv file, for which the template myGDGTdata.csv is provided in the data folder. Make sure that each row has a unique sampleID, otherwise you will get average values of any identically named rows. It outputs two .csv files, which will appear in the results folder if no changes are made to the code. One output file will contain the full results of the 150 best-performing models, and the other will provide a summary with the mean and standard deviation of proxy temperatures.

To get temperatures for your GDGT data, either edit the myGDGTdata.csv file or edit the file names and file paths in the Temp_CalculAllIGAT.R file. Ensure that your file uses these column headers: 
sampleID	gdgt0	gdgt1	gdgt2	gdgt3	cren	creniso

Then, save and run the Temp_CalculAllIGAT.R file. The code can take some time and lots of computer memory for large datasets, so consider running a portion of your data at a time to avoid crashing the computer.

Creating AlllIGATR curves
To create your own AllIGATR curves, start with AllIGATR_create_v0.3.0.R. This code uses the NSGAII algorithm in rmoo to optimize model fit to GDGT data on multiple criteria: temperature fit to surface sediments, distance fit to surface and ancient sediments, and a closed sum close to 100%. The temperature response of each GDGT is modeled as a Gaussian curve, and the minimum, maximum, mu, sigma, and GDGT weighting are determined algorithmically. To reduce computational time, these values are constrained, and the algorithm is seeded with a visually fit model and with previous algorithm output (seed files are in Models/seeds/). This file outputs the files in Models/model_parameters, which each contain thousands of AllIGATR curve fits. Note that the original output files contained many models which were identical to five decimal points. These were manually deleted.

To identify the best-performing models, use AllIGATR_model_curves.R. This code ranks the models by performance in temperature and distance fit, and limits the GDGT-3 maximum value to better reflect ancient sediments. The model fits for the top 150 models from each fitting scenario are saved as Models/all_model_fits_top150.csv. The fits for the recommended fitting scenario are saved as noRS_low2over3_model_fits.csv. The GDGT distributions from the top 150 models of each fitting scenario are calculated at 0.1C intervals between -5 and 60C and saved in Models/model_gdgts/.

The AllIGATR models from the recommended fitting scenario are then checked for high temperature outliers using AllIGATR_high_temp_outliers.R. This file calculates the model temperatures of surface and ancient sediments and plots histograms of temperature prediction. Cutoff temperatures were chosen based on these, and any model which predicts temperatures above these thresholds for 100 or more sediments were flagged as high temperature outlier models.
