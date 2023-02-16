# Identify your working directory for saving outputs of interest
root <- "~/GitHub/DeltaSmelt_LCM_test"
setwd(root)

hydro_root <- file.path(root,"Hydro model output")
salinity_root <- file.path(root,"X2_salinity_analysis")
salinity_original_root <- file.path(root,"X2_salinity_analysis_original")


library(tidyverse)
library(stringr)
library(lubridate)

# Load X2 data from CalSim3
NAA_data <- read.csv(file.path(hydro_root,"NAA_CalSim3_data.csv")) 
D1641_data <- read.csv(file.path(hydro_root,"D1641_CalSim3_data.csv")) %>%
  mutate(X2_current=X2_current+1)

# Combine X2 data
x2_data <- NAA_data %>% select(Date,X2_current) %>% mutate(Scenario="NAA") %>%
  bind_rows((D1641_data %>% select(Date,X2_current) %>% mutate(Scenario="D1641"))) 

x2_data <- na.omit(x2_data) %>% rename(X2=X2_current) %>% mutate(Month=month(Date))

# Create X2 data frame for all relevant regions
x2_data_expanded <- crossing(x2_data, Region=c("NW Suisun","SW Suisun","NE Suisun","SE Suisun","Confluence", "Suisun Marsh"))

# Load Brian Crawford's final salinity-X2 model
salX2mod <- readRDS(file.path(salinity_original_root,"Model_sal_X2.Rdata")) 

# Use the CSAMP X2-Salinity model to convert CalSim3 X2 values to salinity
x2_data_expanded$salinity<-predict(salX2mod,x2_data_expanded)


#Export output file for  model input
write.csv(x2_data_expanded,file.path(salinity_root,"converted_salinity_data.csv"),row.names=F)

