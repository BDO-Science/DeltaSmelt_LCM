#Create a separate base data to allow inclusion of 2015-2016 data not present in the CSAMP dataset
root <- "~/GitHub/DeltaSmelt_LCM"
setwd(root)

hydro_root <- file.path(root,"Hydro model output")
salinity_root <- file.path(root,"X2_salinity_analysis")
salinity_original_root <- file.path(root,"X2_salinity_analysis_original")

library(reshape)
library(tidyverse)
library(stringr)
library(lubridate)
library(conflicted)
library(deltamapr)
require(discretewq)
conflict_prefer("rename", "dplyr")

#Use same datasets as CSAMP Delta Smelt SDM
wq_data <-  wq(Sources = c("EMP", "STN", "FMWT", "EDSM", "DJFMP",
                           "SKT", "20mm", "Suisun", 
                           "Baystudy", "USBR", "USGS_SFBS")) %>%
  filter(!is.na(Latitude)&!is.na(Longitude)) %>%
  st_as_sf(coords=c("Longitude", "Latitude"),crs=st_crs("WGS84"))

#Delta Smelt IBM layer
IBM_layer <- R_DSIBM

#Join IBM layer with wq data
wq_data_IBM <- st_transform(wq_data, crs = st_crs(IBM_layer))
wq_data_IBM<- st_join(wq_data_IBM,IBM_layer)
st_geometry(wq_data_IBM) <- NULL # remove geometry, coerce to data.frame

#Summarize salinity by region, month, and year
wq_data_sum <- wq_data_IBM %>% group_by(Year,Month,SUBREGION) %>% summarise(Salinity=mean(Salinity,na.rm = T)) %>% rename(year=Year,month=Month,region=SUBREGION) %>%
  filter(year %in% c(1995:2016)&!is.na(region)) %>% rename(sal_base=Salinity) %>% ungroup() %>%
  #Add filled in data that is somehow present in the CSAMP version for 1997 January NE Suisun
  add_row(year = 1997, month=1, region = "NE Suisun", sal_base=0.217891455)

#Export data
write.csv(wq_data_sum,file.path(salinity_root,"base_salinity_data.csv"),row.names=T)
