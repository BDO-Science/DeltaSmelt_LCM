# Identify your working directory for saving outputs of interest
root <- "~/GitHub/DeltaSmelt_LCM"
setwd(root)

zoop_root <- file.path(root,"Zooplankton_input_conversion_for_LCM")
zoop_original_root <- file.path(root,"Zooplankton_input_original")

library(reshape)
library(tidyverse)
library(lubridate)
library(sf)
library(deltamapr)
library(readxl)


IMBR_LCM_comparison <- ggplot()+
  geom_sf(data=WW_Delta,fill="cadetblue1", color="cadetblue1")+
  geom_sf(data=R_EDSM_Regions_1617P1,fill=NA,color="black")+
  geom_sf(data=R_DSIBM,fill=NA, color="red")+
  geom_sf_label(data=R_DSIBM, aes(label = SUBREGION),size=2.5,color="red")+
  geom_sf_label(data=R_EDSM_Regions_1617P1, aes(label = Region),size=2.5,color="black")+
  theme_bw()
  
#Print out the map
tiff(filename=file.path(zoop_root,"Figure_Map_compare_IBMR_LCM.tiff"), units="in",type="cairo", bg="white", height=10, 
     width=12, res=300, compression="lzw")
IMBR_LCM_comparison
dev.off()
  
# SW Suisun will be used to represent "Far West" region in the Delta Smelt LCM
# Confluence and all other downstream regions will be used to represent "West" region in the Delta Smelt LCM
  
  
######
# Check the stations used in the Delta Smelt LCM, check just 1994 and on
zoop_data_original <- read.csv(file.path(zoop_original_root,"ZooMysid_74_19_df.csv")) %>% 
  mutate(Date=as.Date(Date)) %>% dplyr::filter(year(Date)>1993) 

zoop_data_stations <- unique(zoop_data_original[c("Station")])

# Read station list table from Arthur Barros
zoop_station_list <- read_xlsx(file.path(zoop_root,"Station_Lookup.xlsx")) %>% rename(Station=StationNZ)

# Join lat and long data
zoop_data_stations <- zoop_data_stations %>% left_join(zoop_station_list %>% select(Station, Core, lat_decimal, lon_decimal)) %>%
  st_as_sf(coords=c("lon_decimal", "lat_decimal"), crs=4326, remove=FALSE)

# Plot
zoop_station_map <- ggplot()+
  geom_sf(data=WW_Delta,fill="cadetblue1", color="cadetblue1")+
  geom_sf(data=R_EDSM_Regions_1617P1,fill=NA,color="black")+
  geom_sf(data=R_DSIBM,fill=NA, color="red")+
  geom_sf(data=zoop_data_stations, color="purple")+
  theme_bw()
zoop_station_map

# Print out the map
tiff(filename=file.path(zoop_root,"Figure_Map_zoop_stations.tiff"), units="in",type="cairo", bg="white", height=10, 
     width=12, res=300, compression="lzw")
zoop_station_map
dev.off()

# Compare the region overlap and count number of stations
zoop_join<- st_join(st_as_sf(zoop_data_stations),st_transform(R_DSIBM,crs=4326))
zoop_join<- st_join(zoop_join,st_transform(R_EDSM_Regions_1617P1,crs=4326))

# Far West = SW Suisun
# West = 2 stations in Suisun Marsh, 2 stations in the Confluence, 2 in SE Suisun, 1 in NW Suisun
