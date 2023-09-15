# Identify your working directory for saving outputs of interest
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
conflict_prefer("rename", "dplyr")

# Load X2 data from CalSim3
EXP1_data <- read.csv(file.path(hydro_root,"EXP1_CalSim3_data.csv")) 
EXP3_data <- read.csv(file.path(hydro_root,"EXP3_CalSim3_data.csv")) 
NAA_data <- read.csv(file.path(hydro_root,"NAA_CalSim3_data.csv")) 
Alt1_data <- read.csv(file.path(hydro_root,"Alt1_CalSim3_data.csv")) 
Alt2v1woTUCP_data <- read.csv(file.path(hydro_root,"Alt2v1woTUCP_CalSim3_data.csv")) 
Alt2v1wTUCP_data <- read.csv(file.path(hydro_root,"Alt2v1wTUCP_CalSim3_data.csv")) 
Alt2v2noTUCP_data <- read.csv(file.path(hydro_root,"Alt2v2noTUCP_CalSim3_data.csv")) 
Alt2v3noTUCP_data <- read.csv(file.path(hydro_root,"Alt2v3noTUCP_CalSim3_data.csv")) 
Alt4_data <- read.csv(file.path(hydro_root,"Alt4_CalSim3_data.csv")) 

# Combine X2 data
x2_data <- NAA_data %>% select(Date,X2_current) %>% mutate(Scenario="NAA") %>%
  bind_rows((EXP1_data %>% select(Date,X2_current) %>% mutate(Scenario="EXP1"))) %>%
  bind_rows((EXP3_data %>% select(Date,X2_current) %>% mutate(Scenario="EXP3"))) %>%
  bind_rows((Alt1_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt1"))) %>%
  bind_rows((Alt2v1woTUCP_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt2v1woTUCP"))) %>%
  bind_rows((Alt2v1wTUCP_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt2v1wTUCP"))) %>%
  bind_rows((Alt2v2noTUCP_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt2v2noTUCP"))) %>%
  bind_rows((Alt2v3noTUCP_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt2v3noTUCP"))) %>%
  bind_rows((Alt4_data %>% select(Date,X2_current) %>% mutate(Scenario="Alt4")))

x2_data <- na.omit(x2_data) %>% rename(X2 = X2_current) %>% mutate(Month=month(Date))

# Create X2 data frame for all relevant regions
x2_data_expanded <- crossing(x2_data, Region=c("NW Suisun","SW Suisun","NE Suisun","SE Suisun","Confluence", "Suisun Marsh"))

# Load Brian Crawford's final salinity-X2 model
salX2mod <- readRDS(file.path(salinity_original_root,"Model_sal_X2.Rdata")) 

###### Rescale X2
# Ensure format is the same with Brian Crawford's original script, so those lines are copy and pasted below
#Import and create X2 dataframe
X2.IBMR <- read.csv(file.path(salinity_original_root,"X2_baseline.csv"))
X2.dat <- melt(X2.IBMR,id.var=c("Year"))
colnames(X2.dat)[2:3] <- c("Month","X2")
levels(X2.dat$Month) <- 1:12
X2.dat$Month <- as.integer(X2.dat$Month)

X2.sal.dat <- expand.grid(Year = 1995:2014, Month = 1:12, Region = c("Confluence","East Delta","Lower Sacramento","Lower San Joaquin", "NE Suisun","NW Suisun",        
                                                                     "Sacramento River","SE Suisun","South Delta","Suisun Marsh","SW Suisun","Yolo Bypass")) # Create dataframe
regions <- levels(X2.sal.dat$Region)
X2.sal.dat$Region <- factor(X2.sal.dat$Region, levels=regions[c(12,7,9,2,3,4,1,8,5,10,11,6)]) # Reorder region names
X2.sal.dat <- left_join(x=X2.sal.dat, y=X2.dat[,c(1:ncol(X2.dat))], by=c("Year","Month")) # Join X2 data

sal_base <- read.csv(file.path(salinity_original_root,"sal_baselineIBMR.csv"), header=F) # From IBMR data. dimensions: 12 subregions (rows), 240 year-months from Jan 1995 to Dec 2014 (columns)

years <- 20
months <- 12
n.strata <- 12
X2.sal.dat$sal <- NA
dim(sal_base)
for (s in 1:n.strata){
  for (y in 1:years){
    for (m in 1:months){ 
      X2.sal.dat$sal[which(X2.sal.dat$Year==y+1994 & X2.sal.dat$Month==m & X2.sal.dat$Region==levels(X2.sal.dat$Region)[s])] <- sal_base[s,(12*(y-1)+m)]
    }
  }
}

X2.sal.dat$Month <- X2.sal.dat$Month - 6.5
range(X2.sal.dat$Month)
sum(is.na(X2.sal.dat$sal))

# Convert to z-score prior to model predict
x2_data_expanded <- x2_data_expanded %>% mutate(X2_original=X2, X2=(X2-mean(X2.sal.dat$X2))/sd(X2.sal.dat$X2))


# Use the CSAMP X2-Salinity model to convert CalSim3 X2 values to salinity
x2_data_expanded$salinity<-predict(salX2mod,x2_data_expanded)

# Ensure that there will be no negative salinity values and use the minimum value in Sam's conversion table
x2_data_expanded$salinity <- ifelse(x2_data_expanded$salinity<0.1,0.1,x2_data_expanded$salinity)

# Finalize data format
x2_data_expanded <- x2_data_expanded %>%
  mutate(year=year(Date)) %>% rename(region=Region, month=Month) %>% select(-Date,-X2,-X2_original) %>%
  spread(Scenario,salinity) 

# Rename column names to sal_
colnames(x2_data_expanded)[4:ncol(x2_data_expanded)] <- paste("sal", colnames(x2_data_expanded)[4:ncol(x2_data_expanded)] , sep = "_")


#Export output file for  model input
write.csv(x2_data_expanded,file.path(salinity_root,"converted_salinity_data.csv"),row.names=F)

