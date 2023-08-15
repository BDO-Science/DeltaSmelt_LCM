# Identify your working directory for saving outputs of interest
root <- "~/GitHub/DeltaSmelt_LCM"
setwd(root)

output_root <- file.path(root,"Hydro model output")

# The following libraries need to be installed and loaded
# NOTE: You also need to have HEC-DSSVue installed on your computer
# See: https://www.hec.usace.army.mil/software/hec-dssvue/downloads.aspx

library(tidyverse)
library(stringr)
library(lubridate)
library(rJava)


#############
#Read DSS file

# The following function for is used for turning CalSim time stamps into R dates. 

from_time_stamp <- function(x) {
  day_ref <- as.Date("1899-12-30")
  return(day_ref+x/1440)
}



# Run this workaround if your R session crashes when running .jinit() - below
# This issue occurs in base R versions 4.2 and later
# In lieu of this workaround, you can also install a patched R version 
# E.g., https://cran.r-project.org/bin/windows/base/rpatched.html

replacement <- function(category = "LC_ALL") {
  
  if (identical(category, "LC_MESSAGES"))
    return("")
  
  category <- match(category, .LC.categories)
  if (is.na(category)) 
    stop("invalid 'category' argument")
  .Internal(Sys.getlocale(category))
  
}
base <- asNamespace("base")
environment(replacement) <- base
unlockBinding("Sys.getlocale", base)
assign("Sys.getlocale", replacement, envir = base)
lockBinding("Sys.getlocale", base)


# This code establishes your java connection with HEC-DSSVue

# Specify your own location for 'HEC-DSSVue'
dss_location <- "C:\\Program Files\\HEC\\HEC-DSSVue\\" 

# Specify your own location for the 'jar' sub-folder
# This identifies all possible java executables that meet be needed
jars <- c(list.files("C:\\Program Files\\HEC\\HEC-DSSVue\\jar")) 

jars <- paste0(dss_location, "jar/", jars)

# Specify your own location for the 'lib' sub-folder
libs <- "-Djava.library.path=C:\\Program Files\\HEC\\HEC-DSSVue\\lib\\"

.jinit(classpath = jars, parameters = libs)


##########
# Function to assemble the dataset

# Identify the DSS file you want to access with dss_input

dss_data_pull_LCME<-function(dss_input="D:\\2023-01-06 - CalSim3 example file for ReROC\\CalSim3_2040MED_120722_DRAFT_wDWRv705update_wCCdraftBC\\DSS\\output\\CS3_L2020_DV_2021_ext_2040MED"){
  # Open the DSS file through rJava
  dssFile <- .jcall("hec/heclib/dss/HecDss", "Lhec/heclib/dss/HecDss;",   method="open", dss_input)
  #Delta Outflow (OUT in Dayflow)
  java.NDOI_MIN <- dssFile$get("/CALSIM/NDOI_MIN/FLOW//1MON/L2020A/") 
  java.NDOI_ADD <- dssFile$get("/CALSIM/NDOI_ADD/FLOW//1MON/L2020A/") 
  OUTFLOW=data.frame(Date=java.NDOI_MIN$times %>% from_time_stamp,OUTFLOW=java.NDOI_MIN$values+java.NDOI_ADD$values)
  #Old and Middle River flow (OMR)
  java.OMR <- dssFile$get("/CALSIM/C_OMR014/CHANNEL//1MON/L2020A/") 
  OMR=data.frame(Date=java.OMR$times %>% from_time_stamp,OMR=java.OMR$values)
  #X2 (previous month)
  java.X2 <- dssFile$get("/CALSIM/X2_PRV/X2-POSITION-PREV//1MON/L2020A/") 
  X2=data.frame(Date=java.X2$times %>% from_time_stamp,X2_prev=java.X2$values)
  
  final_data_frame= OUTFLOW %>% left_join(OMR) %>% left_join(X2) %>%
    mutate(X2_current=lead(X2_prev,n=1))
  return(final_data_frame)
}


#Use the function to create data frame
EXP1_data <- dss_data_pull_LCME(dss_input="C:\\Users\\bmahardja\\Documents\\2023-08-14 - DSS files for 2021 ROC\\_Reclamation_LTO2021\\BA\\EXP1\\CalSim3\\2022MED\\CalSim3_EXP1_2022MED_rev09d_071323_dynGWSW\\DSS\\output\\LTO_BA_EXP1_2022MED_DV")
EXP3_data <- dss_data_pull_LCME(dss_input="C:\\Users\\bmahardja\\Documents\\2023-08-14 - DSS files for 2021 ROC\\_Reclamation_LTO2021\\BA\\EXP3\\CalSim3\\CalSim3_EXP3_2022MED_rev09e_072523_dynGWSW\\DSS\\output\\LTO_BA_EXP3_2022MED_DV")
NAA_data <- dss_data_pull_LCME(dss_input="C:\\Users\\bmahardja\\Documents\\2023-08-14 - DSS files for 2021 ROC\\_Reclamation_LTO2021\\BA\\NAA\\CalSim3\\2022MED\\CalSim3_NAA_2022MED_07192023\\DSS\\output\\CS3_NAA_2022MED_071923_L2020ADV")
Alt2_wTUCP_data <- dss_data_pull_LCME(dss_input="C:\\Users\\bmahardja\\Documents\\2023-08-14 - DSS files for 2021 ROC\\LTO_Alt2PA_07252023_wTUCP\\DSS\\output\\LTO_Alt2PA_wTUCP_0724_tucp1100_NMBin3_DV")
Alt2_woTUCP_data <- dss_data_pull_LCME(dss_input="C:\\Users\\bmahardja\\Documents\\2023-08-14 - DSS files for 2021 ROC\\LTO_Alt2PA_07252023_woTUCP\\DSS\\output\\LTO_Alt2PA_woTUCP_0725_DV")

#Export DSS output files for Delta Smelt LCM model input
write.csv(EXP1_data,file.path(output_root,"EXP1_CalSim3_data.csv"),row.names=F)
write.csv(EXP3_data,file.path(output_root,"EXP3_CalSim3_data.csv"),row.names=F)
write.csv(NAA_data,file.path(output_root,"NAA_CalSim3_data.csv"),row.names=F)
write.csv(Alt2_wTUCP_data,file.path(output_root,"Alt2_wTUCP_CalSim3_data.csv"),row.names=F)
write.csv(Alt2_woTUCP_data,file.path(output_root,"Alt2_woTUCP_CalSim3_data.csv"),row.names=F)

##### Calculate flow input for the Delta Smelt LCME

EXP1_data <- EXP1_data %>% mutate(scenario="EXP1")
EXP3_data <- EXP3_data %>% mutate(scenario="EXP3")
NAA_data <- NAA_data %>% mutate(scenario="NAA")
Alt2_wTUCP_data <- Alt2_wTUCP_data %>% mutate(scenario="Alt2_wTUCP")
Alt2_woTUCP_data <- Alt2_woTUCP_data %>% mutate(scenario="Alt2_woTUCP")

combined_data <- bind_rows(EXP1_data,EXP3_data,NAA_data,Alt2_wTUCP_data,Alt2_woTUCP_data)

#Summer Delta Outflow covariate
data_DeltaOutflow <- combined_data %>% mutate(Calendar_Year=year(Date),Month=month(Date)) %>% 
  #Filter June-August per Smith et al. 2021
  filter(Month %in% c(6:8)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Calendar_Year %in% c(1995:2016)) %>%
  #sum flow by month
  mutate(SumOutflow_Jun_Aug=case_when(Month==6 ~ OUTFLOW*30,
                                            Month==7 ~OUTFLOW*31,
                                            Month==8 ~OUTFLOW*31)) %>%
  #summarize by year
  group_by(Calendar_Year,scenario) %>% summarise(DeltaOutflow_Jun_Aug=sum(SumOutflow_Jun_Aug))

#Apr-May OMR covariate
data_OMR_Apr_May <- combined_data %>% mutate(Calendar_Year=year(Date),Month=month(Date)) %>% 
  #Filter Apr-May per Smith et al. 2021
  filter(Month %in% c(4:5)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Calendar_Year %in% c(1995:2016)) %>%
  #summarize by year
  group_by(Calendar_Year,scenario) %>% summarise(OMR_Apr_May=mean(OMR))

#June OMR covariate
data_OMR_June <- combined_data %>% mutate(Calendar_Year=year(Date),Month=month(Date)) %>% 
  #Filter June per Smith et al. 2021
  filter(Month %in% c(6)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Calendar_Year %in% c(1995:2016)) %>%
  rename(OMR_June=OMR) %>% select(Calendar_Year,scenario,OMR_June)

#Dec-Jan OMR covariate
data_OMR_Dec_Jan <- combined_data %>% 
  #Use water year instead for this one
  mutate(Water_Year=ifelse(month(Date)>9,year(Date)+1,year(Date)),Month=month(Date)) %>% 
  #Filter Dec-Jan per Smith et al. 2021
  filter(Month %in% c(12,1)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Water_Year %in% c(1995:2016)) %>%
  group_by(Water_Year,scenario) %>% summarise(OMR_Dec_Jan=mean(OMR))

#February OMR covariate
data_OMR_Feb <- combined_data %>% mutate(Calendar_Year=year(Date),Month=month(Date)) %>% 
  #Filter February per Smith et al. 2021
  filter(Month %in% c(2)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Calendar_Year %in% c(1995:2016)) %>%
  rename(OMR_Feb=OMR) %>% select(Calendar_Year,scenario,OMR_Feb)

#March OMR covariate
data_OMR_Mar <- combined_data %>% mutate(Calendar_Year=year(Date),Month=month(Date)) %>% 
  #Filter March per Smith et al. 2021
  filter(Month %in% c(3)) %>% 
  #Filter year 1995-2016 for Delta Smelt LCM
  filter(Calendar_Year %in% c(1995:2016)) %>%
  rename(OMR_Mar=OMR) %>%  select(Calendar_Year,scenario,OMR_Mar)

##Merge datasets
data_flowinput_LCM <- data_DeltaOutflow %>% left_join(data_OMR_Apr_May) %>% left_join(data_OMR_June) %>% left_join(data_OMR_Feb) %>%
  left_join(data_OMR_Mar)

#Export data
write.csv(data_flowinput_LCM,file.path(output_root,"FlowData_2022ROC_EffectsAnalysis_CalendarYear.csv"),row.names=F)
write.csv(data_OMR_Dec_Jan,file.path(output_root,"FlowData_2022ROC_EffectsAnalysis_WaterYear.csv"),row.names=F)
