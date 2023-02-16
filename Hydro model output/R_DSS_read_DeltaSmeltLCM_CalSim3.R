# Identify your working directory for saving outputs of interest
root <- "~/GitHub/DeltaSmelt_LCM_test"
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
NAA_data <- dss_data_pull_LCME()
#D1641_data <- dss_data_pull_LCME(dss_input="")


#Export DSS output files for Delta Smelt LCM model input
write.csv(NAA_data,file.path(output_root,"NAA_CalSim3_data.csv"),row.names=F)

