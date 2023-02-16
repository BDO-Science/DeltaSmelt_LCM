#####################################################################.
####### Code purpose: Estimate relationships between X2 and salinity; predict new salinity from X2 mgmt
####### Code developed by Compass Resource Management, for use in CSAMP Delta Smelt Structured Decision Making process.
####### Last updated: 2 Feb 2023
#####################################################################.
setwd("C:/Users/bcraw/Compass Resource Management Ltd/CSAMP Delta Smelt SDM - Documents/3. Mgmt Actions/6. Summer-Fall Flow and X2/IBMR Modeling/Data analysis/X2_salinity_analysis")
rm(list=ls())

library("plyr")
library("dplyr")
library("reshape")
library("stringr")
library("lubridate")
library("tidyr")
library("ggplot2")
library("gridExtra")
library("grid")
library("AICcmodavg")

#######################################################################
#Import and create X2 dataframe
X2.IBMR <- read.csv("X2_baseline.csv")
X2.dat <- melt(X2.IBMR,id.var=c("Year"))
colnames(X2.dat)[2:3] <- c("Month","X2")
levels(X2.dat$Month) <- 1:12
X2.dat$Month <- as.integer(X2.dat$Month)


#Create full dataframe (years, months, subregions)
X2.sal.dat <- expand.grid(Year = 1995:2014, Month = 1:12, Region = c("Confluence","East Delta","Lower Sacramento","Lower San Joaquin", "NE Suisun","NW Suisun",        
                                                                      "Sacramento River","SE Suisun","South Delta","Suisun Marsh","SW Suisun","Yolo Bypass")) # Create dataframe
regions <- levels(X2.sal.dat$Region)
X2.sal.dat$Region <- factor(X2.sal.dat$Region, levels=regions[c(12,7,9,2,3,4,1,8,5,10,11,6)]) # Reorder region names
X2.sal.dat <- left_join(x=X2.sal.dat, y=X2.dat[,c(1:ncol(X2.dat))], by=c("Year","Month")) # Join X2 data

sal_base <- read.csv("sal_baselineIBMR.csv", header=F) # From IBMR data. dimensions: 12 subregions (rows), 240 year-months from Jan 1995 to Dec 2014 (columns)

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

#X2.sal.dat$Month <- as.factor(X2.sal.dat$Month)
X2.sal.dat$Month <- X2.sal.dat$Month - 6.5
range(X2.sal.dat$Month)
sum(is.na(X2.sal.dat$sal))

# Fit regression models: X2
X2.mean <- mean(X2.sal.dat$X2)
X2.sd <- sd(X2.sal.dat$X2)
X2.sal.dat$X2 <- (X2.sal.dat$X2-X2.mean)/X2.sd

### Have to run through above lines of code before evaluating any action ###
#############################################################################

#X2.sal.dat$X2 <- scale(X2.sal.dat$X2)
fit0 <- glm(sal ~ 1,data=X2.sal.dat,family=Gamma(link="log"))
fit1 <- glm(sal ~ X2,data=X2.sal.dat,family=Gamma(link="log"))
fit2 <- glm(sal ~ X2 + Region,data=X2.sal.dat,family=Gamma(link="log"))
fit3 <- glm(sal ~ X2 + Month + Region,data=X2.sal.dat,family=Gamma(link="log"))
fit4 <- glm(sal ~ Month + Region*X2,data=X2.sal.dat,family=Gamma(link="log"))
fit5 <- glm(sal ~ Month + I(Month^2) + Region*X2 + I(X2^2),data=X2.sal.dat,family=Gamma(link="log"))
summary(fit5) # display results
lrModels <- list(fit0, fit1, fit2, fit3, fit4, fit5)
lrNames <- c("Null","fit1","fit2","fit3","fit4","fit5")
aicWt <- aictab(cand.set=lrModels, modnames=lrNames, sort=TRUE, c.hat=1, second.ord=TRUE)

backwards1 <- step(fit4,trace=1)

1-(fit5$deviance/fit5$null.deviance)

# Assess fit of best model
plot(fit5)
car::vif(fit5)
jpeg(file="Fig_X2_sal_mod_diagnostics.jpg",width = 8, height = 8, units = "in", res=600)
par(mfrow = c(2,2))
plot(fit5)
dev.off()
# Save best model
saveRDS(fit5,"model_sal_X2.Rdata")
# use best model (fit5) to predict salinity values
X2.sal.dat$sal.pred <- predict(fit5, newdata=X2.sal.dat, type="response")
X2.sal.dat$Month <- X2.sal.dat$Month + 6.5

# Plot X2 vs. salinity observations and predictions by subregion
X2.sal.dat$sal.pred[which(X2.sal.dat$sal.pred>max(X2.sal.dat$sal, na.rm=T))] <- NA
X2.sal.dat$X2 <- (X2.sal.dat$X2*11.06966)+73.82204
p.sal <- ggplot(X2.sal.dat, aes(x=X2,y=sal)) + 
  geom_point(size=0.4, alpha=0.5) +
  geom_line(data=X2.sal.dat, aes(x = X2, y = sal.pred, group=Month), size=0.7, color="red", alpha=0.5) +
  labs(x= "X2",y="Salinity (PSU)") + 
  scale_x_continuous(breaks=c(seq(44,89,length.out=6))) +
  theme_classic() + theme(panel.background = element_rect(color="black"),
                          axis.text.y=element_text(size=11),axis.text.x=element_text(size=11),axis.title=element_text(size=12,face="bold"),
                          legend.position=c(0.1,0.9),legend.title=element_blank(), plot.margin=unit(c(0.3,0.6,0.4,0.4),"cm"))
p.sal.sub <- p.sal + facet_wrap(~Region)
jpeg(file="Fig_X2 vs salinity_subregionIBMR.jpg",width = 8, height = 8, units = "in", res=600)
p.sal.sub
dev.off()


#######################################################################
# Estimating relationships salinity effects from X2 changes in other actions
#######################################################################
# Predict changes from actions/scenarios that change X2
salX2mod <- readRDS("Model_sal_X2.Rdata") # Load best salinity-X2 model
sal_base <- read.csv("sal_baselineIBMR.csv", header=F) # From IBMR data. dimensions: 12 subregions (rows), 240 year-months from Jan 1995 to Dec 2014 (columns)
years <- 20
months <- 12
n.strata <- 12
region.names <- levels(X2.sal.dat$Region)
### Get list of X2 actions/scenarios
file.list <- dir("C:/Users/bcraw/Compass Resource Management Ltd/CSAMP Delta Smelt SDM - Documents/3. Mgmt Actions/6. Summer-Fall Flow and X2/IBMR Modeling/Data analysis/X2_salinity_analysis/X2_inputs")

### The remaining code chunk predicts new salinity values for each X2 mgmt scenario. Need to specify the X2 scenario, one at a time in the next line, and run the remaining lines of code for each X2 scenario
scenario <- file.list[1]
X2_mgmt <- read.csv(paste0("X2_inputs/",scenario), header=F)   # dimensions: 20 years (rows: 1995:2014), 12 months (columns)

X2_mgmt[!is.na(X2_mgmt)] <- (X2_mgmt[!is.na(X2_mgmt)]-X2.mean)/X2.sd # Rescale X2

sal_mgmt <- sal_base #Create storage space for salinity results

# Check if X2 mgmt occurs for a given month and year (if not, no values are adjusted). 
# Use model to predict new salinity values)
for (s in 1:n.strata){
  for (y in 1:years){
    for (m in 1:months){ 
      sal_mgmt[s,(12*(y-1)+m)] <- if(is.na(X2_mgmt[y,m])) {NA} else { 
        predict(salX2mod, data.frame(Month=m-6.5, Region=region.names[s], X2=X2_mgmt[y,m]), type="response")}
    }
  }
}

# Replace any salinity values predicted beyond observed range with maximum observed salinity
sal_mgmt[sal_mgmt>max(sal_base,na.rm=T)] <- max(sal_base,na.rm=T)
# Save mgmt effects for IBMR data inputs
write.csv(sal_mgmt, paste0("sal_outputs/sal_mgmt_",str_sub(scenario,4)))
          
          


