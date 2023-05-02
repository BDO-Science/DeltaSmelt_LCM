### Create clean data files from the EMP Zooplankton data:
# K. Newman, 20 July 2011
# 15 Nov 2011
# 21 Nov 2011- looking at life stages within some taxa (Eurytemora affinis)
# 30 Nov 2011- using the "standardized" space-time series plotting
# 14 Jan 2012- bring in the 20 mm CB net samples
# 23 Mar 2012- extracting Physical data from Zoop surveys for Pete Smith
# 10 Aug 2012- creating clean file from CB net, including Lat and Long
# 16 Aug 2012- adding copepodids and a few others
# 12 Jul 2013- adding AVERNAL (to match Rose, Kimmerer, Edwards, + Bennett paper)
# 13 Jul 2013- just using the Rose et al taxa and
#             converting counts to Biomass using Wim's values
# 21 Mar 2014- re-examining the Biomass weights...


# Apr 8 - switching to using log of the geometric mean (as advised by Steve Slater,
#			4/8/13 email) instead of arithmetic mean
# Apr 9 - adding in EMP's SINOCAL with a weight of 3.4 (per Steve Slater's advice
#			and Wim's "bug weights" spreadsheet
# Apr 15 - changing how I deal with limnoithona - use the actual limnoithona fields
#			instead of the "Tot.limno" field that was apparently added (?)
# Aug 5 - change cleaning order: use taxon cut-off years to replace 0's with NA's
#           in individual data sets first, then merge Zooplankton and Mysid
# Feb 16, 16 - Using updated data sets through 2014.
# v6(lm): Updating data through 2016. The previous CB file was CPUV and I used 
#         carbon weights from Wim to calculte BPUV. Now I have a CB BPUV file 
#         from April, so use that instead.
# v7(WES): Updating data through 2019, using new files from Arthur Barros (CDFW)
# v8 (lp): Jan 4, 2021. Updating code to address issues with new data called in v7.
# -New data has different field names, different numbers of fields.
# -keep_col right below line 222 wants non-existent fields.
# -Number of rows of data by year is the same in new and old.
# -New data does not have Index==2 stations.
# -1994 NJ_BPUV_Mar0May0 is the problematic one from an LCM perspective.

###############################################################################
#setwd("//IFW8LODI-FS/Common/Statshare/DSLCM")
library(zoo)   # to use na.approx(); nice because it acts on data frames directly

## Load utility file:
predictor_root <- file.path("Input_Data","Predictor_Variables")
source(file="DataCleaner_Utility_v3.r")

## Display extra output (e.g. frequency tables) in R:
verbose <- TRUE

###################################################################################
## Save file names and zooplankton names/carbon weights:
station_filename <- "ZPStations.csv"


## EMPCBMatricesMASTMay2017.csv is a copy of the "1972-2014 CB BPUE Matrix" worksheet 
## in the file EMPCBMatricesMASTMay2017.xlsx:
#zoo_filename <- "1972-2019CBBPUEMatrix.csv"
zoo_filename <- "EMPCBMatricesMASTMay2017.csv"
zoo_taxon_filename <- "CB.taxon.cutoffs.csv"


## EMPMysidMatricesMASTMay2017.csv is a copy of the MysidBPUEMatrix1972-2016 worksheet 
## in the file EMPMysidMatricesMASTMay2017.xlsx:
#mysid_filename <- "1972-2019MysidBPUEMatrix.csv"
mysid_filename <- "EMPMysidMatricesMASTMay2017.csv"
mysid_taxon_filename <- "Mysids.taxon.cutoffs.csv"



## Select zooplankton species:
naup.names <- c("COPNAUP", "EURYNAUP", "OTHCOPNAUP", "PDIAPNAUP", "SINONAUP")

juv.names <- c("CALJUV", "EURYJUV", "OTHCALJUV", "PDIAPJUV", "SINOCALJUV", 
                "ASINEJUV", "ACARJUV", "DIAPTJUV", "TORTJUV")

adult.names <- c("EURYTEM", "OTHCALAD", "PDIAPFOR", "SINOCAL", "AVERNAL", "LIMNOSPP", 
                  "LIMNOSINE", "LIMNOTET")

clad.names <- c("BOSMINA", "DAPHNIA", "DIAPHAN", "OTHCLADO")



## Select mysid species:
mysid.names <- c("Hyperacanthomysis.longirostris","Neomysis.mercedis")



###################################################################################
## Read in station data frame and remove stations outside of the four main regions.
## Then read in zooplankton data frame and merge with station df.

## Create file names for reading in data:
station_file <- file.path(raw_root,"ZooMysid_raw",station_filename)

zoo_file <- file.path(raw_root,"ZooMysid_raw",zoo_filename)
zoo_taxon_file <- file.path(raw_root,"ZooMysid_raw",zoo_taxon_filename)

mysid_file <- file.path(raw_root,"ZooMysid_raw",mysid_filename)
mysid_taxon_file <- file.path(raw_root,"ZooMysid_raw",mysid_taxon_filename)


## Read in station data frame:
station_df_orig <- read.csv(station_file, header=TRUE, stringsAsFactors=FALSE)
station_df <- subset(station_df_orig, Region %in% region.list)

cat("Removing stations from the station data frame that are outside\nof the",
	"4 main regions (Far West, West, North, South) ...\n\n",
	nrow(station_df_orig) - nrow(station_df), "stations removed.\n\n")


## Read in zooplankton and mysid data frames:
zoo_df <- read.csv(zoo_file, header=TRUE, stringsAsFactors=FALSE)
zoo_df$Date <- as.Date(zoo_df$Date,"%m/%d/%Y")		# use name Date in both df's

mysid_df <- read.csv(mysid_file, header=TRUE, stringsAsFactors=FALSE)
mysid_df$Date <- as.Date(mysid_df$Date,"%m/%d/%Y")


## Use the Taxon Cutoff files to distinguish 0's from NA's. Between the start and end
## years, use 0. Outside of this range, use NA.
cat("Distinguishing between 0's and NA's ...\n")

# Zooplankton:
zoo_cutoffs <- read.csv(zoo_taxon_file, header=TRUE, stringsAsFactors=FALSE)
zoo_cutoffs[is.na(zoo_cutoffs$End.Year),"End.Year"] <- 2019

for(i in seq(nrow(zoo_cutoffs))) {
	species <- zoo_cutoffs[i,"Taxon"]
	start.yr <- zoo_cutoffs[i,"Start.Year"]
	end.yr <- zoo_cutoffs[i,"End.Year"]
	ok.years <- (start.yr <= zoo_df$Year & zoo_df$Year <= end.yr)
	cat("",species,"fraction of 0's kept =",sum(ok.years)/length(ok.years),"\n")
	zoo_df[!ok.years,species] <- NA  # everything before start year and after end year
}

# Mysid:
mysid_cutoffs <- read.csv(mysid_taxon_file, header=TRUE, stringsAsFactors=FALSE)
mysid_cutoffs[is.na(mysid_cutoffs$End.Year),"End.Year"] <- 2019

for(i in seq(nrow(mysid_cutoffs))) {
	species <- mysid_cutoffs[i,"Taxon.Name"]
	start.yr <- mysid_cutoffs[i,"Start.Year"]
	end.yr <- mysid_cutoffs[i,"End.Year"]
	ok.years <- (start.yr <= mysid_df$Year & mysid_df$Year <= end.yr)
	cat("",species,"fraction of 0's kept =",sum(ok.years)/length(ok.years),"\n")
	mysid_df[!ok.years,species] <- NA  # everything before start year and after end year
}
cat("\n")


## Merge by common columns Year, Date, Station, Survey, SurveyRep, SurveyCode, Index,  
## and TimeStart. Avoid using Chl.a because of rounding differences. It would be nice  
## if zoo_df had a TimeStart column ...
## Use the name ZooMysid_df for the combined data frame:

ZooMysid_df <- merge(zoo_df, mysid_df, by=c("Year","Date","Station","Survey","SurveyRep",
                      "SurveyCode","Index","TimeStart"), all=TRUE, sort=FALSE)
ZooMysid_df$Month <- factor(months(ZooMysid_df$Date), levels=month.name, ordered=TRUE)

station_match <- match(ZooMysid_df$Station, station_df$Station)
ZooMysid_df$Region <- station_df[station_match,"Region"]  # replaces the original Region field
ZooMysid_df$Lat <- station_df[station_match,"lat.decimal"]
ZooMysid_df$Lon <- -station_df[station_match,"long.decimal"]

## Per the recommendation in the document "ReadMeZooplanktonStudyMatricesJune2015.doc",
## only use core/index stations 1 - 2 and surveys 3 - 11 for long-term analysis. And maybe
## remove data prior to 1974 in order to keep core station 2.
#ZooMysid_df <- subset(ZooMysid_df, Index %in% 1:2 & Survey %in% 3:11 & Year >= 1974)

# Trying this: include all for comparison
#ZooMysid_df <- subset(ZooMysid_df, Index %in% 0:2 & Survey %in% 3:11 & Year >= 1974)

###stop and choose above about if Index =0 or not is included.

## Sort by date and station:
ZooMysid_df <- ZooMysid_df[order(ZooMysid_df$Date, ZooMysid_df$Station), ]


## Check that there is only one record per date-station:
tempDS <- with(ZooMysid_df, paste(Date,Station,sep="."))
if(length(unique(tempDS)) < length(tempDS)) {
	cat("Warning: multiple records per date-station.\n\n")
} else {
	rownames(ZooMysid_df) <- tempDS
}


## Delete records without station code matches in the Station Coordinate File:
## This will automatically remove any records outside of the 4 main regions.
ok <- ZooMysid_df$Station %in% station_df$Station
not.ok <- !ok

if(sum(not.ok) > 0){
	## Investigate records without matches:
	cat("Deleting records without station code matches in the station data frame ...\n")

	if(verbose) {
		my_df <- as.data.frame(table(ZooMysid_df$Station[not.ok]))
		colnames(my_df) <- c("Station","Num_of_records_deleted")
		print(my_df)
	}
	cat("\n")
	ZooMysid_df <- ZooMysid_df[ok,]
}


## Zooplankton: add carbon weights across species to get "total" BPUV measures:
cat("Calculating zooplankton and mysid 'total' BPUV measures by Date-Station ...\n\n")

naupBPUV <- apply(ZooMysid_df[ ,naup.names], 1, sum.na)
juvBPUV <- apply(ZooMysid_df[ ,juv.names], 1, sum.na)
adultBPUV <- apply(ZooMysid_df[ ,adult.names], 1, sum.na)
cladBPUV <- apply(ZooMysid_df[ ,clad.names], 1, sum.na)

ZooMysid_df$NJ_BPUV <- naupBPUV + juvBPUV
ZooMysid_df$JA_BPUV <- juvBPUV + adultBPUV
ZooMysid_df$JAC_BPUV <- juvBPUV + adultBPUV + cladBPUV
ZooMysid_df$NJAC_BPUV <- naupBPUV + juvBPUV + adultBPUV + cladBPUV
ZooMysid_df$AC_BPUV <- adultBPUV + cladBPUV


## Mysid: values are already carbon-weighted, but add to get a "total" BPUV measure:
ZooMysid_df$M_BPUV <- apply(ZooMysid_df[ ,mysid.names], 1, sum.na)
ZooMysid_df$JACM_BPUV <- with(ZooMysid_df, JAC_BPUV + M_BPUV)
ZooMysid_df$NJACM_BPUV <- with(ZooMysid_df, NJAC_BPUV + M_BPUV)
ZooMysid_df$ACM_BPUV <- with(ZooMysid_df, AC_BPUV + M_BPUV)



## Keep select columns:
cat("Keeping select columns, and removing the rest ...\n\n")
keep_col <- c("Year","Date","Month","Survey","SurveyRep","TimeStart","Station",
							"Region","Lat","Lon","Secchi.x","Chl.a.x","TempSurf.x",
							naup.names,juv.names,adult.names,clad.names,
							"NJ_BPUV","JA_BPUV","JAC_BPUV","NJAC_BPUV",
							mysid.names,"M_BPUV","JACM_BPUV","NJACM_BPUV","ACM_BPUV")
max(abs(ZooMysid_df$TempSurf.x-ZooMysid_df$TempSurf.y),na.rm=T)
max(abs(ZooMysid_df$Chl.a.x-ZooMysid_df$Chl.a.y),na.rm=T)

# keep_col <- c("Year","Date","Month","Survey","SurveyRep","TimeStart","Station",
#               "Region","Lat","Lon","Secchi.x","Chl_a.x","Temperature.x",
#               naup.names,juv.names,adult.names,clad.names,
#               "NJ_BPUV","JA_BPUV","JAC_BPUV","NJAC_BPUV",
#               mysid.names,"M_BPUV","JACM_BPUV","NJACM_BPUV","ACM_BPUV")
#max(abs(ZooMysid_df$Temperature.x-ZooMysid_df$Temperature.y),na.rm=T)
#max(abs(ZooMysid_df$Chl_a.x-ZooMysid_df$Chl_a.y),na.rm=T)
setdiff(keep_col,names(ZooMysid_df))

ZooMysid_df <- ZooMysid_df[ ,keep_col]

names(ZooMysid_df) <- sub(".x","",names(ZooMysid_df), fixed=TRUE)


## Sort:
ZooMysid_df <- with(ZooMysid_df, ZooMysid_df[order(Date,Survey,Station,SurveyRep), ])


## Display some tables, etc.:
if(verbose){
	table(ZooMysid_df$Year, ZooMysid_df$Month)
	table(ZooMysid_df$Survey, ZooMysid_df$Month)
	names(ZooMysid_df)
	dim(ZooMysid_df)
}


###################################################################################

## Save ZooMysid DS files:
ZooMysid_74_19_df <- ZooMysid_df
write.csv(ZooMysid_74_19_df, file=file.path(clean_food_root,"ZooMysid_74_19_df.csv"), 
					row.names=FALSE)

ZooMysid_df<-read.csv("ZooMysid_74_19_df.csv")

##############################################################################
##############################################################################
## Use linear interpolation model ...

BPUV <- c("NJ_BPUV","JA_BPUV","JAC_BPUV","NJAC_BPUV","M_BPUV","JACM_BPUV","NJACM_BPUV",
          "ACM_BPUV")
ZooMysid_YMR <- split.average(ZooMysid_df,splitBy="ymr",average="median",fields=BPUV)
# other average options: "mean", "geo", "log.geo", "sum"


## Use this for imputation and plotting:
ZooMysid_YMR[ ,paste0(BPUV,"_orig")] <- ZooMysid_YMR[ ,BPUV]
ZooMysid_YMR[ ,paste0(BPUV,"_col")] <- regCol(ZooMysid_YMR$Region)

for(r in region.list) {
	rSub <- subset(ZooMysid_YMR, Region==r)
	Time <- 1:nrow(rSub)	# assuming df contains all possible years and months

	for(s in BPUV) {
		# Use common rownames between ZooMysid_YMR and rSub to impute:
		isNA <- is.na(rSub[,s])
		rnames <- rownames(rSub)[isNA]

		impute_all <- na.approx(rSub[,s],x=Time,xout=Time,na.rm=FALSE)
		ZooMysid_YMR[rnames,s] <- impute_all[isNA]
		ZooMysid_YMR[rnames,paste0(s,"_col")] <- "orange"
	}
}


## Dates relative to January 1, 1970:
xrange <- c(as.numeric(as.Date("1970-01-01")), as.numeric(as.Date("2011-01-01")))
mySeq <- seq(xrange[1], xrange[2], by=365)

plot_BPUV <- function(field, main) {
	for(r in region.list) {
		rSub <- subset(ZooMysid_YMR, Region == r)

		if(!plot3) dev.new(width=13,height=4.5,units="in")
		yrange <- range(rSub[ ,field], na.rm=TRUE)
		par(mfrow=c(1,1),mar=c(4.5,4.5,3,2),oma=c(0,0,0,1),cex.axis=1.5,cex.lab=1.5,
				cex.main=1.5)
		plot(rSub$pDate, rSub[,field], xlim=xrange, ylim=yrange, xaxt="n",
				type="b", lwd=1, main=paste(main,r), xlab="Date",
				ylab="Median by YMR", col=regCol(r))
		points(rSub$pDate, rSub[,field], col=rSub[,paste0(field,"_col")], lwd=2)
		axis(side=1, at=mySeq, labels=as.Date(mySeq, origin="1970-01-01"))
		legend("topright",region.list,col=regCol(region.list), lwd=rep(2,4), lty=rep(1,4),
				cex=1.4, pch=rep(1,4))
		legend("bottomleft","interpolated",col="orange",pch=1, cex=1.4, lwd=2, lty=NA)
	}
}

plot3 <- FALSE
if(plot3) pdf("ZooMysid_median_lin_interp.pdf", width=18, height=10, onefile=TRUE)

plot_BPUV("NJ_BPUV","nauplii + juvenile BPUV")

plot_BPUV("JA_BPUV","juvenile + adult BPUV")
plot_BPUV("JAC_BPUV","juvenile + adult + cladocerans BPUV")

plot_BPUV("JACM_BPUV","juvenile + adult + clad + mysid BPUV")

if(plot3) dev.off()


##############################################################################

## Save ZooMysid YMR files. First remove extra fields:
ZooMysid_74_19_df_median <- ZooMysid_YMR[ ,c("Year","Month","Region",BPUV)]

write.csv(ZooMysid_74_19_df_median, 
					file=file.path(clean_food_root,"ZooMysid_74_19_df_median.csv"),
					row.names=FALSE)

##############################################################################
## Output info on when/where imputation was done (for documentation purposes).
regAbbrev <- function(x) {
	abbVec <- c("Far West"="FW","West"="W","South"="S","North"="N")
	return(abbVec[as.character(x)])
}

makeImputeTable <- function(useDF,field) {
	before <- useDF[ ,paste0(field,"_orig")]
	after <- useDF[ ,field]

	mySub <- useDF[is.na(before) & !is.na(after), ]
	myVec <- with(mySub, ULS(Region, list(Month,Year), drop=TRUE, FUN=function(x) {
					paste(regAbbrev(x), collapse=",") } ))
          
  # months_f <- factor(extractName(names(myVec),1), levels=month.name, ordered=TRUE)
  # months_f <- droplevels(months_f)
	myDF <- data.frame("Month"=extractName(names(myVec),1),
                      "Year"=extractName(names(myVec),2),
                      "Reg"=myVec, stringsAsFactors = FALSE)
	myDF[myDF$Reg=="FW,W,N,S","Reg"] <- "All"
	finalDF <- reshape(myDF, v.names="Reg", idvar="Year", timevar="Month",
                      direction="wide")
	finalDF[is.na(finalDF)] <- 0

  names(finalDF) <- sub("Reg.","",names(finalDF),fixed=TRUE)
  my_months <- names(finalDF)[names(finalDF) != "Year"]
  my_months_f <- factor(my_months, levels=month.name, ordered=TRUE)
  my_months_f <- droplevels(my_months_f)
  
  finalDF <- finalDF[ ,c("Year", levels(my_months_f))]
	print(finalDF)
}

a <- makeImputeTable(useDF=ZooMysid_YMR, field="NJ_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="JA_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="JAC_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="NJAC_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="M_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="JACM_BPUV")
makeImputeTable(useDF=ZooMysid_YMR, field="NJACM_BPUV")


#write.csv(a,paste0(clean_food_R,"zoop_impute_table.csv"))


# # ##############################################################################
# # ##############################################################################
# # ## Use ZooMysid_YMR to create an array that can be used by the ADMB input data
# # ##   creation file easily. These data are also normalized (de-meaned and
# # ##   scaled):

# # BPUV <- c("NJ_BPUV","NJAC_BPUV","JACM_BPUV","NJACM_BPUV")
# # useYears <- as.character(with(ZooMysid_YMR, seq(min(Year),max(Year),by=1)))

# # # In ZooMysid_YMR, Month changes fastest, then Year, then Region, so create
# # #   the array like this: [month,year,region,BPUV measure]
# # ZooMysid_array <- array(0,dim=c(12,length(useYears),4,length(BPUV)),
												# # dimnames=list("Month"=month.name,"Year"=useYears,
																			# # "Region"=region.list,"Type"=BPUV))
# # for(sp in BPUV) {
	# # ZooMysid_array[ , , ,sp] <- ZooMysid_YMR[,sp]
# # }

# # ## Permute to this order: Month, Region, Year, Composition
# # ZooMysid_array <- aperm(ZooMysid_array, c("Month","Region","Year","Type"))

# # # Normalize using the mean and sd calculated over ALL the factors (year,
# # #   month, region, and BPUV measure)
# # myMean <- mean(ZooMysid_array, na.rm=TRUE)
# # mySD <- sd(ZooMysid_array, na.rm=TRUE)
# # ZooMysid_array_norm <- (ZooMysid_array - myMean)/mySD

# # range(ZooMysid_array_norm, na.rm=TRUE)

# # # Also try just scaling it, so that food values always remain >= 0:
# # mySD <- sd(ZooMysid_array, na.rm=TRUE)
# # ZooMysid_array_scaled <- ZooMysid_array/mySD

# # range(ZooMysid_array_scaled, na.rm=TRUE)



# # ##############################################################################

# # # # ## Save zooplankton ADMB array files:
# # # # ZooMysid_74_16_array_norm_median <- ZooMysid_array_norm
# # # # dump("ZooMysid_74_16_array_norm_median",file=paste0(clean_food_R,"ZooMysid_74_16_array_norm_median.R"))
# # # # # write.csv(ZooMysid_74_16_df_array_norm_median,paste0(clean_food_root,"ZooMysid_74_16_df_array_norm_median.csv"),
			# # # # # row.names=FALSE)

# # # # ZooMysid_74_16_array_scaled_median <- ZooMysid_array_scaled
# # # # # write.csv(ZooMysid_74_16_df_array_scaled_median,paste0(clean_food_root,"ZooMysid_74_16_df_array_scaled_median.csv"),
# # # # dump("ZooMysid_74_16_array_scaled_median",file=paste0(clean_food_R,"ZooMysid_74_16_array_scaled_median.R"))
			# # # # # row.names=FALSE)

# # ##############################################################################




