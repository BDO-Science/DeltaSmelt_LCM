### Utility Functions for Delta Smelt Modeling
#
# Written by K. Newman, Nov 2011 - Jun 2013
# Functions added by L. Mitchell, Jan 2014
#
# July 28, 2014: Fixed typo in logit (formerly said log(1/(1-x))
###################################################################################

## Define directory:

abund_root <- file.path("Input_Data","Population_Abundance_Estimates","Code_And_Data",
												"0_Data_Cleaning")

raw_root <- file.path(predictor_root,"Data_Raw")
clean_root <- file.path(predictor_root,"Output")

raw_cdec_root <- file.path(raw_root,"CDEC_Temperature_raw")
raw_iss_root <- file.path(raw_root,"ISS_Beachseine_raw")
raw_length_root <- file.path(raw_root,"Length_raw")
raw_mysid_root <- file.path(raw_root,"Mysid_raw")
raw_phys_root <- file.path(raw_root,"PhysicalVariable_raw")
raw_salv_root <- file.path(raw_root,"Salvage_raw")
raw_zoop_root <- file.path(raw_root,"Zooplankton_raw")
raw_fish_root <- file.path(abund_root,"Data_Raw","Catch_Data_Sets")
raw_station_root <- file.path(abund_root,"Data_Raw","Station_Data_Sets")

clean_cdec_root <- file.path(clean_root,"CDEC_Temperature")
clean_length_root <- file.path(clean_root,"FishSurveyLength")
clean_food_root <- file.path(clean_root,"Food")
clean_iss_root <- file.path(clean_root,"ISS_Beachseine")
clean_LSZ_root <- file.path(clean_root,'LSZ_MacWilliams')
clean_physical_root <- file.path(clean_root,"PhysicalVariable")
clean_salv_root <- file.path(clean_root,"Salvage")
clean_spawn_root <- file.path(clean_root,"Spawning")
clean_fish_root <- file.path(abund_root,"Output")
# # # clean_entrain_R <- file.path(clean_root,"Entrainment_R_objects")

###################################################################################

month.list <- month.name
region.list <- c("Far West","West","North","South")


## 'Permanent' region colors for consistent plotting. Region is sometimes
## a character, sometimes a factor, so always convert to character first.
## The default argument is all regions; this facilitates making legends.
regCol <- function(x=region.list) {
	colVec <- c("Far West"="black","West"="red","South"="blue","North"="green")
	return(colVec[as.character(x)])
}

expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))


mytitle <- function(x) mtext(x,side=3,outer=TRUE)



## A function that "strips" out "compacted" location coordinates
## and translates the results into a "decimal" measure:
wgs.converter <- function(x) {
	degrees <- trunc(x)
	minutes <- trunc(100*(x-degrees))
	temp    <- 100*x - trunc(100*x)
	seconds <- (x*100-trunc(x*100))*60

	out <- degrees+minutes/60+seconds/3600
	return(out)
}

## Convert conductivity to salinity, and vice versa:
conductivity.to.salinity <- function(x)  -100*log(1-x/178500)
salinity.to.conductivity <- function(y) 178500*(1-exp(-0.01*y))


ratio.calc <- function(x,y,out.opt="ratio") {
	# 17 Jan 2012
	#--note allowing missing values to differ between x and y
	#---also sloppy with std error calcs if uneven number of obs'ns
	n <- length(x)
	if(length(y) != n) stop("error in ratio.sd, unequal x and y lengths")
	x.bar <- mean(x,na.rm=TRUE)
	y.bar <- mean(y,na.rm=TRUE)
	ratio <- x.bar/y.bar
	v.ratio <- (1/y.bar)^2*var(x,na.rm=TRUE)/n + (x.bar^2/y.bar^4)*var(y,na.rm=TRUE)/n
	s.ratio <- sqrt(v.ratio)
	if(out.opt=="ratio") {
		out <- list(ratio=ratio)
	} else {
		if(out.opt=="se") {
			out <- list(s.ratio=s.ratio)
		} else  out <- list(ratio=ratio,s.ratio=s.ratio)
	}
	return(out)
}

#-- Cochram p 155, eq'n 6.13
ratio.se <- function(y,x) {
 n <- length(y)
 ratio <- sum(y)/sum(x)
 s.y.x <- cov(x,y)
 ratio.var <- (1/(n*mean(x)^2))*(var(y) + ratio^2*var(x) - 2*ratio*s.y.x)
 out <- sqrt(ratio.var)
 return(out)
}


#--list ratio function
ratio.list <- function(y,x,z,trim.opt=0) {
 y.split   <- split(y,z)
 x.split   <- split(x,z)
 ratio     <- unlist(lapply(y.split,mean,trim=trim.opt))/
              unlist(lapply(x.split,mean,trim=trim.opt))

 n         <- length(ratio)
 temp      <- match(z,names(ratio))
 error     <- y - ratio[temp]*x
 sq.error  <- (error)^2
 ssq.error <- unlist(lapply(split(sq.error,z),sum))
 n.vec     <- unlist(lapply(y.split,length))
 no.rep    <- n.vec==1
 ratio.var <- rep(NA,n)
 ratio.var[!no.rep] <- ssq.error[!no.rep]/(n.vec[!no.rep]-1)
 ratio.se  <- sqrt(ratio.var)
 out       <- list(ratio=ratio,ratio.se=ratio.se)
 return(out)
}

#tc <- c(10,8,0,1,12,2)
#tv <- c(2000,2300,2200,1990,2210,2309)
#tz <- c(1,1,1,1,2,2)
#ratio.list(tc,tv,tz)

#tr <- sum(tc[1:4])/sum(tv[1:4])
#sqrt(sum((tc[1:4]-tr*tv[1:4])^2)/3)
#tr <- sum(tc[5:6])/sum(tv[5:6])
#sqrt(sum((tc[5:5]-tr*tv[5:6])^2)/1)

decimal.to.degrees.converter <- function(x) {
 degrees <- trunc(x)
 rem     <- x-degrees
 seconds <- rem*3600
 minutes <- trunc(seconds/60)
 seconds <- seconds-minutes*60
  out     <- list(degrees=degrees,minutes=minutes,seconds=seconds)
 return(out)
}

space.time.plots <- function(data.set,group.list,region.list,data.region,
 data.year.month,yr.min,yr.max,data.label="",
  robust=FALSE,robust.label="",fn=NULL,pdf.it=FALSE,cex.opt=0.7,draw.opt="b") {

  year.list <- yr.min:yr.max
  year.month.region.group.array  <-
      array(NA,dim=c(yr.max-yr.min+1,12,length(region.list),length(group.list)),
      dimnames=list(yr.min:yr.max,month.list,region.list,group.list))

  if(pdf.it) pdf(file=fn)
  par(mfrow=c(1,1),ask=FALSE)
  my.cex <- cex.opt
  for(i in 1:length(group.list)) {
   temp.counts     <- data.set[,group.list[i]]
   for(j in 1:4) {
    region.opt   <-  region.list[j]
    temp.grps    <- split(temp.counts[data.region==region.opt],
      data.year.month[data.region==region.opt])
   if(robust) {
       temp.central <- unlist(lapply(temp.grps,median,na.rm=TRUE))
     } else {
       temp.central <- unlist(lapply(temp.grps,mean,na.rm=TRUE))
     }
   if(j==1) {
        ylim.lower <- min(temp.central,na.rm=TRUE)
        ylim.upper <- max(temp.central,na.rm=TRUE)
     } else {
        ylim.lower <- min(c(ylim.lower,temp.central),na.rm=TRUE)
        ylim.upper <- max(c(ylim.upper,temp.central),na.rm=TRUE)
     }
  }
  cat(group.list[i],ylim.lower,ylim.upper,"\n")

 for(j in 1:4) {
     year.month.mat <- matrix(NA,nrow=yr.max-yr.min+1,ncol=12,
      dimnames=list(yr.min:yr.max,month.list))
    region.opt   <- region.list[j]
    temp.grps    <- split(temp.counts[data.region==region.opt],
      data.year.month[data.region==region.opt])

    if(robust) {
     temp.central <-  unlist(lapply(temp.grps,median,na.rm=TRUE) )
     }  else {
     temp.central <-  unlist(lapply(temp.grps,mean,na.rm=TRUE) )
    }
    x.year       <- substr(names(temp.central),start=1,stop=4)
    x.month      <- substring(names(temp.central),first=6)
    row.pos      <- match(x.year,dimnames(year.month.mat)[[1]])
    col.pos      <- match(x.month,dimnames(year.month.mat)[[2]])
    year.month.mat[cbind(row.pos,col.pos)] <- temp.central
    temp.central <- ts(data=as.vector(t(year.month.mat)),
      start=as.numeric(yr.min),frequency=12)
   if(j==1) {
     plot(temp.central,ylab=paste(robust.label,group.list[i]),
       main=paste(data.label,robust.label,"\n",group.list[i]),
         ylim=c(ylim.lower,ylim.upper),
       type='n',cex=my.cex)
    }
    points(temp.central,col=j,pch=j,cex=my.cex,type=draw.opt)
    cat("done with group",i,"region",j,"\n")
    year.month.region.group.array[,,j,i] <- year.month.mat

   }
  legend(year.list[9],0.8*ylim.upper,legend=region.list,pch=1:4,col=1:4)
 }
 par(mfrow=c(1,1),ask=FALSE)
 if(pdf.it) dev.off()
 return(year.month.region.group.array)
}

catch.ts.plot <- function(y,year.subset,date.subset,chop.quantile=0.95,label) {
  plot(date.subset,y,ylim=c(0,quantile(y,prob=chop.quantile )),
     type='n',main=label,xlab="Year",ylab="",axes=FALSE)
  axis(side=2,tick=TRUE)
  #axis(side=1,tick=FALSE)
  uniq.years <- unique(year.subset)
 for(i in 1:length(uniq.years)) {
   interval <- year.subset == uniq.years[i]
   lines(loess(y[interval]~date.subset[interval]))
   axis(side=1,at=date.subset[interval][1],labels=uniq.years[i]%%100)
   #lines(date.subset[interval],y[interval],type='b')
 }
}

#--- Code to create flat file summarized over Year-Month-Region
space.time.file.creation <- function(data.set,region.list,env.list,species.list,data.region,
 data.year.month,yr.min,yr.max,data.label="",
  robust=FALSE,robust.label="",fn=NULL) {

  num.regions <- length(region.list)
  num.species <- length(species.list)
  num.rows <- (yr.max-yr.min +1)*12*num.regions
  num.cols <- 1+1+1+length(env.list)+1+ num.species
  out <- matrix(NA,nrow=num.rows,ncol=num.cols,
   dimnames=list(NULL,
   c("Year","Month","Region",env.list,"Volume",
    paste(species.list,rep("n",num.species),rep("sd",num.species)))))

 return(out)
 }

#--- UTM to Lat-Long
utm.to.lat.long <- function(northing,easting,zone=10,cm.easting=500000) {
 #From www.uwgb.edu/dutchs/UsefulData/UTMFormulas.htm  (and spreadsheet)
 #i.e., this is purely a copy of Professor Steven Dutch's work
 #Ken Newman, USFWS, 5 December 2011

 if(zone>0) long0 <- 6*zone-183 else stop("Invalid zone value")

 k0 <- 0.9996
 a  <- 6378137           #equatorial radius
 b  <- 6356752.3142      #polar radius
 e  <- sqrt(1-b^2/a^2)    #eccentricity in earth's elliptical cross-section

 x  <-  cm.easting - easting
 y  <- northing

 # Meridional Arc
 M  <- y/k0

 #Footprint Latitude
 mu <- M/(a*(1-e^2/4-3*e^4/64-5*e^6/4^4 - 6*e^8/4^5))
 e1 <- (1-sqrt(1-e^2))/(1+sqrt(1-e^2))
 #J1 to J4 are apparently approximations, truncated?
 J1 <- (3*e1/2-27*e1^3/32)
 J2 <- (21*e1^2/16 - 55*e1^4/32)
 J3 <- 151*e1^3/96
 J4 <- 1097*e1^4/512
 fp <- mu + J1*sin(2*mu)+J2*sin(4*mu) + J3*sin(6*mu) + J4*sin(8*mu)

 #Latitude and Longitude
 eprime2 <- (e*a/b)^2
 C1 <- eprime2*cos(fp)^2
 T1 <- tan(fp)^2
 R1 <- a*(1-e^2)/(1-e^2*sin(fp)^2)^(3/2)
 N1 <- a/(1-e^2*sin(fp)^2)^(1/2)
 D <- x/(N1*k0)

 #Coefficients for Calculating Latitude
 Q1 <- N1*tan(fp)/R1
 Q2 <- D^2/2
 Q3 <- (5+3*T1+10*C1-4*C1^2-9*eprime2)*D^4/24
 Q4 <- (61+90*T1+298*C1+45*T1^2-3*C1^2-252*eprime2)*D^6/720

 #Coefficients for Calculating Longitude
 Q5 <- D
 Q6 <- (1+2*T1+C1)*D^3/6
 Q7 <- (5-2*C1+28*T1-3*C1^2+8*eprime2 + 24*T1^2)*D^5/120

 lat <- (fp-Q1*(Q2+Q3+Q4))*180/pi
 long <- long0 - ((Q5-Q6+Q7)/cos(fp))*180/pi
 out <- list(lat=lat,long=long)
 return(out)
}

#utm.to.lat.long(northing=4273583,easting=631504,zone=10,cm.easting=500000)
#  $lat  # [1] 38.6010039781
# $long  # [1] -121.489817545

#utm.to.lat.long(4386575,586762)
# $lat   # [1] 39.6245449326
# $long  # [1] -121.989096590

# Age 0 length cut-offs by month, from Steve Slater, 14 Dec 2011
age0.cutoffs <- c(40,50,50,50,50,60, 65, 70, 75, 80,80,80)
names(age0.cutoffs) <- month.list


linear.interpolate <- function(x) {
 #simple linear interpolation for missing values in a numerical vector,
 # find days bracketing missing days
 # calculate the line from bracket days' values
 n <- length(x)
 out <- x
 for(i in 1:n) {
  if(is.na(x[i])) {
   if(i > 1 & i < n) {
    before <- i-1; after <- i+1
    while(is.na(x[after]) & after <=n)  { after <- after + 1}
    slope <- (x[after]-x[before])/(after-before)
    intercept <- x[before];#cat(intercept,slope,"\n")
    x[i:(after-1)] <- out[i:(after-1)] <- intercept + slope*((i-before):(after-1-before))
    #out[i] <- mean(c(out[before],x[after]),na.rm=TRUE)
    } else {    #if 1st value is missing
      if(i==1) {
       j  <- 2
       while(is.na(x[j])) j <- j+1
       out[i] <- x[j]
     } else  { #if last value is missing
       out[n] <- out[n-1]
      }
    }
   }
  }
  return(out)
 }

 # test function
 #x <- c(10,8,NA,NA,3)
 #x <- c(20,10,12,NA,NA,NA,20)
 #linear.interpolate(x)

 sal.ec <- function (ec, temp=-100)
    {
    #    Calculates salinity from EC, assumed to be at 25C, or cond. at given temp
    #  ec is a vector of ec values in either mS/cm or uS/cm
    #  Converts EC in mS/cm to salinity.
    #     Algorithm taken from Matlab routines by CSIRO, based on a conductivity ratio Rt
    #    which is the ratio of sample conductivity to conductivity at S<-25, T<-15

    # SEAWATER Library
    # Version 2.0.1   22-Apr-1998
    #                 * Phil.Morgan@marine.csiro.au *

    #  The conductivity ratio of the sample at S, 25C, 0 pressure is calculated first from
    #  the conductivity at 25, 15, 0
    #     Using the algorithm in program
    # Algorithm: Salinity <- f ( Rt )
    #        Rt = COnd (Unknown S, T=25, P=0) / Cond (S=35, T=15, P=0)
    #            =  Cond (S, 25, 0) /  ( Cond (35, 15, 0)  *  ( Cond (35,25,0) / Cond (35,15,0) ) )
    #
    #
    # Extension for S < 2 due to Hill (1978) added May 2007 - from Alan Jassby
    # First step: cond (35,15,0) is "given" in sw_c3515.m in above program
    c3515 <- 42.914

    # Next rt is calculated as C(35,T,0)/C(35,15,0) from program sw_salrt.m
    # Eqn (3) p.7 Unesco.

    if (length(temp)>1 & length(temp) !=length(ec))
              stop ("Inputs must be same length or temp must be length 1")
    tt <- ifelse (temp < -10,  25, temp)
    c0 <-  0.6766097
    c1 <-  2.00564e-2
    c2 <-  1.104259e-4
    c3 <- -6.9698e-7
    c4 <-  1.0031e-9

    rt <- c0 + (c1 + (c2 + (c3 + c4 * tt) * tt) * tt) * tt #Cond ratio at 25C to that at 15 C

    # Then the ratio Rt is calculated (the conductivity ratio of the sample)
    # and salinity is  calculated from sw_sals.m

    if(max(na.omit(ec)) > 100) ec <- ec/1000
    Rt <- ec / ( c3515 * rt )


    a0 <-  0.0080
    a1 <- -0.1692
    a2 <- 25.3851
    a3 <- 14.0941
    a4 <- -7.0261
    a5 <-  2.7081

    b0 <-  0.0005
    b1 <- -0.0056
    b2 <- -0.0066
    b3 <- -0.0375
    b4 <-  0.0636
    b5 <- -0.0144

    k <-  0.0162
    x <- 400*Rt
    y <- 100*Rt

    Rtx <- sqrt(Rt)
    delT <- tt - 15
    ft <- delT / (1 + k*delT)
    delS <- ft * ( b0 + (b1 + (b2+ (b3 + (b4 + b5 * Rtx) * Rtx) * Rtx)* Rtx) * Rtx)
    pss78 <- a0 + (a1 + (a2 + (a3 + (a4 + a5 * Rtx) * Rtx) * Rtx) * Rtx) * Rtx
    pss78 <- pss78 + delS
    S <-  ifelse(pss78 <= 2, pss78, pss78 - a0/(1+1.5 * x + x^2) - ft*b0/(1 + y^.5 + y^1.5))
    S
}


#-- determining logit function coefficients given particular x and y values

logit.coef <- function(p.low,p.high,x.low,x.high) {
 y.low <-  log(p.low/(1-p.low))
 y.high <- log(p.high/(1-p.high))
 slope <- (y.high-y.low)/(x.high-x.low)
 intercept <- y.high-slope*x.high
 out <- list(b0=intercept,b1=slope)
 return(out)
}

#out <- logit.coef(p.low=0.1,p.high=0.9,x.low=8,x.high=10)

#x.seq <- seq(0,12,by=0.2)
#y.seq <- expit(out$b0,out$b1,x.seq)
#plot(x.seq,y.seq,type='l',xlab="Conductivity (1000s)",ylab="",main="Pr(move|Conductivity)")


#--- pairs panel function correlation with NA's ok

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

#---  Matrix version of convolution of rmultinom
rmultinom.convolution <- function(X,P) {
 # X is a matrix, n.components x n.particles, of n's
 # P is a 3-D array, n.components x n.components x n.particles, of probabilities
 #   where the probabilities sum to 1 by row in first two dimensions
 n.components <- dim(X)[[1]]
 n.outcomes   <- dim(P)[[2]]
 n.particles  <- dim(P)[[3]]

 output       <- matrix(data=0,nrow=n.components,ncol=n.particles)

 remainder.n <- X
 cumul.prob  <- matrix(data=0,nrow=n.components,ncol=n.particles)

 #run through the columns of the Probability matrix, using conditional
 #  binomials to create the binomaial
 for(i in 1:(n.components-1)) {
    unadj.prob <- P[,i,]                     #extract column from Prob matrix for next location
    adjust.prob <- unadj.prob/(1-cumul.prob) #adjust probability sequentially
    temp <- rbinom(n=n.components*n.particles,size=as.vector(remainder.n),
                prob=as.vector(adjust.prob))
    temp <- matrix(temp,nrow=n.components,ncol=n.particles)
    cumul.prob <- cumul.prob+unadj.prob
    remainder.n <- remainder.n - temp
    output[i,] <-  apply(temp,2,sum)
 }
 output[4,] <-  apply(remainder.n,2,sum)
 return(output)
}


# Steve Slater's Age Assignment based on Length and Month
dsm.age.key <- function(mths,lengths) {
 x <- c(40, 50, 50, 50, 50, 60, 65, 70, 75, 80, 80, 80)
 mth.match <- match(mths,1:12)
 ageclass <- rep(NA,length(lengths))
 ageclass[lengths < x[mth.match]] <- 0
 ageclass[lengths >= x[mth.match]] <- 1
 return(ageclass)
}

# x <- c(8,2,3,12,7,NA)
# y <- c(90,40,50,60,55,NA)
# dsm.age.key(mths=x,lengths=y)


wtd.mean.sd <- function(n.vec,y.vec,na.opt=TRUE) {
  y <- rep(y.vec,n.vec)
  y.bar <- mean(y,na.rm=na.opt)
  y.sd  <- sd(y,na.rm=na.opt)
  out <- list(y.bar=y.bar,y.sd=y.sd)
  return(out)
 }


#  12 Dec 2012, solving for desired logistic model coeff given
# 2 desired probabilities for 2 covariate values
# logit.coef <- function(p1,p2,x1,x2) {
# y1 <- log(p1/(1-p1))
# y2 <- log(p2/(1-p2))
# b1 <- (y1-y2)/(x1-x2)
# b0 <- y2-b1*x2
# out <- list(b0=b0,b1=b1)
# return(out)
# }


 #-- Zero inflated Poisson (ZIP) mle
 zip.mle <- function(theta,x) {
   # 9 Jan 2013  - 10 Jan 2013
   prob0  <- theta[1]
   lambda <- theta[2]
   n <- length(x)
   indicator.0 <- x==0
   log.like <- sum(log( prob0*indicator.0 + (1-prob0)*dpois(x,lambda)))
    of <- -log.like
   return(of)
   }

 testing <- FALSE
 if(testing) {
   n <- 500
   prob0   <- 0.3; lambda <- 7
   zeros   <- rbinom(n=1,size=n,prob=prob0)
   notzero <- rpois(n=(n-zeros),lambda=lambda)
   x       <- c(rep(0,zeros),notzero)

   initial.par <- c(prob0=0.5,lambda=7)
   out <- optim(par=initial.par,fn=zip.mle,x=x,method="L-BFGS-B",
    lower=c(0,0),upper=c(1,Inf))
   print(out$par)
}

  zip.offset.mle <- function(theta,y,myoffset,covariates) {
   # 9 Jan 2013
   #y is a vector of length n, the observations
   #myoffset is a vector of length n
   #covariates is an n x p matrix with p covariates

   prob0           <- theta[1]
   beta.vec        <- cbind(theta[-1])
   log.intensity   <- cbind(1,covariates) %*% beta.vec
   log.lambda      <- log(myoffset) + log.intensity
   lambda          <- exp(log.lambda)
   print(as.vector(lambda))
   indicator.0 <- y==0
   log.like <- sum(log(prob0*indicator.0 + (1-prob0)*dpois(y,lambda)))
   print(log.like)
   of <- -log.like
   return(of)
   }

if(testing) {
 n        <- 30
 prob0    <- 0.3
 myoffset <- rep(10,n)
 b0       <- 0.1
 b1       <- 1.5
 covariates <- runif(n,0,1)
 log.intensity  <- b0 + b1*covariates
 my.offset      <- 10
 lambda         <- my.offset * exp(log.intensity)

 zeros   <- rbinom(n=1,size=n,prob=prob0)
 notzero <- rpois(n=(n-zeros),lambda=lambda)
 y   <- c(rep(0,zeros),notzero)

  #how well can Poisson params be estimated?
 y.trunc      <- y[-c(1:zeros)]
 x.trunc      <- covariates[-c(1:zeros)]
 offset.trunc <- rep(my.offset,n-zeros)

 m.nonzero <- glm(y.trunc ~ x.trunc, offset=log(offset.trunc),family=poisson)
 print(cbind(c(b0,b1),coef(m.nonzero)))

 initial.par <- c(prob0=0.2,b0=0.5,b1=0.4)
 out <- optim(par=initial.par,fn=zip.offset.mle,y=y,myoffset=my.offset,
   covariates=covariates,method="L-BFGS-B",lower=c(0,-Inf,-Inf),upper=c(1,Inf,Inf))
 print(cbind(c(prob0,b0,b1),out$par))
}

# 21 Feb 2013, multi.expit function
multi.expit <- function(z1,z2) {
 p1 <- exp(z1)/(1+exp(z1)+exp(z2))
 return(p1)
}

# 13 July 2013, geometric mean function
geo.mean <- function(x,zeros=0,na.rm=FALSE) {
 # the user can arbitrarily set the value for a "zero"
 if(na.rm) x <- x[!is.na(x)]
 if(zeros != 0) x[x==0] <- zeros

 n <- length(x)
 out <- prod(x)^(1/n)
 return(out)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

#--- 31 July 2013, multi-logit coefficient estimation
# ...this is NOT working...why?
#multi.logit.coef <- function(y,x) {
# z   <- cbind(logit(y))
# out <- solve(x) %*% z
# return(out)
#}

#y <- c(0.1,0.2,0.8)
#x <- matrix(data=c(1, -1000,  5,
#                   1, -5000, 10,
#                   1, -5000, 15),nrow=3,ncol=3,byrow=TRUE)

#solve(x) %*% x

#out <- multi.logit.coef(y=y,x=x)
#-- check results
#temp <- x %*% out
#expit(sum(x[1,]*out))
#expit(sum(x[2,]*out))
#expit(sum(x[3,]*out))
#out <- lm(logit(y) ~ x[,2] + x[,3])
#fitted(out)
#summary(out)



#out <- logit.coef(p.low=0.1,p.high=0.9,x.low=5,x.high=10)
#expit(1*out$b0 + 5*out$b1)
#expit(1*out$b0 + 10*out$b1)

#--- Linear regression coeffs
lin.reg.coef <- function(x1,x2,y1,y2) {
  b0 <- (y2*x1-y1*x2)/(x1-x2)
  b1 <- (y1-b0)/x1
  out <- list(b0=b0,b1=b1)
  return(out)
}

#y1 <- log(600)
#y2 <- log(2000)
#x1 <- 65
#x2 <- 90
#out <- lin.reg.coef(x1=x1,x2=x2,y1=y1,y2=y2)
#exp(out$b0+out$b1*x1)
#exp(out$b0+out$b1*x2)





#---- Tide Variable extraction
tide.var.extraction <- function(tide.file,survey.file) {
    tide.file$Date    <- as.Date(tide.file$Date,format="%m/%d/%Y")
    tide.file$Date.Time <- paste(tide.file$Date,tide.file$TimeStart)

    fish <- survey.file
    fish.Date.Time <- paste(fish$Date,fish$TimeStart)

   #Match the records in Tide file with those in Survey file
   x <- match(tide.file$Date.Time,fish.Date.Time)

   if(sum(is.na(x)) > 0) {
	 # Remove the records without matches:
     tide.file <- tide.file[!is.na(x),]

		# Redefine x based on the updated tide.file:
		x <- match(tide.file$Date.Time,fish.Date.Time)
   }

   # if(sum(is.na(x)) > 0) {
     # cat("Some records in Tide not found in Fish File",sum(is.na(x)),"\n")
     # print(tide.file[is.na(x),])
     # stop()
   # }

   # fish$TideStage[x]           <- tide.file$TideStage
   # fish$HighType[x]            <- (tide.file$"HighType")
   # fish$"Time-to-High-Min"[x]  <- tide.file$"Time.to.High.Min"
   # fish$LowType[x]             <- (tide.file$"LowType")
   # fish$"Time-to-Low-Min"[x]   <- tide.file$"Time.to.Low.Min"
   # fish$"TideVelocity"[x]        <- tide.file$"TideVelocity"
   # fish$EbbType[x]          <- (tide.file$"Ebb.Type")
   # fish$"Time-to-Ebb-Min"[x]   <- tide.file$"Time.to.Ebb.Min"
   # fish$FloodType[x]           <- (tide.file$"FloodType")
   # fish$"Time-to-Flood-Min"[x] <- tide.file$"Time.to.Flood.Min"
   # fish$"Time-to-Slack-Min"[x] <- tide.file$"Time.to.Slack.Min"

   fish[x,"TideStage"]           <- tide.file$TideStage
   fish[x,"HighType"]           <- (tide.file$"HighType")
   fish[x,"Time-to-High-Min"]  <- tide.file$"Time.to.High.Min"
   fish[x,"LowType"]             <- (tide.file$"LowType")
   fish[x,"Time-to-Low-Min"]   <- tide.file$"Time.to.Low.Min"
   fish[x,"TideVelocity"]        <- tide.file$"TideVelocity"
   fish[x,"EbbType"]         <- (tide.file$"Ebb.Type")
   fish[x,"Time-to-Ebb-Min"]   <- tide.file$"Time.to.Ebb.Min"
   fish[x,"FloodType"]           <- (tide.file$"FloodType")
   fish[x,"Time-to-Flood-Min"] <- tide.file$"Time.to.Flood.Min"
   fish[x,"Time-to-Slack-Min"] <- tide.file$"Time.to.Slack.Min"

   return(fish)
}




############################################################
## Functions added by L. Mitchell (2013 - 2014):


## Define a sum function that automatically uses na.rm=TRUE,
## but returns NA (instead of 0) if x contains all NA's:
sum.na <- function(x) {
	out <- NA
	if(!all(is.na(x))) {out <- sum(x,na.rm=TRUE)}
	return(out)
}


## Define a median function that automatically uses na.rm=TRUE,
## but returns NA if x contains all NA's. This function is redundant, run
## median(rep(NA,5),na.rm=TRUE), for example, but keeping it here for
## consistency:
median.na <- function(x) {
	out <- NA
	if(!all(is.na(x))) {out <- median(x,na.rm=TRUE)}
	return(out)
}


## Combine multiple tows at a given Date-Station into one record:
combine_tows <- function(catch_df_Old, sum_col, mean_col) {
	cat("Combining tows with the same date-station values ...\n")
  cat("Summing over the following columns:\n")
  print(sum_col)
	cat("Averaging over the following columns:\n")
  print(mean_col)

	## Create new data frame with only one instance of each date-station value.
	## Keep the first record from that date-station.
	catch_df_Old_split <- split(catch_df_Old, catch_df_Old$DS)
	catch.New <- do.call("rbind", lapply(catch_df_Old_split,
		function(x) {
			ret <- x[1, ]

			for(col_name in sum_col) {
				ret[ ,col_name] <- sum(x[ ,col_name])
			}

      for(col_name in mean_col) {
        if(all(is.na(x[ ,col_name]))) {
          ret[ ,col_name] <- NA
        } else {
          ret[ ,col_name] <- mean(x[ ,col_name], na.rm=TRUE)
        }
      }

			return(ret)
		})
	)

	# ind_1 <- match(uniq.DS,catch_df_Old$DS)	 # finds FIRST match of arg1 IN arg2
	# catch.New <- catch_df_Old[ind_1,]  	   	 # retrieve 1st record of duplicates
	# rownames(catch.New) <- catch.New$DS	     # rownames for convenience

	# ## Sum delta smelt catches and volume wrt date-station:
	# for(i in c("delta.smelt","delta.smelt.age0","delta.smelt.age1","Volume")) {
		# sum_list <- unlist(lapply(split(catch_df_Old[,i],catch_df_Old$DS), sum))
		# catch.New[names(sum_list),i] <- sum_list
	# }


	## Create and fill-in individual tow columns:
	cat("\nSaving individual tow columns in catch data frame ...\n\n")

	## catch_df_Old MUST be sorted by DS. If Tow numbers are available, sort
	## by those too. Create column "new.tow" in case Tow numbers aren't
	## available or aren't consecutive:
	catch_df_Old <- catch_df_Old[order(catch_df_Old$DS, catch_df_Old$Tow), ]
	catch_df_Old$new.tow <- unlist(lapply(split(catch_df_Old$DS,catch_df_Old$DS),
                          function(x) { seq(length(x))} ))
	for(i in 1:max(catch_df_Old$new.tow)) {
    for(col_name in c(sum_col, mean_col)) {
      new_col_name <- paste0(col_name,".tow",i)

      # Want non-matches to be NA:
      catch.New[ ,new_col_name] <- NA

      # Get all ith tow records:
      tow_df <- subset(catch_df_Old, new.tow == i)
      catch.New[tow_df$DS, new_col_name] <- tow_df[ ,col_name]
    }

		# name1 <- paste0("delta.smelt.tow",i)
		# name2 <- paste0("delta.smelt.age0.tow",i)
		# name3 <- paste0("delta.smelt.age1.tow",i)
		# name4 <- paste0("Volume.tow",i)
		# name5 <- paste0("EstimatedTowDepth_ft.tow",i)
		# name6 <- paste0("Inland_silverside.tow",i)
		# name7 <- paste0("Striped_bass_age0.tow",i)
		# name8 <- paste0("Striped_bass_age1_plus.tow",i)
		# name9 <- paste0("Striped_bass_all.tow",i)
		# name10 <- paste0("Longfin_Smelt.tow",i)

		# # Want non-matches to be NA:
		# catch.New[,c(name1,name2,name3)] <- NA

		# # Get all ith tow records:
		# tow_df <- subset(catch_df_Old, new.tow==i)
		# catch.New[tow_df$DS, name1] <- tow_df$delta.smelt
		# catch.New[tow_df$DS, name2] <- tow_df$delta.smelt.age0
		# catch.New[tow_df$DS, name3] <- tow_df$delta.smelt.age1
		# catch.New[tow_df$DS, name4] <- tow_df$Volume
		# catch.New[tow_df$DS, name5] <- tow_df$EstimatedTowDepth_ft
		# catch.New[tow_df$DS, name6] <- tow_df$Inland_silverside
		# catch.New[tow_df$DS, name7] <- tow_df$Striped_bass_age0
		# catch.New[tow_df$DS, name8] <- tow_df$Striped_bass_age1_plus
		# catch.New[tow_df$DS, name9] <- tow_df$Striped_bass_all
		# catch.New[tow_df$DS, name10] <- tow_df$Longfin_Smelt

	}

	catch.New[ ,c("Tow")] <- NA		# Leave this column(s) for consistency
	return(catch.New)
}




get_invalid_boolean <- function(x, additional_invalid_list=list()) {
	## Element "%in%" can contain a vector, but the others
	## should be scalars. This doesn't check that, though.
	## Automatically check for NAs.

	z <- is.na(x)
	if("%in%" %in% names(additional_invalid_list)) {
		z <- z | (x %in% additional_invalid_list[["%in%"]])
	}

	if("<" %in% names(additional_invalid_list)) {
		z <- z | (x < additional_invalid_list[["<"]])
	}

	if(">" %in% names(additional_invalid_list)) {
		z <- z | (x > additional_invalid_list[[">"]])
	}

	if("<=" %in% names(additional_invalid_list)) {
		z <- z | (x <= additional_invalid_list[["<="]])
	}

	if(">=" %in% names(additional_invalid_list)) {
		z <- z | (x >= additional_invalid_list[[">="]])
	}

	return(z)
}


## Impute by mean specified mean value.
doImpute <- function(dataframe, impute_var, split_vec, additional_invalid_list=list(),
											verbose=TRUE) {
	## dataframe: data frame with values that need to be imputed. A filled-in version of
	##		dataframe gets returned.
	## impute_var: string giving the name of the column where imputation is to be carried out
	##								(only one column at a time)
	## split_vec: string vector of column names. Split dataframe by these columns joined
	## 								together, then calculate the mean for each unique combination.
	## additional_invalid_list: vector giving other values (and ranges of values) that should be
	##								replaced as if they were NA's. NA's are automatically imputed.
	##								See the function get_invalid_boolean().
	##
	## Note that this imputes using MEAN values.
	## Returns a copy of dataframe with imputed values.

	x <- dataframe[ ,impute_var]
	invalid_index <- get_invalid_boolean(x, additional_invalid_list) #(is.na(x) | x %in% additional_invalid_list)
	orig_invalid_num <- sum(invalid_index)

	if(orig_invalid_num > 0) {
		dataframe$split_col <- ""
		for(i in seq(length(split_vec))) {
			dataframe$split_col <- paste(dataframe$split_col, dataframe[ ,split_vec[i]], sep="")
		}

		invalid_data <- dataframe[invalid_index, ]
		valid_data <- dataframe[!invalid_index, ]
		valid_data$tmp <- valid_data[ ,impute_var]
		if(nrow(valid_data) == 0) {
			if(verbose) cat(impute_var,"\n")
			cat("No data to use for imputation!\n")
			return(dataframe)
		}

		means <- aggregate(tmp ~ split_col, data=valid_data, FUN=mean)
		row.names(means) <- means$split_col

		if(verbose) cat(impute_var,"\n")
		cat(" << Impute by",paste(split_vec, collapse="-"),">>\n")
		cat(" Initial number of invalid values: ", orig_invalid_num, "\n")
		cat(" Unique initial invalid values: ", unique(invalid_data[ ,impute_var]), "\n")
		for(i in which(invalid_index)) {
			this_record_id <- dataframe[i,"split_col"]
			replacement_val <- means[this_record_id,"tmp"]

			if(!is.null(replacement_val) & !is.na(replacement_val)) {
				dataframe[i,impute_var] <- replacement_val
			}
		}

		new_x <- dataframe[ ,impute_var]
		new_invalid_num <- sum(get_invalid_boolean(new_x, additional_invalid_list))
		cat(" Current number of invalid values: ", new_invalid_num, "\n")
		cat(" Number filled in: ", orig_invalid_num - new_invalid_num, "\n\n")
	}

	dataframe$split_col <- NULL
	return(dataframe)
}




## Impute fields in the catch dataframe using mean values.
impute_catch_fields <- function(dataframe, impute_var, additional_invalid_list=list()) {
	## dataframe: data frame with values that need to be imputed. A filled-in version of
	##		dataframe gets returned.
	## impute_var: name of column to be imputed over, as a string (only one column
	##		at a time)
	## additional_invalid_list: numeric vector giving other values that should be
	##		replaced as if they were NA's.

	## This function is like doImpute(), except it automatically tries some pre-defined
	## split_vec combinations.
	##
	## Impute with mean values. Try this order:
	## 1) Date-Station (if there are multiple tows)
	## 2) Date-SubRegion
	## 3) Date-Region
	## 4) Year-Month-SubRegion
	## 5) Year-Month-Region
	## 6) Month-Region (historic month-region average across years; 2nd to last resort)
	## 7) Year-Region (across all months in that year; very last resort)


	## Temporarily rename column for convenience. Hopefully it won't overwrite an existing column.
	dataframe$my_var <- dataframe[ ,impute_var]


	#### Initial check for invalid values:
	## All invalid (no data to impute with):
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	if(all(invalid_index)) {
		cat(impute_var,"- ALL INVALID VALUES. Cannot Impute. \n\n")
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## No imputation needed:
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Date-Station mean:
	dataframe <- doImpute(dataframe, impute_var, c("Date","Station"), additional_invalid_list,
													verbose=TRUE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Date-SubRegion mean:
	dataframe <- doImpute(dataframe, impute_var, c("Date","SubRegion"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Date-Region mean:
	dataframe <- doImpute(dataframe, impute_var, c("Date","Region"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Year-Month-SubRegion mean:
	dataframe <- doImpute(dataframe, impute_var, c("Year","Month","SubRegion"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Year-Month-Region mean:
	dataframe <- doImpute(dataframe, impute_var, c("Year","Month","Region"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Month-Region mean:
	dataframe <- doImpute(dataframe, impute_var, c("Month","Region"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	## Impute with Year-Region mean:
	dataframe <- doImpute(dataframe, impute_var, c("Year","Region"), additional_invalid_list,
													verbose=FALSE)
	invalid_index <- get_invalid_boolean(dataframe$my_var, additional_invalid_list) # with(dataframe, is.na(my_var) | my_var %in% additional_invalid_list)
	invalid_num <- sum(invalid_index)
	if(invalid_num == 0) {
		dataframe$my_var <- NULL
		return(dataframe)
	}

	dataframe$my_var <- NULL
	return(dataframe)
}







## Create a linear model for turbidity as a function of secchi for imputation
## (especially in the case of 20mm where lots of turbidity values are missing
## over long periods of time).
impute_turbidity_with_secchi <- function(dataframe, additional_invalid_list=list(), show_plot=TRUE) {
	invalid_index <- get_invalid_boolean(dataframe$Turbidity, additional_invalid_list)
	orig_invalid_num <- sum(invalid_index)

	if(orig_invalid_num > 0) {
		invalid_data <- dataframe[invalid_index, ]
		valid_data <- dataframe[!invalid_index, ]
		if(nrow(valid_data) == 0) {
			if(verbose) {
				cat("Turbidity\n")
			}
			cat("No data to use for imputation!\n")
			return(dataframe)
		}

		#means <- aggregate(Turbidity ~ Secchi, data=valid_data, FUN=mean)
		#my_linear_approx_fun <- approxfun(x=means$Secchi, y=means$Turbidity)

		#means <- filter(valid_data$Turbidity, c(1,1,1)/3, sides=2)
		#my_linear_approx_fun <- approxfun(x=valid_data$Secchi, y=means)

		# my_loess <- loess(Turbidity ~ Secchi, data=valid_data, span=0.1)
		my_linearized_fit <- lm(log(Turbidity) ~ log(Secchi), data=valid_data)


		if(verbose) cat("Turbidity\n")
		cat(" << Impute using linear secchi model >>\n")
		cat(" Initial number of invalid values: ", orig_invalid_num, "\n")
		cat(" Unique initial invalid values: ", unique(invalid_data[ ,"Turbidity"]), "\n")
		if(show_plot) {
			dev.new()
			xmax <- max(valid_data$Secchi, na.rm=TRUE) * 1.2
			ymax <- max(valid_data$Turbidity, na.rm=TRUE) * 1.2
			plot(valid_data$Secchi, valid_data$Turbidity, xlim=c(0,xmax), ylim=c(0,ymax))
		}
		for(i in which(invalid_index)) {
			this_secchi <- dataframe[i,"Secchi"]
			#replacement_val <- my_linear_approx_fun(this_secchi)
			# replacement_val <- exp(predict(my_loess, this_secchi))
			replacement_val <- exp(predict(my_linearized_fit, data.frame(Secchi=this_secchi)))

			if(!is.na(replacement_val)) {
				dataframe[i,"Turbidity"] <- replacement_val
				if(show_plot) {	points(this_secchi, replacement_val, col="red", pch=0) }
			}
		}

		new_invalid_index <- get_invalid_boolean(dataframe$Turbidity, additional_invalid_list)
		new_invalid_num <- sum(new_invalid_index)
		cat(" Current number of invalid values: ", new_invalid_num, "\n")
		cat(" Number filled in: ", orig_invalid_num - new_invalid_num, "\n\n")

	}

	return(dataframe)
}






## Create a linear model for CondBott as a function of CondSurf for imputation.
impute_condbott_with_condsurf <- function(dataframe, additional_invalid_list=list(), show_plot=TRUE) {
	invalid_index <- get_invalid_boolean(dataframe$CondBott, additional_invalid_list)
	orig_invalid_num <- sum(invalid_index)

	if(orig_invalid_num > 0) {
		invalid_data <- dataframe[invalid_index, ]
		valid_data <- dataframe[!invalid_index, ]
		if(nrow(valid_data) == 0) {
			if(verbose) cat("CondBott\n")
			cat("No data to use for imputation!\n")
			return(dataframe)
		}

		lm_fit <- lm(CondBott ~ CondSurf, data=valid_data)

		if(verbose) { cat("CondBott\n") }
		cat(" << Impute using linear CondSurf model >>\n")
		cat(" Initial number of invalid values: ", orig_invalid_num, "\n")
		cat(" Unique initial invalid values: ", unique(invalid_data[ ,"CondBott"]), "\n")
		if(show_plot) {
			dev.new()
			lims <- max(c(valid_data$CondSurf, valid_data$CondBott), na.rm=TRUE) * 1.2
			plot(valid_data$CondSurf, valid_data$CondBott, xlim=c(0,lims), ylim=c(0,lims))
		}
		for(i in which(invalid_index)) {
			this_condsurf <- dataframe[i,"CondSurf"]
			replacement_val <- predict(lm_fit, data.frame(CondSurf=this_condsurf))

			if(!is.na(replacement_val)) {
				dataframe[i,"CondBott"] <- replacement_val
				if(show_plot) {	points(this_condsurf, replacement_val, col="red", pch=0) }
			}
		}

		new_invalid_index <- get_invalid_boolean(dataframe$CondBott, additional_invalid_list)
		new_invalid_num <- sum(new_invalid_index)
		cat(" Current number of invalid values: ", new_invalid_num, "\n")
		cat(" Number filled in: ", orig_invalid_num - new_invalid_num, "\n\n")
	}

	return(dataframe)
}








## Impute missing ages using average proportion of age 0's, trying all the
## averages defined by the function impute_catch_fields.
## For FMWT also apply Dave Contreras' percent estimates (23 Aug 2012).
## This should be run after the catch and length files have been merged.
impute_age <- function(dataframe) {

	subdata <- subset(dataframe, delta.smelt > 0)
	subdata$prop_age0 <- with(subdata, delta.smelt.age0/delta.smelt)

	subdata <- impute_catch_fields(subdata, impute_var="prop_age0", additional_invalid_list=list())
	missing_rownames <- row.names(dataframe)[is.na(dataframe$delta.smelt.age0)]

	if(survey == "FMWT") {	# survey used globally
		for(rn in missing_rownames) {
			this_month <- dataframe[rn,"Month"]

			## Use Dave Contreras' estimates when there are no data-based estimates.
			if(!rn %in% row.names(subdata)) {
				if(this_month %in% month.name[1:5]) {
					use_prop_age0 <- 0
				} else if(this_month %in% month.name[9]) {
					use_prop_age0 <- 0.9
				} else if(this_month %in% month.name[10:12]) {
					use_prop_age0 <- 1
				}
			} else {
				use_prop_age0 <- subdata[rn,"prop_age0"]
			}

			## Also use Dave Contreras' estimates when the empirical estimates in Sept-Nov are < 0.9:
			if(this_month %in% month.name[9:11] && use_prop_age0 < 0.9) {
				use_prop_age0 <- 0.9
			}

			dataframe[rn,"delta.smelt.age0"] <- round(use_prop_age0 * dataframe[rn,"delta.smelt"])
			dataframe[rn,"delta.smelt.age1"] <- dataframe[rn,"delta.smelt"] -
																											dataframe[rn,"delta.smelt.age0"]
		}
	} else {	# the other surveys
		if(any(!missing_rownames %in% row.names(subdata))) {
			cat("Error in impute_age\n\n")
			return(c())
		}

		dataframe[missing_rownames,"delta.smelt.age0"] <-
			round(subdata[missing_rownames,"prop_age0"] * dataframe[missing_rownames,"delta.smelt"])
		dataframe[missing_rownames,"delta.smelt.age1"] <-
			dataframe[missing_rownames,"delta.smelt"] - dataframe[missing_rownames,"delta.smelt.age0"]

		# for(rn in missing_rownames) {
			# use_prop_age0 <- subdata[rn,"prop_age0"]

			# dataframe[rn,"delta.smelt.age0"] <- round(use_prop_age0 * dataframe[rn,"delta.smelt"])
			# dataframe[rn,"delta.smelt.age1"] <- dataframe[rn,"delta.smelt"] -
																											# dataframe[rn,"delta.smelt.age0"]
		# }
	}

	return(dataframe)
}






## Impute missing lengths using average lengths by Year-Month or by
## Month if no Year-Month average is available:
## This should be run after the catch and length files have been merged.

impute_length <- function(dataframe) {
	## Do the imputations in a temporary data frame that only contains records for which
	## there were non-zero catches of delta smelt, then plug those values back in to the
	## main data frame as necessary.
	## Note that if delta.smelt.age0 is 0, Age0_L_bar SHOULD STAY NA instead of being imputed.
	## Similarly for delta.smelta.age1.

	## If there are no missing values to be imputed, don't display any messages.

	## Age 0:
	subdata_0 <- subset(dataframe, delta.smelt.age0 > 0)
	subdata_0 <- impute_catch_fields(subdata_0, impute_var="Age0_L_bar")
	missing_bool <- with(dataframe, delta.smelt.age0 > 0 & is.na(Age0_L_bar))
	missing_rownames <- row.names(dataframe)[missing_bool]

	if(sum(missing_bool) > 0) {
		if(any(!missing_rownames %in% row.names(subdata_0))) {
			cat("Error in impute_length\n\n")
			return(c())
		}

		dataframe[missing_rownames,"Age0_L_bar"] <- subdata_0[missing_rownames,"Age0_L_bar"]
		cat("Imputed", sum(missing_bool), " missing values of Age0_L_bar.\n\n")
	}


	## Do the same thing for Age 1:
	subdata_1 <- subset(dataframe, delta.smelt.age1 > 0)
	subdata_1 <- impute_catch_fields(subdata_1, impute_var="Age1_L_bar")
	missing_bool <- with(dataframe, delta.smelt.age1 > 0 & is.na(Age1_L_bar))
	missing_rownames <- row.names(dataframe)[missing_bool]

	if(sum(missing_bool) > 0) {
		if(any(!missing_rownames %in% row.names(subdata_1))) {
			cat("Error in impute_length\n\n")
			return(c())
		}

		dataframe[missing_rownames,"Age1_L_bar"] <- subdata_1[missing_rownames,"Age1_L_bar"]
		cat("Imputed", sum(missing_bool), " missing values of Age1_L_bar.\n\n")
	}

	return(dataframe)
}



## This takes a numeric vector and returns a string showing the elements of the
## 	vector in a compact way, e.g., c(1,2,3,7,9,10) becomes "1:3, 7, 9:10". It was
## 	designed for integer vectors, but it will technically work for non-integer
## 	vectors as well.
compact_vec <- function(x) {
	if(length(x) %in% c(0,1)) return(x)

	x <- sort(unique(x))
	n <- length(x)
	start_of_run <- 1 - c(0, diff(x) == 1)
	#viewDF <- data.frame(x, start_of_run)

	ind.start <- which(start_of_run == 1)
	n1 <- length(ind.start)
	ind.stop <- c(ind.start[2:n1]-1,n)
	DF <- data.frame("Ind.start"=ind.start, "Ind.stop"=ind.stop)

	vec <- seq(nrow(DF))
	clean_vec <- lapply(vec, function(i) {
				ind1 <- DF[i,1]
				ind2 <- DF[i,2]
				if(ind1 == ind2) {
					return(x[ind1])
				} else {
					return(paste(x[ind1],x[ind2],sep=":"))
				}	})

	return(paste(clean_vec, sep="",collapse=", "))
}




## Customized function for extracting parts of string vectors.
## Example: extractName("2007.January",1) will return "2007".
## This is useful for breaking down list names after using split().
extractName <- function(strVec, keep, sep=".") {
	temp <- sapply(strVec, function(xx) {
		          paste(unlist(strsplit(xx,sep,fixed=TRUE))[keep],collapse=sep)	})
	return(temp)
}




## Split, lapply, then unlist. Has the same "drop" default as split().
## Can also be used in cases where splitting isn't needed.
ULS <- function(useData, splitBy, FUN, unList=TRUE, drop=FALSE, ...) {
	X <- split(useData, splitBy, drop)
	if(unList) {
		unlist(lapply(X, FUN, ...))
	} else {
		lapply(X, FUN, ...)
	}
}



## Split a data frame by Year-Month-Region, Year-Month, or Year, calculate
## the specified averages of the specified fields, then put everything back in a
## new data frame containing all 12 months, all 4 regions, and every year in
## the RANGE of the original data frame. Use NA whenever there are no data
## available to calculate an average (instead of NaN).
## 12/3/14 LM: I added "sum" as an additional option for the argument "average".
##   Should probably rename this function but won't right now.
split.average <- function(useDF, splitBy, average, fields) {
	years <- with(useDF, seq(min(Year),max(Year),1))

	## Split the original data frame by YMR, YM, or Y:
	if(splitBy == "ymr") {
		useList <- with(useDF, split(useDF, list(Year,Month,Region), drop=TRUE))
		newDF <- expand.grid(Month=month.name,Year=years,Region=region.list)
		rownames(newDF) <- with(newDF, paste(Year,Month,Region,sep="."))
	} else if(splitBy == "ym") {
		useList <- with(useDF, split(useDF, list(Year,Month), drop=TRUE))
		newDF <- expand.grid(Month=month.name,Year=years)
		rownames(newDF) <- with(newDF, paste(Year,Month,sep="."))
	} else if(splitBy == "y") {
		useList <- with(useDF, split(useDF, list(Year), drop=TRUE))
		newDF <- expand.grid(Year=years)
		rownames(newDF) <- newDF$Year
	}

	## Calculate the average of the specified fields and save in newDF:
	for(sp in fields) {
	    # Make sure entry is NA if there are no data with which to calculate average.
		# It is not necessary to do this first, but leave this for clarity.
		newDF[,sp] <- NA

		temp <- sapply(useList, function(x) {
					allNA <- all(is.na(x[,sp]))
					if(allNA) {
						return(NA)
					} else {
						switch(average,
							"mean" = mean(x[,sp], na.rm=TRUE),
							"median" = median(x[,sp], na.rm=TRUE),
							"geo" = geo.mean(x[,sp],zeros=0.01,na.rm=TRUE),
							"log.geo" = log(geo.mean(x[,sp],zeros=0.01,na.rm=TRUE)),
							"sum" = sum(x[,sp]))
					}
		})
		newDF[names(temp),sp] <- temp
	}

	## Include a pseudo Date column (Month/1/Year) for easy plotting:
	if(splitBy %in% c("ymr","ym")) {
		newDF$pDate <- with(newDF, paste(match(Month,month.name),1,Year,sep="/"))
		newDF$pDate <- as.Date(newDF$pDate,format="%m/%d/%Y")
		newDF$pDate <- as.numeric(newDF$pDate)
	}

	return(newDF)
}





myReshape <- function(dat, row_var, col_var, fill_var) {
	# This uses reshape() and assumes going from long to wide, one variable only.
	# This was designed to create Year x Month matrices that could be dropped directly
	#   into a dat file, using, for example write.table().
	# After using this function, double check row/column names and/or retrieve things by
	#   row and columns names in order to ensure that they are in the correct order.

	tmp <- subset(dat, select=c(row_var, col_var, fill_var))
	tmp <- reshape(tmp, v.names=fill_var, idvar=row_var, timevar=col_var, direction="wide")
	rownames(tmp) <- as.character(tmp[ ,row_var])
	colnames(tmp) <- gsub(paste0(fill_var,"."), "", colnames(tmp))

	tmp <- tmp[ ,names(tmp) != row_var]
	return(tmp)
}






interp.fun <- function(position,dataset,var.name,window.size=5,verbose=FALSE) {
	## Simplistic interpolation of missing turbidity values.
	## Written by K. Newman.

	n <- dim(dataset)[1]
	no.luck <- TRUE
	block <- window.size
	while(no.luck) {
		if((position-block) <= 0 |  (position+block) > n) {
			print("Warning: ran out of data")
		}
		available <-  dataset[[var.name]][(position-block):(position+block)]
		imputed   <-  mean(available,na.rm=TRUE)
		if(is.na(imputed)) {
			block <- block+window.size
		} else {
			no.luck <- FALSE
			if(verbose) {
				cat("block=",block,"---imputed",imputed,"\n")
			}
		}
	}
	return(imputed)
}



